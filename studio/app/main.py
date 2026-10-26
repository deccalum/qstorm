from __future__ import annotations

import asyncio
import math
from collections.abc import Coroutine
from dataclasses import dataclass
from typing import Any

import js # pyright: ignore[reportMissingImports]
import numpy as np
from pyodide.ffi import create_proxy, to_js

import component
import mixer_fader
from component import Component
from library import Summary, Table, Value, browse, extraction_index, integer, present, real, source_csv
from mixer_fader import FULL, channels
from plot import Mark, bars, fit, ink, marks
from relay import cli, level, mixer, pause, play, read, read_roots, stft, stop, waveform
from view import View, cells, definitions, fields, markers
from waveform import (Atomization, Overlay, Overview, Scene, Spectrum, Stft, Ticks, Waveform, Wave, bounds, height,
                      millisecond, palette)


ZOOM = 2
WHEEL = 0.002
LINE_PX = 16
DRAG_PX = 3
FOLD_PX = 8
MS_PER_MINUTE = 60000
WINDOW_MARGIN = 0.25
MINIMUM_FRAMES = 16
NEIGHBOURS = 6
SHAPE_FLOOR_DB = 6.0
SHAPE_HEADROOM_DB = 3.0

@dataclass
class Select:
    ms: float
    x: float
    moved: bool

@dataclass
class Pan:
    x: float

@dataclass
class Drag:
    window: int
    column: str
    moved: bool

@dataclass
class Seek:
    ms: float

def element(identifier: str) -> Any:
    found = js.document.getElementById(identifier)
    assert found is not None
    return found
def wrap(value: int, count: int) -> int:
    return value % count
def every(selector: str) -> list[Any]:
    return list(js.document.querySelectorAll(selector))
def datum(node: Any, name: str) -> str | None:
    value = node.getAttribute(f"data-{name}")
    return value if isinstance(value, str) else None

Gesture = Select | Pan | Drag | Seek
class Studio:
    def __init__(self, sources: Table, factors: Table, presets: Table) -> None:
        self.sources = sources
        self.factors = factors
        self.presets = presets
        self.wave = Wave(element("wave"), element("window-layer"))
        self.wave.loaded = self.redraw
        self.spectrum = Spectrum(element("stft"))
        self.ticks = Ticks(element("ticks"))
        self.overview = Overview(element("overview"))
        self.ink = ink()
        self.second = True
        self.source_index: int | None = None
        self.source_id = 0
        self.roots: Table | None = None
        self.summary: dict[int, Summary] = {}
        self.deleted = False
        self.windowed = False
        self.atomized = False
        self.query = ""
        self.sort = ("id", 1)
        self.expanded = False
        self.scroll_top = 0.0
        self.editing = False
        self.checked: set[tuple[int, int]] = set()
        self.hidden_column: set[str] = set()
        self.source_directory = ""
        self.rescanned: float | None = None
        self.tab = "components"
        self.showing_stft = True
        self.scene: Scene | None = None
        self.report: Component | None = None
        self.shape: Waveform | None = None
        self.busy = ""
        self.gesture: Gesture | None = None
        self.panning = False
        self.play_clock = 0.0
        self.paused_at: float | None = None
        self.master = FULL
        self.balance = 0
        self.channel_count = 0
        self.mono = True
        self.looping = False
        self.metering = False
        self.following = False
        self.tasks: set[asyncio.Task] = set()
        self.proxies: list[object] = []
        self.keys = {
            " ":          (self.play_toggle, ()),
            "Escape":     (self.clear,       ()),
            "ArrowLeft":  (self.step_through, (-1,)),
            "ArrowRight": (self.step_through, (1,)),
            "ArrowUp":    (self.zoom,        (1/ZOOM,)),
            "ArrowDown":  (self.zoom,        (ZOOM,)),
            "Delete":     (self.reject,      ()),
            "-":          (self.zoom,        (ZOOM,)),
            "+":          (self.zoom,        (1/ZOOM,)),
            "0":          (self.center,      ()),
            "w":          (self.zoom_window, ()),
            "l":          (self.loop,        ()),
            "f":          (self.follow,      ()),
            "m":          (self.mono_toggle, ())}
        self.ctrl_keys = {
            "ArrowLeft":  (self.step_source, (-1,)),
            "ArrowRight": (self.step_source, (1,))}
    def spawn(self, coroutine: Coroutine) -> None:
        task = asyncio.ensure_future(coroutine)
        self.tasks.add(task)
        task.add_done_callback(self.tasks.discard)
    def on(self, target: Any, kind: str, handler: object, options: dict | None = None) -> None:
        proxy = create_proxy(handler)
        self.proxies.append(proxy)
        if options is None:
            target.addEventListener(kind, proxy)
            return
        target.addEventListener(kind, proxy, to_js(options, dict_converter=js.Object.fromEntries))
    def listen(self) -> None:
        canvas = self.wave.canvas
        self.on(js.document, "keydown", self.key)
        self.on(js.document, "mousedown", self.press)
        self.on(element("list"), "click", self.chosen)
        self.on(element("analysis"), "scroll", self.analysis_scroll)
        self.on(js.window.matchMedia("(prefers-color-scheme: light)"), "change", self.scheme_change)
        self.on(element("edit-table"), "click", self.edit_toggle)
        self.on(element("edit-done"), "click", self.edit_toggle)
        self.on(element("delete-checked"), "click", lambda event: self.spawn(self.reject_checked()))
        self.on(element("restore-checked"), "click", lambda event: self.spawn(self.restore_checked()))
        self.on(element("columns"), "click", self.columns_toggle)
        self.on(element("column-list"), "click", self.column_chosen)
        self.on(element("export"), "click", self.export_csv)
        self.on(element("rescan"), "click", lambda event: self.spawn(self.rescan()))
        self.on(element("reset"), "click", self.reset_ask)
        self.on(element("reset-cancel"), "click", self.reset_ask)
        self.on(element("reset-run"), "click", lambda event: self.spawn(self.reset_extraction()))
        self.on(js.document.querySelector(".browser"), "click", self.browser_click)
        self.on(element("query"), "input", self.query_change)
        self.on(element("atom-table"), "click", self.component_chosen)
        self.on(js.document.querySelector(".panel.atomization .head"), "click", self.tab_chosen)
        self.on(js.document.querySelector(".layer-key"), "click", self.layer_chosen)
        self.on(element("factor"), "change", self.factor_change)
        self.on(element("factor-fields"), "change", self.factor_edit)
        self.on(element("copy-factor"), "click", lambda event: self.spawn(self.copy_factor()))
        self.on(element("ard"), "change", self.ard_change)
        self.on(element("center"), "click", self.center)
        self.on(element("zoom-out"), "click", lambda event: self.zoom(ZOOM))
        self.on(element("zoom-window"), "click", self.zoom_window)
        self.on(element("zoom-in"), "click", lambda event: self.zoom(1 / ZOOM))
        self.on(element("loop"), "click", self.loop)
        self.on(element("follow"), "click", self.follow)
        self.on(element("compact"), "click", lambda event: self.look("compact"))
        self.on(element("half"), "click", lambda event: self.look("half"))
        self.on(element("logarithmic"), "click", self.decibel)
        self.on(element("show-stft"), "click", self.stft_toggle)
        self.on(element("floor"), "change", self.floor)
        self.on(element("play"), "click", lambda event: self.spawn(self.play_toggle()))
        self.on(element("stop"), "click", self.halt)
        self.on(element("previous"), "click", lambda event: self.step(-1))
        self.on(element("next"), "click", lambda event: self.step(1))
        self.on(element("previous-source"), "click", lambda event: self.spawn(self.step_source(-1)))
        self.on(element("next-source"), "click", lambda event: self.spawn(self.step_source(1)))
        self.on(element("reject"), "click", self.reject)
        self.on(element("markers"), "focusin", self.field_focus)
        self.on(element("markers"), "change", self.field_change)
        self.on(element("extract"), "click", lambda event: self.spawn(self.extract()))
        self.on(element("atomize"), "click", lambda event: self.spawn(self.atomize()))
        self.on(element("analyze"), "click", lambda event: self.spawn(self.analyze()))
        self.on(element("solo"), "click", lambda event: self.spawn(self.solo()))
        self.on(element("master"), "input", self.master_change)
        self.on(element("pan"), "input", self.balance_change)
        self.on(element("mono"), "click", self.mono_toggle)
        self.on(canvas, "wheel", self.wheel, {"passive": False})
        self.on(canvas, "pointerdown", self.pointer_down)
        self.on(canvas, "pointermove", self.pointer_move)
        self.on(canvas, "pointerup", self.pointer_up)
        self.on(self.ticks.canvas, "pointerdown", self.ticks_down)
        self.on(self.ticks.canvas, "pointermove", self.ticks_move)
        self.on(self.ticks.canvas, "pointerup", self.ticks_up)
        self.on(self.overview.canvas, "pointerdown", self.overview_down)
        self.on(self.overview.canvas, "pointermove", self.overview_move)
        self.on(self.overview.canvas, "pointerup", self.overview_up)
        mixer_fader.install(self.on)
    def key(self, event: Any) -> None:
        if event.target.closest("input, select"): return
        action = (self.ctrl_keys if event.ctrlKey else self.keys).get(event.key)
        if action is None: return
        event.preventDefault()
        method, argument = action
        result = method(*argument)
        if isinstance(result, Coroutine): self.spawn(result)
    def press(self, event: Any) -> None:
        if event.target.closest("button"): event.preventDefault()
    def chosen(self, event: Any) -> None:
        button = event.target.closest("button")
        if not button: return
        if datum(button, "sort") is not None:
            name = datum(button, "sort")
            self.sort = (name, -self.sort[1] if self.sort[0] == name else 1)
            self.browser()
            return
        if datum(button, "source") is None: return
        index = int(button.dataset.source)
        if button.classList.contains("window") and self.editing:
            self.check_window(int(button.dataset.window))
            return
        if button.classList.contains("window"):
            self.spawn(self.open_window(index, int(button.dataset.window)))
            return
        self.spawn(self.open(index))
    def check_window(self, window: int) -> None:
        scene = self.scene
        assert scene is not None and scene.roots is not None
        key = extraction_index(scene.roots, scene.roots.row[window])
        if key in self.checked: self.checked.discard(key)
        else: self.checked.add(key)
        self.browser()
    def edit_toggle(self, event: object = None) -> None:
        self.editing = not self.editing
        self.checked.clear()
        if self.editing: self.expanded = True
        self.browser()
    def columns_toggle(self, event: object = None) -> None:
        listed = element("column-list")
        listed.hidden = not listed.hidden
    def column_chosen(self, event: Any) -> None:
        button = event.target.closest("button")
        if not button: return
        name = button.dataset.column
        if name in self.hidden_column: self.hidden_column.discard(name)
        else: self.hidden_column.add(name)
        self.browser()
    def export_csv(self, event: object = None) -> None:
        text = source_csv(self.sources, self.summary, self.order(), self.hidden_column, self.second)
        blob = js.Blob.new(to_js([text]), to_js({"type": "text/csv"}, dict_converter=js.Object.fromEntries))
        address = js.URL.createObjectURL(blob)
        link = js.document.createElement("a")
        link.href = address
        link.download = "source.csv"
        link.click()
        js.URL.revokeObjectURL(address)
    def reset_ask(self, event: object = None) -> None:
        confirm = element("reset-confirm")
        confirm.hidden = not confirm.hidden
        element("reset").hidden = not confirm.hidden
    def scheme_change(self, event: object) -> None:
        self.ink = ink()
        for surface in (self.wave, self.spectrum, self.ticks, self.overview): surface.palette = palette(surface.canvas)
        self.spectrum.painted = None
        self.spectrum.component_painted = None
        self.changed()
    def analysis_scroll(self, event: Any) -> None:
        stage = event.currentTarget
        top = stage.scrollTop
        if abs(top - self.scroll_top) < FOLD_PX: return
        if top < stage.scrollHeight - stage.clientHeight: element("chrome").classList.toggle("hidden", top > self.scroll_top)
        self.scroll_top = top
    def browser_click(self, event: Any) -> None:
        button = event.target.closest("button")
        if not button: return
        view = datum(button, "view")
        unit = datum(button, "unit")
        flag = datum(button, "filter")
        source = datum(button, "source")
        if view is not None:
            self.expanded = view == "expanded"
            self.browser()
            return
        if unit is not None:
            self.second = unit == "sec"
            self.browser()
            return
        if flag == "deleted":
            self.show_deleted()
            return
        if flag == "windowed": self.windowed = not self.windowed
        if flag == "atomized": self.atomized = not self.atomized
        if flag is not None:
            self.browser()
            return
        if source is not None and button.closest(".fade"): self.spawn(self.open(int(source)))
    def query_change(self, event: Any) -> None:
        self.query = str(event.target.value).strip().lower()
        self.browser()
    def factor_change(self, event: object) -> None:
        self.factor()
    def factor_edit(self, event: Any) -> None:
        number = event.target
        self.spawn(self.edit_factor(number.dataset.column, number.valueAsNumber))
    def ard_change(self, event: object) -> None:
        self.ard()
    def decibel(self, event: object = None) -> None:
        self.look("decibel")
    def floor(self, event: Any) -> None:
        self.wave.look.floor = float(event.target.value)
        self.spectrum.floor = self.wave.look.floor
        self.redraw()
    def stft_toggle(self, event: object = None) -> None:
        self.showing_stft = not self.showing_stft
        element("spectrum-plot").hidden = not self.showing_stft
        self.redraw()
    def show_deleted(self, event: object = None) -> None:
        self.deleted = not self.deleted
        if self.scene is not None:
            roots = self.roots
            assert roots is not None
            self.scene.roots = self.kept(roots)
            self.choose(None)
        self.changed()
    def loop(self, event: object = None) -> None:
        self.looping = not self.looping
        self.redraw()
    def follow(self, event: object = None) -> None:
        self.following = not self.following
        self.redraw()
    def field_focus(self, event: Any) -> None:
        if datum(event.target, "field") is None: return
        self.focus_field(int(event.target.dataset.field))
    def field_change(self, event: Any) -> None:
        scene = self.scene
        assert scene is not None and scene.window is not None
        number = event.target
        self.spawn(self.commit(scene.window, number.dataset.column, number.valueAsNumber))
    def master_change(self, event: Any) -> None:
        self.master = int(event.target.value)
        self.send_mixer()
    def balance_change(self, event: Any) -> None:
        self.balance = int(event.target.value)
        self.send_mixer()
    def mono_toggle(self, event: object = None) -> None:
        self.mono = not self.mono
        self.send_mixer()
        self.redraw()
    def send_mixer(self) -> None:
        step = channels(self.balance, self.channel_count)
        self.spawn(mixer(self.master, step, self.mono))
    def status(self, text: str) -> None:
        self.busy = text
        element("status").textContent = text
    def source_cell(self, name: str) -> Value:
        index = self.source_index
        assert index is not None
        row = self.sources.row[index]
        return row[self.sources.name.index(name)]
    def factor_id(self) -> int:
        select = element("factor")
        return self.factors.key_cells(self.factors.row[select.selectedIndex])[0]
    def ard_id(self) -> int:
        select = element("ard")
        return self.presets.key_cells(self.presets.row[select.selectedIndex])[0]
    async def survey(self) -> None:
        every_root = await read("root")
        extractions = await read("extraction")
        atomizations = await read("atomization")
        source_of_extraction = {}
        newest: dict[int, int] = {}
        for row in extractions.row:
            extraction = extractions.key_cells(row)[0]
            source_id = integer(row[extractions.name.index("source_id")])
            source_of_extraction[extraction] = source_id
            newest[source_id] = max(newest.get(source_id, extraction), extraction)
        windows: dict[int, int] = {}
        rejected = every_root.name.index("rejected")
        for row in every_root.row:
            extraction = integer(row[every_root.name.index("extraction")])
            source_id = source_of_extraction[extraction]
            if extraction == newest[source_id] and not row[rejected]: windows[source_id] = windows.get(source_id, 0) + 1
        atomized: dict[int, set[tuple[int, int]]] = {}
        for row in atomizations.row:
            extraction = integer(row[atomizations.name.index("extraction")])
            source_id = source_of_extraction[extraction]
            atomized.setdefault(source_id, set()).add((extraction, integer(row[atomizations.name.index("index")])))
        for row in self.sources.row:
            source_id = self.sources.key_cells(row)[0]
            self.summary[source_id] = Summary(windows.get(source_id, 0), len(atomized.get(source_id, set())))
    async def open(self, index: int) -> None:
        if self.busy: return
        self.halt()
        source_id = self.sources.key_cells(self.sources.row[index])[0]
        self.status(f"{source_id} decoding")
        width = self.wave.canvas.width
        raw = await waveform(source_id, -math.inf, math.inf, width)
        shown = Waveform.parse(raw)
        whole = Waveform.parse(await waveform(source_id, -math.inf, math.inf, max(self.overview.canvas.width, 1)))
        roots = await read_roots(source_id)
        self.roots = roots
        self.source_index = index
        self.source_id = source_id
        self.checked.clear()
        self.channel_count = integer(self.source_cell("channels"))
        self.send_mixer()
        view = View(shown.x0, shown.x1, MINIMUM_FRAMES*shown.frame_ms)
        self.scene = Scene(
            view=view, source_id=source_id, waveform=shown,
            shown=(view.x0, view.x1, width),
            roots=self.kept(roots), window=None, field=None,
            selection=None, play_from=None, play_until=None,
            playhead=None, atomization=None, component=None, stft=None,
            reconstruction=None, residual=None, overview=whole, overlay=None)
        self.report = None
        self.shape = None
        for surface in (self.wave, self.spectrum, self.ticks, self.overview): surface.scene = self.scene
        self.status("")
        self.changed()
    async def step_source(self, direction: int) -> None:
        order = self.order()
        if not order: return
        position = order.index(self.source_index) if self.source_index in order else -1
        await self.open(order[wrap(position + direction, len(order))])
    async def open_window(self, index: int, window: int) -> None:
        if index != self.source_index: await self.open(index)
        scene = self.scene
        assert scene is not None and scene.roots is not None
        roots = scene.roots
        self.choose(window)
        x0, x1 = bounds(roots.ms(roots.row[window]))
        scene.view.focus(x0, x1, WINDOW_MARGIN)
        self.changed()
    def root_key(self, window: int) -> tuple[int, int, int]:
        scene = self.scene
        assert scene is not None and scene.roots is not None
        row = scene.roots.row[window]
        extraction = integer(row[scene.roots.name.index("extraction")])
        index = integer(row[scene.roots.name.index("index")])
        revision = integer(row[scene.roots.name.index("revision")])
        return extraction, index, revision
    async def latest_atomization(self, window: int) -> Atomization | None:
        extraction, index, revision = self.root_key(window)
        runs = await read("atomization", extraction)
        made = []
        for row in runs.row:
            same = row[runs.name.index("index")] == index and row[runs.name.index("revision")] == revision
            if same: made.append(integer(row[runs.name.index("id")]))
        if not made: return None
        return await Atomization.fetch(max(made))
    async def load_report(self) -> None:
        scene = self.scene
        self.report = None
        if scene is None or scene.atomization is None or scene.window is None: return
        atomization = scene.atomization
        found = await component.fetch(atomization, self.root_key(scene.window)[0])
        if scene.atomization is atomization: self.report = found
    def component_index(self) -> int | None:
        scene = self.scene
        if scene is None or scene.atomization is None or scene.component is None: return None
        return int(scene.atomization.index[scene.component])
    def component_chosen(self, event: Any) -> None:
        button = event.target.closest("button")
        scene = self.scene
        if not button or scene is None or scene.atomization is None or datum(button, "component") is None: return
        k = int(button.dataset.component)
        position = list(scene.atomization.index).index(k)
        self.halt()
        replay = datum(button, "play") is not None
        scene.component = position if replay else None if scene.component == position else position
        self.place_overlay()
        self.changed()
        if replay: self.spawn(self.play_toggle())
    def tab_chosen(self, event: Any) -> None:
        button = event.target.closest("button")
        if not button or datum(button, "tab") is None: return
        self.tab = button.dataset.tab
        self.atomization()
    def layer_chosen(self, event: Any) -> None:
        button = event.target.closest("button")
        if not button or datum(button, "layer") is None or self.busy: return
        self.spectrum.layer = button.dataset.layer
        self.spawn(self.reload_layer())
    def place_overlay(self) -> None:
        scene = self.scene
        if scene is None: return
        scene.overlay = None
        k = self.component_index()
        if k is None or self.report is None: return
        each = self.report.part[k]
        if each["adsr_x0"] is None or each["spectral_low"] is None: return
        scene.overlay = Overlay(real(each["adsr_x0"]), real(each["adsr_x1"]), real(each["spectral_low"]),
                                real(each["spectral_high"]), f"k{k}", None)
        assert scene.atomization is not None
        self.spawn(self.load_overlay_layer(scene.overlay, scene.atomization.id, k))
    async def load_overlay_layer(self, overlay: Overlay, atomization_id: int, k: int) -> None:
        raw = await cli("render_component", atomization_id, k)
        overlay.layer = Stft.parse(raw)
        self.redraw()
    async def sound(self, offset: float) -> None:
        scene = self.scene
        assert scene is not None
        k = self.component_index()
        if k is not None:
            assert self.report is not None
            await play(integer(self.report.part[k]["audio"]), offset)
            return
        if scene.window is not None:
            roots = scene.roots
            assert roots is not None
            await play(integer(roots.row[scene.window][roots.name.index("audio")]), offset)
            return
        await play(integer(self.source_cell("audio")), offset)
    async def solo(self) -> None:
        if self.scene is None or self.scene.component is None: return
        self.halt()
        await self.play_toggle()
    def changed(self) -> None:
        self.browser()
        self.source()
        self.factor()
        self.ard()
        self.window()
        self.atomization()
        self.redraw()
    async def load_window(self, window: int) -> None:
        scene = self.scene
        assert scene is not None
        extraction, index, revision = self.root_key(window)
        element("stft-key").textContent = f"stft {extraction}/{index} loading"
        raw = await stft(extraction, index, revision)
        atomization = await self.latest_atomization(window)
        if scene.window != window: return
        scene.stft = Stft.parse(raw)
        scene.atomization = atomization
        scene.component = None
        scene.overlay = None
        scene.reconstruction = None
        scene.residual = None
        await self.load_report()
        await self.load_shape(window)
        await self.reload_layer()
        self.changed()
    async def load_shape(self, window: int) -> None:
        scene = self.scene
        assert scene is not None and scene.roots is not None
        x0, x1 = bounds(scene.roots.ms(scene.roots.row[window]))
        width = max(element("shape").getBoundingClientRect().width, 1)
        self.shape = Waveform.parse(await waveform(scene.source_id, x0, x1, round(width)))
    def kept(self, roots: Table) -> Table:
        if self.deleted: return roots
        rejected = roots.name.index("rejected")
        row = []
        for each in roots.row:
            if not each[rejected]: row.append(each)
        return Table(roots.column, row, roots.setting)
    def store(self, roots: Table) -> None:
        self.roots = roots
        scene = self.scene
        assert scene is not None
        scene.roots = self.kept(roots)
    def order(self) -> list[int]:
        found = []
        for i, row in enumerate(self.sources.row):
            source_id = self.sources.key_cells(row)[0]
            counted = self.summary.get(source_id, Summary(0, 0))
            if self.windowed and counted.windows == 0: continue
            if self.atomized and counted.atomized == 0: continue
            if self.query and not any(self.query in str(value).lower() for value in row if isinstance(value, (int, str))): continue
            found.append(i)
        name, direction = self.sort
        def value(i: int) -> tuple[int, Any]:
            source_id = self.sources.key_cells(self.sources.row[i])[0]
            if name == "windows": return (0, self.summary[source_id].windows)
            if name == "atomized": return (0, self.summary[source_id].atomized)
            cell = self.sources.row[i][self.sources.name.index(name)]
            if cell is None: return (1, 0)
            return (0, cell if isinstance(cell, (int, float)) else str(cell).lower())
        found.sort(key=value, reverse=direction < 0)
        return found
    def browser(self) -> None:
        element("browser-compact").hidden = self.expanded
        element("browser-filter").hidden = not self.expanded
        element("browser-expanded").hidden = not self.expanded
        element("edit-table").classList.toggle("active", self.editing)
        element("edit-bar").hidden = not self.editing
        element("checked-count").textContent = f"{len(self.checked)} selected"
        element("rescan-age").textContent = self.rescan_age()
        self.column_menu()
        view = element("browser-view")
        view.dataset.view = "compact" if self.expanded else "expanded"
        view.classList.toggle("expanded", self.expanded)
        element("browser-view-name").textContent = "collapse" if self.expanded else "expand"
        for button in every("[data-unit]"): button.classList.toggle("active", (button.dataset.unit == "sec") == self.second)
        for button in every("[data-filter]"):
            on = {"deleted": self.deleted, "windowed": self.windowed, "atomized": self.atomized}[button.dataset.filter]
            button.classList.toggle("active", on)
        order = self.order()
        position = order.index(self.source_index) + 1 if self.source_index in order else "–"
        element("position").textContent = f"{position} / {len(order)}"
        element("export-count").textContent = f"{len(order)} rows"
        if self.expanded:
            chosen = self.scene.window if self.scene else None
            roots = self.scene.roots if self.scene else None
            browse(element("list"), self.sources, self.summary, order, self.source_index, roots, chosen, self.second, self.sort,
                   self.hidden_column, self.editing, self.checked)
            windows = sum(self.summary[self.sources.key_cells(self.sources.row[i])[0]].windows for i in order)
            atomized = sum(1 for i in order if self.summary[self.sources.key_cells(self.sources.row[i])[0]].atomized)
            element("totals").textContent = f"{len(order)} of {len(self.sources.row)} · {windows} windows · {atomized} atomized"
            return
        self.neighbours(order)
    def column_menu(self) -> None:
        holder = element("column-list")
        holder.replaceChildren()
        total = 0
        for i, column in enumerate(self.sources.column):
            if not self.sources.shown(i, False): continue
            total += 1
            button = js.document.createElement("button")
            button.type = "button"
            button.className = "menu-item"
            button.dataset.column = column.name
            button.classList.toggle("active", column.name not in self.hidden_column)
            name = js.document.createElement("span")
            name.textContent = column.name
            switch = js.document.createElement("span")
            switch.className = "switch"
            button.append(name, switch)
            holder.append(button)
        element("column-count").textContent = f"{total - len(self.hidden_column)} of {total}"
    def rescan_age(self) -> str:
        if self.rescanned is None: return ""
        minutes = (js.Date.now() - self.rescanned) / MS_PER_MINUTE
        return f"{minutes:.0f} min ago"
    def neighbours(self, order: list[int]) -> None:
        position = order.index(self.source_index) if self.source_index in order else -1
        current = element("current-source")
        current.textContent = str(self.source_id) if self.source_index is not None else "pick a recording"
        left = order[max(0, position - NEIGHBOURS):position] if position >= 0 else []
        right = order[position + 1:position + 1 + NEIGHBOURS] if position >= 0 else order[:2*NEIGHBOURS]
        for holder, part in ((element("left-sources"), left), (element("right-sources"), right)):
            holder.replaceChildren()
            for i in part:
                source_id = self.sources.key_cells(self.sources.row[i])[0]
                counted = self.summary[source_id]
                button = js.document.createElement("button")
                button.type = "button"
                mark = js.document.createElement("i")
                mark.className = "atomized" if counted.atomized else "windowed" if counted.windows else ""
                button.append(mark, str(source_id))
                button.dataset.source = str(i)
                holder.append(button)
        stats = element("source-stats")
        stats.replaceChildren()
        if self.source_index is None: return
        counted = self.summary[self.source_id]
        scene = self.scene
        pair = [("windows", str(counted.windows)), ("atomized", str(counted.atomized))]
        if scene is not None: pair.append(("total", self.duration(scene.view.src_x1 - scene.view.src_x0)))
        if scene is not None and scene.window is not None and scene.roots is not None:
            x0, x1 = bounds(scene.roots.ms(scene.roots.row[scene.window]))
            pair.append((f"window {scene.window + 1}", self.duration(x1 - x0)))
        for name, value in pair:
            each = js.document.createElement("span")
            label = js.document.createElement("span")
            label.className = "name"
            label.textContent = " " + name
            shown = js.document.createElement("span")
            shown.className = "value"
            shown.textContent = value
            each.append(shown, label)
            stats.append(each)
    def duration(self, ms: float) -> str:
        return f"{ms / 1000:.1f} s" if self.second else f"{ms:.0f} ms"
    def source(self) -> None:
        if self.source_index is None:
            element("source-id").textContent = ""
            element("source-detail").textContent = ""
            return
        element("source-id").textContent = str(self.source_id)
        element("source-detail").textContent = self.source_detail()
    def source_detail(self) -> str:
        codec = self.source_cell("original_codec")
        if codec is None: return ""
        part = [str(self.source_cell("original_container")), str(codec)]
        depth = self.source_cell("original_bit_depth")
        if depth is not None: part.append(f"{integer(depth)} bit")
        rate = self.source_cell("original_bit_rate")
        if rate is not None: part.append(f"{integer(rate) / 1000:.0f} kbit/s")
        part.append(f"{real(self.source_cell('hz')):g} kHz")
        layout = self.source_cell("original_channel_layout")
        channels = f"{integer(self.source_cell('channels'))} ch"
        part.append(channels if layout is None else f"{channels} {layout}")
        return " · ".join(part)
    def factor(self) -> None:
        select = element("factor")
        chosen = max(select.selectedIndex, 0)
        select.replaceChildren()
        name = self.factors.name.index("name")
        column = self.factors.column[name]
        setting = self.factors.setting[name]
        for row in self.factors.row:
            option = js.document.createElement("option")
            dna_id = self.factors.key_cells(row)[0]
            option.textContent = f"{dna_id} {present(row[name], column, setting, False, False)}"
            select.append(option)
        select.selectedIndex = chosen
        row = self.factors.row[chosen]
        fields(element("factor-fields"), self.factors, row, self.factors.real_name(), None)
    def ard(self) -> None:
        select = element("ard")
        chosen = max(select.selectedIndex, 0)
        select.replaceChildren()
        name = self.presets.name.index("name")
        column = self.presets.column[name]
        setting = self.presets.setting[name]
        for row in self.presets.row:
            option = js.document.createElement("option")
            option.textContent = present(row[name], column, setting, False, False)
            select.append(option)
        select.selectedIndex = chosen
        definitions(element("ard-fields"), cells(self.presets, chosen, False, [""]))
    def window(self) -> None:
        key = element("window-key")
        scene = self.scene
        roots = scene.roots if scene else None
        reject = element("reject")
        if scene is None or roots is None or scene.window is None:
            reject.disabled = True
            reject.classList.remove("active")
            element("window-index").textContent = "– / –" if roots is None else f"– / {len(roots.row)}"
            element("window-range").textContent = ""
            element("markers").replaceChildren()
            fit(element("shape"))
            key.textContent = ""
            return
        window = scene.window
        row = roots.row[window]
        rejected = row[roots.name.index("rejected")]
        reject.disabled = bool(self.busy)
        reject.classList.toggle("active", bool(rejected))
        element("window-index").textContent = f"{window + 1} / {len(roots.row)}"
        x0, x1 = bounds(roots.ms(row))
        element("window-range").textContent = f"{millisecond(x0)} – {millisecond(x1)} ms · {x1 - x0:.0f} ms"
        markers(element("markers"), roots, row, scene.field)
        self.draw_shape(roots, row)
        shown = []
        for name, text in cells(roots, window, True, [""]): shown.append(f"{name} {text}")
        key.textContent = " · ".join(shown)
    def draw_shape(self, roots: Table, row: list[Value]) -> None:
        canvas = element("shape")
        context, ratio = fit(canvas)
        shown = self.shape
        if shown is None: return
        floor = real(row[roots.name.index("floor")])
        top = real(row[roots.name.index("peak_y")])
        y0 = floor - SHAPE_FLOOR_DB
        y1 = top + SHAPE_HEADROOM_DB
        with np.errstate(divide="ignore"):
            level = 20*np.log10(shown.peak())
        share = np.clip((level - y0) / (y1 - y0), 0.02, 0.98)
        bars(context, share, self.ink.faint, (2*ratio, 0, canvas.width - 4*ratio, canvas.height), ratio)
        y = canvas.height*(1 - (floor - y0) / (y1 - y0))
        context.setLineDash(js.Array.new(3*ratio, 3*ratio))
        context.strokeStyle = self.ink.strong
        context.lineWidth = ratio
        context.beginPath()
        context.moveTo(0, y)
        context.lineTo(canvas.width, y)
        context.stroke()
        context.setLineDash(js.Array.new())
        context.fillStyle = self.ink.mid
        context.textBaseline = "bottom"
        context.fillText("floor", canvas.width - 32*ratio, y - 2*ratio)
        x0, x1 = bounds(roots.ms(row))
        found = []
        for column in roots.ms_name():
            ms = real(row[roots.name.index(column)])
            name = column.removesuffix("_x")
            label = f"peak {top:.1f}" if name == "peak" else name
            found.append(Mark((ms - x0) / (x1 - x0), label, name == "peak"))
        marks(context, found, (0, 0, canvas.width, canvas.height), self.ink, ratio)
    def atomization(self) -> None:
        scene = self.scene
        atomization = scene.atomization if scene else None
        found = self.report
        for button in every("[data-tab]"): button.classList.toggle("active", button.dataset.tab == self.tab)
        element("analyze").disabled = atomization is None or bool(self.busy)
        table = element("atom-table")
        if scene is None or atomization is None or found is None or scene.stft is None:
            element("atomization-key").textContent = ""
            table.replaceChildren()
            for name in ("run", "shape-facts", "band-residual", "trace-axis"): element(name).replaceChildren()
            fit(element("trace"))
            self.clear_detail()
            return
        k = self.component_index()
        window = f"window {scene.window + 1}" if scene.window is not None else ""
        element("atomization-key").textContent = f"#{atomization.id} · {window} · K {len(found.order)} of {int(real(found.run['factor_K']))} · floor {atomization.floor:.3g}"
        component.summary(found, scene.stft, self.ink)
        if self.tab == "components": component.components(table, found, scene.stft, k, self.ink)
        if self.tab == "groups": component.groups(table, found, k)
        component.detail(found, k, atomization, scene.stft, self.ink)
    def clear_detail(self) -> None:
        element("detail-title").textContent = "component"
        element("detail-sub").textContent = ""
        element("detail-sections").replaceChildren()
        element("solo").disabled = True
        element("detail-w-range").textContent = ""
        element("detail-h-range").textContent = ""
        fit(element("detail-w"))
        fit(element("detail-h"))
    def redraw(self) -> None:
        for surface in (self.wave, self.spectrum, self.ticks, self.overview): surface.request()
        scene = self.scene
        ready = scene is not None and not self.busy
        element("extract").disabled = not ready
        element("atomize").disabled = scene is None or bool(self.busy) or scene.window is None
        element("loop").classList.toggle("active", self.looping)
        element("follow").classList.toggle("active", self.following)
        element("mono").classList.toggle("active", self.mono)
        element("show-stft").classList.toggle("active", self.showing_stft)
        for button in every("[data-layer]"):
            button.classList.toggle("active", button.dataset.layer == self.spectrum.layer)
            button.disabled = button.dataset.layer != "source" and (scene is None or scene.atomization is None)
        if scene is None: return
        view = scene.view
        frames_per_pixel = self.wave.frames_per_pixel()
        decimals = 2 if frames_per_pixel < 10 else 0
        element("range").textContent = f"{millisecond(view.src_x0)} – {millisecond(view.src_x1)} ms · {frames_per_pixel:.{decimals}f} frames / px"
        element("view-range").textContent = f"{millisecond(view.x0)} – {millisecond(view.x1)} ms"
        look = self.wave.look
        element("compact").classList.toggle("active", look.compact)
        element("half").classList.toggle("active", look.half)
        element("logarithmic").classList.toggle("active", look.decibel)
        shown = self.spectrum.shown(scene)
        if shown is not None:
            layer = {"source": "stft", "reconstruction": "WH", "residual": "V − WH"}[self.spectrum.layer if shown is not scene.stft else "source"]
            element("stft-key").textContent = f"{layer} · {len(shown.band)} bands · hop {shown.step:.1f} ms"
        self.transport()
    def transport(self) -> None:
        scene = self.scene
        assert scene is not None
        playing = scene.playhead is not None and self.paused_at is None
        if scene.playhead is None: scene.play_from, scene.play_until = self.play_span()
        play_button = element("play")
        play_button.textContent = "❚❚" if playing else "▶"
        position = scene.playhead if scene.playhead is not None else scene.play_from
        clock = element("clock")
        clock.replaceChildren()
        clock.append(millisecond(position))
        until = js.document.createElement("span")
        until.textContent = f" / {millisecond(scene.play_until)} ms"
        clock.append(until)
        selection = "none"
        if scene.selection:
            x0, x1 = bounds(scene.selection)
            selection = f"{millisecond(x0)} – {millisecond(x1)} ms"
        element("target").textContent = f"{self.play_target()} · selection {selection}"
    def center(self, event: object = None) -> None:
        if self.scene is None: return
        self.scene.view.whole()
        self.redraw()
    def zoom_window(self, event: object = None) -> None:
        scene = self.scene
        if scene is None or scene.window is None or scene.roots is None: return
        roots = scene.roots
        x0, x1 = bounds(roots.ms(roots.row[scene.window]))
        scene.view.focus(x0, x1, WINDOW_MARGIN)
        self.redraw()
    def look(self, option: str) -> None:
        look = self.wave.look
        setattr(look, option, not getattr(look, option))
        self.redraw()
    def reject(self, event: object = None) -> None:
        if self.scene is None or self.scene.window is None: return
        self.spawn(self.reject_window(self.scene.window))
    async def reject_window(self, window: int) -> None:
        scene = self.scene
        assert scene is not None and scene.roots is not None
        roots = scene.roots
        rejected = roots.row[window][roots.name.index("rejected")]
        await self.commit(window, "rejected", 0 if rejected else 1)
        if not self.deleted: self.choose(None)
        await self.survey()
        self.changed()
    def zoom(self, factor: float) -> None:
        if self.scene is None: return
        view = self.scene.view
        view.zoom(factor, (view.x0 + view.x1) / 2)
        self.redraw()
    def step(self, direction: int) -> None:
        scene = self.scene
        if scene is None or scene.roots is None or not scene.roots.row: return
        roots = scene.roots
        start = -1 if direction > 0 else 0
        current = start if scene.window is None else scene.window
        chosen = wrap(current + direction, len(roots.row))
        self.choose(chosen)
        x0, x1 = bounds(roots.ms(roots.row[chosen]))
        scene.view.focus(x0, x1, WINDOW_MARGIN)
        self.changed()
    async def step_through(self, direction: int) -> None:
        scene = self.scene
        if scene is None: return
        count = len(scene.roots.row) if scene.roots is not None else 0
        current = scene.window
        inside = current is not None and 0 <= current + direction < count
        if inside or (current is None and count > 0):
            self.step(direction)
            return
        await self.step_source(direction)
        scene = self.scene
        if scene is None or scene.roots is None or not scene.roots.row: return
        self.step(direction)
    def choose(self, window: int | None) -> None:
        scene = self.scene
        assert scene is not None
        if scene.window != window: scene.field = None
        if scene.window != window and window is not None: self.spawn(self.load_window(window))
        if window is None:
            scene.stft = None
            scene.reconstruction = None
            scene.residual = None
            scene.atomization = None
            scene.component = None
            scene.overlay = None
            self.report = None
            self.shape = None
        scene.window = window
    def focus_field(self, field: int) -> None:
        scene = self.scene
        assert scene is not None
        scene.field = field
        self.redraw()
    def clear(self) -> None:
        if self.scene is None: return
        self.scene.selection = None
        self.choose(None)
        self.changed()
    async def commit(self, window: int, column: str, value: float) -> None:
        if self.busy: return
        scene = self.scene
        assert scene is not None and scene.roots is not None
        roots = scene.roots
        row = roots.row[window]
        run = integer(row[roots.name.index("extraction")])
        index = integer(row[roots.name.index("index")])
        self.status("editing")
        await cli("edit_root", run, index, column, value)
        self.store(await read_roots(self.source_id))
        scene.window = window
        self.spawn(self.load_window(window))
        self.status("")
        self.changed()
    async def edit_factor(self, column: str, value: float) -> None:
        if self.busy: return
        self.status("editing")
        await cli("edit_dna", self.factor_id(), column, value)
        self.factors = await read("dna")
        self.status("")
        self.changed()
    async def copy_factor(self) -> None:
        if self.busy: return
        self.status("copying")
        await cli("copy_dna", self.factor_id())
        self.factors = await read("dna")
        self.factor()
        element("factor").selectedIndex = len(self.factors.row) - 1
        self.status("")
        self.changed()
    async def extract(self) -> None:
        if self.busy or self.scene is None: return
        self.status("extracting")
        await cli("root_atom", self.source_id, self.factor_id())
        self.sources = await read("source")
        self.store(await read_roots(self.source_id))
        self.choose(None)
        await self.survey()
        self.status("")
        self.changed()
    async def atomize(self) -> None:
        scene = self.scene
        if self.busy or scene is None or scene.window is None: return
        window = scene.window
        extraction, index, revision = self.root_key(window)
        self.status("atomizing")
        await stft(extraction, index, revision)
        await cli("atomize", extraction, index, self.ard_id(), 0)
        scene.atomization = await self.latest_atomization(window)
        scene.component = 0
        scene.reconstruction = None
        scene.residual = None
        await self.load_report()
        self.place_overlay()
        await self.survey()
        self.status("")
        await self.reload_layer()
        self.changed()
    async def analyze(self) -> None:
        scene = self.scene
        if self.busy or scene is None or scene.atomization is None: return
        atomization = scene.atomization
        self.status("analyzing")
        await cli("group_component", atomization.id)
        await self.load_report()
        self.place_overlay()
        self.status("")
        self.changed()
    async def reject_checked(self) -> None:
        if self.busy or self.scene is None: return
        self.status("deleting")
        for extraction, index in self.checked: await cli("edit_root", extraction, index, "rejected", 1)
        await self.reload_edited_roots()
    async def restore_checked(self) -> None:
        if self.busy or self.scene is None: return
        self.status("restoring")
        for extraction, index in self.checked: await cli("restore_root", extraction, index)
        await self.reload_edited_roots()
    async def reload_edited_roots(self) -> None:
        self.checked.clear()
        self.store(await read_roots(self.source_id))
        self.choose(None)
        await self.survey()
        self.status("")
        self.changed()
    async def rescan(self) -> None:
        if self.busy: return
        directory = js.window.prompt("source directory", self.source_directory)
        if directory is None: return
        self.source_directory = directory
        self.status("ingesting")
        await cli("ingest_source", directory)
        self.sources = await read("source")
        if self.source_index is not None:
            identifiers = [self.sources.key_cells(row)[0] for row in self.sources.row]
            self.source_index = identifiers.index(self.source_id)
        await self.survey()
        self.rescanned = js.Date.now()
        self.status("")
        self.changed()
    async def reset_extraction(self) -> None:
        if self.busy: return
        self.reset_ask()
        self.halt()
        self.status("resetting")
        await cli("clear_root_extraction")
        await self.survey()
        if self.scene is not None:
            self.store(await read_roots(self.source_id))
            self.choose(None)
        self.status("")
        self.changed()
    async def reload_layer(self) -> None:
        scene = self.scene
        if scene is None: return
        atomization = scene.atomization
        if atomization is None or self.spectrum.layer == "source":
            self.redraw()
            return
        if self.spectrum.layer == "reconstruction" and scene.reconstruction is None:
            self.status("reconstruction")
            raw = await cli("render_reconstruction", atomization.id)
            self.status("")
            if scene.atomization is atomization: scene.reconstruction = Stft.parse(raw)
        if self.spectrum.layer == "residual":
            stored = await read("residual", atomization.id)
            if not stored.row:
                self.status("residual")
                await cli("subtract_reconstruction", atomization.id)
                self.status("")
            raw = await cli("render_residual", atomization.id)
            if scene.atomization is atomization: scene.residual = Stft.parse(raw)
        self.redraw()
    def play_span(self) -> tuple[float, float]:
        scene = self.scene
        assert scene is not None
        if scene.window is not None:
            roots = scene.roots
            assert roots is not None
            row = roots.row[scene.window]
            return real(row[roots.name.index("adsr_x0")]), real(row[roots.name.index("adsr_x1")])
        if scene.selection: return bounds(scene.selection)
        return scene.view.src_x0, scene.view.src_x1
    def play_origin(self) -> float:
        scene = self.scene
        assert scene is not None
        if scene.window is None: return 0.0
        roots = scene.roots
        assert roots is not None
        return real(roots.row[scene.window][roots.name.index("adsr_x0")])
    def play_target(self) -> str:
        scene = self.scene
        assert scene is not None
        if scene.atomization is not None and scene.component is not None:
            return f"component {int(scene.atomization.index[scene.component])}"
        if scene.window is not None: return f"window {scene.window + 1}"
        return "source"
    async def play_toggle(self) -> None:
        if self.busy or self.scene is None: return
        scene = self.scene
        if scene.playhead is not None:
            await pause()
            self.pause_toggle()
            return
        x0, x1 = self.play_span()
        await self.sound(x0 - self.play_origin())
        self.play_clock = js.performance.now()
        self.paused_at = None
        scene.play_from = x0
        scene.play_until = x1
        scene.playhead = x0
        self.wave.on_frame = self.tick
        self.redraw()
    def pause_toggle(self) -> None:
        if self.paused_at is None:
            self.paused_at = js.performance.now()
            self.quiet()
            self.redraw()
            return
        self.play_clock += js.performance.now() - self.paused_at
        self.paused_at = None
        self.redraw()
    def tick(self) -> bool:
        scene = self.scene
        if scene is None or scene.playhead is None or scene.play_from is None or scene.play_until is None or self.paused_at is not None: return False
        if isinstance(self.gesture, Seek): return True
        position = scene.play_from + js.performance.now() - self.play_clock
        if position >= scene.play_until and self.looping:
            self.play_clock = js.performance.now()
            scene.playhead = scene.play_from
            self.spawn(self.sound(scene.play_from - self.play_origin()))
            self.transport()
            return True
        if position >= scene.play_until:
            self.halt()
            return False
        scene.playhead = position
        if not self.metering: self.spawn(self.meter())
        view = scene.view
        if self.following and not view.x0 <= position < view.x1:
            view.pan(position - view.x0)
            self.redraw()
        self.spectrum.request()
        self.ticks.request()
        self.transport()
        return True
    async def meter(self) -> None:
        self.metering = True
        shown = await level()
        self.metering = False
        scene = self.scene
        playing = scene is not None and scene.playhead is not None and self.paused_at is None
        if not playing or not shown: return
        floor = self.wave.look.floor
        left = height(shown[0], floor)
        right = height(shown[-1], floor)
        style = element("meter").style
        style.setProperty("--l", str(left))
        style.setProperty("--r", str(right))
    def quiet(self) -> None:
        style = element("meter").style
        style.setProperty("--l", "0")
        style.setProperty("--r", "0")
    def halt(self, event: object = None) -> None:
        if self.scene is None or self.scene.playhead is None: return
        self.spawn(stop())
        scene = self.scene
        scene.play_from = None
        scene.play_until = None
        scene.playhead = None
        self.paused_at = None
        self.wave.on_frame = None
        self.quiet()
        self.redraw()
    def wheel(self, event: Any) -> None:
        if self.scene is None: return
        event.preventDefault()
        scale = LINE_PX if event.deltaMode == 1 else 1
        horizontal = event.shiftKey or abs(event.deltaX) > abs(event.deltaY)
        view = self.scene.view
        if horizontal:
            delta = (event.deltaY if event.shiftKey else event.deltaX) * scale
            view.pan(delta*self.wave.ms_per_css_pixel())
        else:
            factor = math.exp(event.deltaY*scale*WHEEL)
            view.zoom(factor, self.wave.ms_at(event.clientX))
        self.redraw()
    def pointer_down(self, event: Any) -> None:
        if self.scene is None: return
        scene = self.scene
        self.wave.canvas.setPointerCapture(event.pointerId)
        if event.button == 1:
            self.gesture = Pan(event.clientX)
            return
        marker = None if self.busy else self.wave.marker_at(event.clientX)
        if marker:
            roots = scene.roots
            assert roots is not None
            self.choose(marker.window)
            scene.field = marker.field
            self.gesture = Drag(marker.window, roots.ms_name()[marker.field], False)
            self.changed()
            return
        self.gesture = Select(self.wave.ms_at(event.clientX), event.clientX, False)
    def pointer_move(self, event: Any) -> None:
        scene = self.scene
        if scene is None: return
        gesture = self.gesture
        if gesture is None:
            marker = None if self.busy else self.wave.marker_at(event.clientX)
            self.wave.canvas.style.cursor = "ew-resize" if marker else ""
            return
        if isinstance(gesture, Pan):
            scene.view.pan((gesture.x - event.clientX) * self.wave.ms_per_css_pixel())
            gesture.x = event.clientX
            self.redraw()
            return
        if isinstance(gesture, Drag):
            roots = scene.roots
            assert roots is not None
            ms = min(max(self.wave.ms_at(event.clientX), scene.view.src_x0), scene.view.src_x1)
            roots.row[gesture.window][roots.name.index(gesture.column)] = ms
            gesture.moved = True
            self.window()
            self.redraw()
            return
        if not isinstance(gesture, Select): return
        if abs(event.clientX - gesture.x) > DRAG_PX: gesture.moved = True
        if not gesture.moved: return
        scene.selection = [gesture.ms, self.wave.ms_at(event.clientX)]
        self.redraw()
    def pointer_up(self, event: Any) -> None:
        gesture = self.gesture
        if gesture is None: return
        self.gesture = None
        if isinstance(gesture, Pan): return
        if isinstance(gesture, Drag):
            scene = self.scene
            if scene is None or scene.roots is None or not gesture.moved: return
            roots = scene.roots
            value = roots.row[gesture.window][roots.name.index(gesture.column)]
            self.spawn(self.commit(gesture.window, gesture.column, real(value)))
            return
        if isinstance(gesture, Select) and not gesture.moved: self.pick(self.wave.ms_at(event.clientX))
        self.changed()
    def ticks_down(self, event: Any) -> None:
        if not self.ticks.on_progress(event.clientX): return
        self.ticks.canvas.setPointerCapture(event.pointerId)
        seek = Seek(self.seek_ms(event.clientX))
        self.gesture = seek
        assert self.scene is not None
        self.scene.playhead = seek.ms
        self.redraw()
    def ticks_move(self, event: Any) -> None:
        gesture = self.gesture
        self.ticks.canvas.style.cursor = "pointer" if self.ticks.on_progress(event.clientX) else ""
        if not isinstance(gesture, Seek) or self.scene is None: return
        gesture.ms = self.seek_ms(event.clientX)
        self.scene.playhead = gesture.ms
        self.redraw()
    def ticks_up(self, event: Any) -> None:
        gesture = self.gesture
        if not isinstance(gesture, Seek): return
        self.gesture = None
        self.seek(gesture.ms)
    def overview_down(self, event: Any) -> None:
        if self.scene is None or self.scene.overview is None: return
        self.overview.canvas.setPointerCapture(event.pointerId)
        self.panning = True
        self.overview_move(event)
    def overview_move(self, event: Any) -> None:
        scene = self.scene
        if not self.panning or scene is None: return
        view = scene.view
        view.pan(self.overview.ms_at(event.clientX) - (view.x0 + view.x1) / 2)
        self.redraw()
    def overview_up(self, event: Any) -> None:
        self.panning = False
    def seek_ms(self, client_x: float) -> float:
        scene = self.scene
        assert scene is not None and scene.play_from is not None and scene.play_until is not None
        ms = self.ticks.ms_at(client_x)
        return min(max(ms, scene.play_from), scene.play_until)
    def seek(self, ms: float) -> None:
        scene = self.scene
        assert scene is not None and scene.play_from is not None and scene.play_until is not None
        now = js.performance.now()
        paused = self.paused_at is not None
        self.play_clock = now - (ms - scene.play_from)
        if paused: self.paused_at = now
        scene.playhead = ms
        self.spawn(self.restart(ms - scene.play_from, paused))
        self.redraw()
    async def restart(self, offset: float, paused: bool) -> None:
        await self.sound(offset)
        if paused: await pause()
    def pick(self, ms: float) -> None:
        scene = self.scene
        assert scene is not None
        roots = scene.roots
        scene.selection = None
        found = None
        if roots is not None:
            for index, row in enumerate(roots.row):
                x0, x1 = bounds(roots.ms(row))
                if x0 <= ms < x1: found = index
        self.choose(found)

async def boot() -> Studio:
    sources = await read("source")
    factors = await read("dna")
    presets = await read("ard")
    studio = Studio(sources, factors, presets)
    await studio.survey()
    studio.listen()
    studio.changed()
    return studio
STUDIO = asyncio.ensure_future(boot())
