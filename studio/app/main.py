from __future__ import annotations

import asyncio
import math
from collections.abc import Coroutine
from dataclasses import dataclass
from typing import Any

import js # pyright: ignore[reportMissingImports]
from pyodide.ffi import create_proxy, to_js

import mixer_fader
from library import Table, Value, integer, listing, present, real
from mixer_fader import FULL, channels
from relay import cli, level, mixer, pause, play, read, stop, waveform
from view import View, cells, definitions, fields
from waveform import Scene, Waveform, Wave, height, millisecond


PAN = 0.25
ZOOM = 2
WHEEL = 0.002
LINE_PX = 16
DRAG_PX = 3
WINDOW_MARGIN = 0.25
MINIMUM_FRAMES = 16

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
def bounds(ms: list[float]) -> tuple[float, float]:
    return min(ms), max(ms)

Gesture = Select | Pan | Drag | Seek
class Studio:
    def __init__(self, sources: Table, factors: Table) -> None:
        self.sources = sources
        self.factors = factors
        self.wave = Wave(element("wave"), element("window-layer"))
        self.wave.loaded = self.redraw
        self.second = True
        self.source_index: int | None = None
        self.source_id = 0
        self.roots: Table | None = None
        self.windows: dict[int, Table] = {}
        self.deleted = False
        self.scene: Scene | None = None
        self.busy = ""
        self.gesture: Gesture | None = None
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
            "ArrowLeft":  (self.pan,         (-PAN,)),
            "ArrowRight": (self.pan,         (PAN,)),
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
        self.on(element("millisecond"), "click", lambda event: self.unit(False))
        self.on(element("second"), "click", lambda event: self.unit(True))
        self.on(element("factor"), "change", self.factor_change)
        self.on(element("center"), "click", self.center)
        self.on(element("zoom-out"), "click", lambda event: self.zoom(ZOOM))
        self.on(element("zoom-window"), "click", self.zoom_window)
        self.on(element("zoom-in"), "click", lambda event: self.zoom(1 / ZOOM))
        self.on(element("loop"), "click", self.loop)
        self.on(element("follow"), "click", self.follow)
        self.on(element("deleted"), "click", self.show_deleted)
        self.on(element("compact"), "click", lambda event: self.look("compact"))
        self.on(element("half"), "click", lambda event: self.look("half"))
        self.on(element("logarithmic"), "click", self.decibel)
        self.on(element("floor"), "change", self.floor)
        self.on(element("play"), "click", lambda event: self.spawn(self.play_toggle()))
        self.on(element("stop"), "click", self.halt)
        self.on(element("previous"), "click", lambda event: self.step(-1))
        self.on(element("next"), "click", lambda event: self.step(1))
        self.on(element("reject"), "click", self.reject)
        self.on(element("axis-x"), "focusin", self.field_focus)
        self.on(element("axis-x"), "change", self.field_change)
        self.on(element("extract"), "click", lambda event: self.spawn(self.extract_all()))
        self.on(element("extract-selection"), "click", lambda event: self.spawn(self.extract_selection()))
        self.on(element("master"), "input", self.master_change)
        self.on(element("pan"), "input", self.balance_change)
        self.on(element("mono"), "click", self.mono_toggle)
        self.on(canvas, "wheel", self.wheel, {"passive": False})
        self.on(canvas, "pointerdown", self.pointer_down)
        self.on(canvas, "pointermove", self.pointer_move)
        self.on(canvas, "pointerup", self.pointer_up)
        mixer_fader.install(self.on)
    def key(self, event: Any) -> None:
        if event.target.closest("input, select"): return
        action = self.keys.get(event.key)
        if action is None: return
        event.preventDefault()
        method, argument = action
        result = method(*argument)
        if isinstance(result, Coroutine): self.spawn(result)
    def press(self, event: Any) -> None:
        if event.target.closest("button"): event.preventDefault()
    def chosen(self, event: Any) -> None:
        button = event.target.closest("button")
        if button is None: return
        index = int(button.dataset.source)
        if button.classList.contains("window"):
            self.spawn(self.open_window(index, int(button.dataset.window)))
            return
        self.spawn(self.open_rec(index))
    def factor_change(self, event: object) -> None:
        self.factor()
    def decibel(self, event: object = None) -> None:
        self.look("decibel")
    def floor(self, event: Any) -> None:
        self.wave.look.floor = float(event.target.value)
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
    async def open(self, index: int) -> None:
        if self.busy: return
        self.halt()
        source_id = self.sources.key_cells(self.sources.row[index])[0]
        self.status(f"{source_id} decoding")
        width = self.wave.canvas.width
        raw = await waveform(source_id, -math.inf, math.inf, width)
        shown = Waveform.parse(raw)
        roots = await read("root", source_id)
        self.roots = roots
        self.source_index = index
        self.source_id = source_id
        self.channel_count = integer(self.source_cell("channels"))
        self.send_mixer()
        self.windows[source_id] = roots
        view = View(shown.start, shown.end, MINIMUM_FRAMES*shown.frame_ms)
        self.scene = Scene(
            view=view, source_id=source_id, waveform=shown,
            shown=(view.t0, view.t1, width),
            roots=self.kept(roots), window=None, field=None,
            selection=None, play_from=None, play_until=None,
            playhead=None)
        self.wave.scene = self.scene
        self.status("")
        self.changed()
    async def open_rec(self, index: int) -> None:
        source_id = self.sources.key_cells(self.sources.row[index])[0]
        if index != self.source_index:
            await self.open(index)
            return
        if source_id in self.windows: del self.windows[source_id]
        else:
            roots = self.roots
            assert roots is not None
            self.windows[source_id] = roots
        self.changed()
    async def open_window(self, index: int, window: int) -> None:
        if index != self.source_index: await self.open(index)
        scene = self.scene
        assert scene is not None and scene.roots is not None
        roots = scene.roots
        self.choose(window)
        lo, hi = bounds(roots.ms(roots.row[window]))
        scene.view.focus(lo, hi, WINDOW_MARGIN)
        self.changed()
    def changed(self) -> None:
        self.library()
        self.source()
        self.factor()
        self.window()
        self.redraw()
    def kept(self, roots: Table) -> Table:
        if self.deleted: return roots
        rejected = roots.name.index("rejected")
        row = []
        for each in roots.row:
            if not each[rejected]: row.append(each)
        return Table(roots.column, row, roots.setting)
    def store(self, roots: Table) -> None:
        self.roots = roots
        if self.source_id in self.windows: self.windows[self.source_id] = roots
        scene = self.scene
        assert scene is not None
        scene.roots = self.kept(roots)
    def library(self) -> None:
        chosen = self.scene.window if self.scene else None
        shown = {}
        for source_id, roots in self.windows.items(): shown[source_id] = self.kept(roots)
        listing(element("list"), self.sources, shown, self.second, self.source_index, chosen)
        element("deleted").classList.toggle("active", self.deleted)
        element("millisecond").classList.toggle("active", not self.second)
        element("second").classList.toggle("active", self.second)
    def source(self) -> None:
        if self.source_index is None: return
        pair = cells(self.sources, self.source_index, True, self.sources.units())
        definitions(element("source-facts"), pair)
    def factor(self) -> None:
        select = element("factor")
        chosen = max(select.selectedIndex, 0)
        select.replaceChildren()
        name = self.factors.name.index("name")
        column = self.factors.column[name]
        setting = self.factors.setting[name]
        for row in self.factors.row:
            option = js.document.createElement("option")
            option.textContent = present(row[name], column, setting, False, False)
            select.append(option)
        select.selectedIndex = chosen
        definitions(element("factor-fields"), cells(self.factors, chosen, False, [""]))
    def window(self) -> None:
        axis_x = element("axis-x")
        axis_y = element("axis-y")
        key = element("window-key")
        scene = self.scene
        roots = scene.roots if scene else None
        reject = element("reject")
        if scene is None or roots is None or scene.window is None:
            reject.disabled = True
            reject.classList.remove("active")
            element("window-index").textContent = ""
            axis_x.replaceChildren()
            definitions(axis_y, [])
            key.textContent = ""
            return
        window = scene.window
        row = roots.row[window]
        rejected = row[roots.name.index("rejected")]
        reject.disabled = bool(self.busy)
        reject.classList.toggle("active", bool(rejected))
        element("window-index").textContent = f"{window + 1} / {len(roots.row)}"
        fields(axis_x, roots, row, scene.field)
        definitions(axis_y, cells(roots, window, True, ["db"]))
        shown = []
        for name, text in cells(roots, window, True, [""]): shown.append(f"{name} {text}")
        key.textContent = "  ".join(shown)
    def redraw(self) -> None:
        self.wave.request()
        scene = self.scene
        ready = scene is not None and not self.busy
        element("extract").disabled = not ready
        element("extract-selection").disabled = scene is None or bool(self.busy) or not scene.selection
        element("loop").classList.toggle("active", self.looping)
        element("follow").classList.toggle("active", self.following)
        element("mono").classList.toggle("active", self.mono)
        if scene is None: return
        view = scene.view
        frames_per_pixel = self.wave.frames_per_pixel()
        decimals = 2 if frames_per_pixel < 10 else 0
        element("range").textContent = f"{millisecond(view.t0)} – {millisecond(view.t1)} ms"
        element("density").textContent = f"{frames_per_pixel:.{decimals}f} frames / pixels"
        look = self.wave.look
        element("compact").classList.toggle("active", look.compact)
        element("half").classList.toggle("active", look.half)
        element("logarithmic").classList.toggle("active", look.decibel)
        if scene.selection:
            lo, hi = bounds(scene.selection)
            element("selection").textContent = f"selection {millisecond(lo)} – {millisecond(hi)} ms"
        else:
            element("selection").textContent = ""
        self.transport()
    def transport(self) -> None:
        scene = self.scene
        assert scene is not None
        playing = scene.playhead is not None and self.paused_at is None
        if scene.playhead is None: scene.play_from, scene.play_until = self.play_span()
        element("play").textContent = "PAUSE" if playing else "PLAY"
        element("target").textContent = self.play_target() if scene.playhead is None else ""
    def center(self, event: object = None) -> None:
        if self.scene is None: return
        self.scene.view.whole()
        self.redraw()
    def zoom_window(self, event: object = None) -> None:
        scene = self.scene
        if scene is None or scene.window is None or scene.roots is None: return
        roots = scene.roots
        lo, hi = bounds(roots.ms(roots.row[scene.window]))
        scene.view.focus(lo, hi, WINDOW_MARGIN)
        self.redraw()
    def unit(self, second: bool) -> None:
        self.second = second
        self.library()
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
        self.changed()
    def pan(self, fraction: float) -> None:
        if self.scene is None: return
        self.scene.view.pan(self.scene.view.span*fraction)
        self.redraw()
    def zoom(self, factor: float) -> None:
        if self.scene is None: return
        view = self.scene.view
        view.zoom(factor, (view.t0 + view.t1) / 2)
        self.redraw()
    def step(self, direction: int) -> None:
        scene = self.scene
        if scene is None or scene.roots is None or not scene.roots.row: return
        roots = scene.roots
        start = -1 if direction > 0 else 0
        current = start if scene.window is None else scene.window
        chosen = wrap(current + direction, len(roots.row))
        self.choose(chosen)
        lo, hi = bounds(roots.ms(roots.row[chosen]))
        scene.view.focus(lo, hi, WINDOW_MARGIN)
        self.changed()
    def choose(self, window: int | None) -> None:
        scene = self.scene
        assert scene is not None
        if scene.window != window: scene.field = None
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
        run = integer(row[roots.name.index("run")])
        index = integer(row[roots.name.index("index")])
        self.status("editing")
        await cli("edit", run, index, column, value)
        self.store(await read("root", self.source_id))
        scene.window = window
        self.status("")
        self.changed()
    async def extract(self, start: float, end: float) -> None:
        scene = self.scene
        assert scene is not None
        self.status("extracting")
        await cli("run_atomextract", self.source_id, self.factor_id(), start, end)
        self.sources = await read("source_info")
        self.store(await read("root", self.source_id))
        self.choose(None)
        self.status("")
        self.changed()
    async def extract_all(self) -> None:
        if self.busy or self.scene is None: return
        await self.extract(self.scene.view.start, self.scene.view.end)
    async def extract_selection(self) -> None:
        if self.busy or self.scene is None or not self.scene.selection: return
        lo, hi = bounds(self.scene.selection)
        await self.extract(lo, hi)
    def play_span(self) -> tuple[float, float]:
        scene = self.scene
        assert scene is not None
        if scene.window is not None:
            roots = scene.roots
            assert roots is not None
            return bounds(roots.ms(roots.row[scene.window]))
        if scene.selection: return bounds(scene.selection)
        return scene.view.t0, scene.view.t1
    def play_target(self) -> str:
        scene = self.scene
        assert scene is not None
        if scene.window is not None: return f"window {scene.window + 1}"
        if scene.selection: return "selection"
        return "view"
    async def play_toggle(self) -> None:
        if self.busy or self.scene is None: return
        scene = self.scene
        if scene.playhead is not None:
            await pause()
            self.pause_toggle()
            return
        start, end = self.play_span()
        await play(self.source_id, start, end)
        self.play_clock = js.performance.now()
        self.paused_at = None
        scene.play_from = start
        scene.play_until = end
        scene.playhead = start
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
            self.spawn(play(self.source_id, scene.play_from, scene.play_until))
            self.transport()
            return True
        if position >= scene.play_until:
            self.halt()
            return False
        scene.playhead = position
        if not self.metering: self.spawn(self.meter())
        view = scene.view
        if self.following and not view.t0 <= position < view.t1:
            view.pan(position - view.t0)
            self.redraw()
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
        if self.wave.on_progress(event.clientX, event.clientY):
            seek = Seek(self.seek_ms(event.clientX))
            self.gesture = seek
            scene.playhead = seek.ms
            self.redraw()
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
            seekable = self.wave.on_progress(event.clientX, event.clientY)
            cursor = ""
            if marker: cursor = "ew-resize"
            if seekable: cursor = "pointer"
            self.wave.canvas.style.cursor = cursor
            return
        if isinstance(gesture, Seek):
            gesture.ms = self.seek_ms(event.clientX)
            scene.playhead = gesture.ms
            self.redraw()
            return
        if isinstance(gesture, Pan):
            scene.view.pan((gesture.x - event.clientX) * self.wave.ms_per_css_pixel())
            gesture.x = event.clientX
            self.redraw()
            return
        if isinstance(gesture, Drag):
            roots = scene.roots
            assert roots is not None
            ms = min(max(self.wave.ms_at(event.clientX), scene.view.start), scene.view.end)
            roots.row[gesture.window][roots.name.index(gesture.column)] = ms
            gesture.moved = True
            self.window()
            self.redraw()
            return
        if abs(event.clientX - gesture.x) > DRAG_PX: gesture.moved = True
        if not gesture.moved: return
        scene.selection = [gesture.ms, self.wave.ms_at(event.clientX)]
        self.redraw()
    def pointer_up(self, event: Any) -> None:
        gesture = self.gesture
        if gesture is None: return
        self.gesture = None
        if isinstance(gesture, Pan): return
        if isinstance(gesture, Seek):
            self.seek(gesture.ms)
            return
        if isinstance(gesture, Drag):
            scene = self.scene
            if scene is None or scene.roots is None or not gesture.moved: return
            roots = scene.roots
            value = roots.row[gesture.window][roots.name.index(gesture.column)]
            self.spawn(self.commit(gesture.window, gesture.column, real(value)))
            return
        if not gesture.moved: self.pick(self.wave.ms_at(event.clientX))
        self.changed()
    def seek_ms(self, client_x: float) -> float:
        scene = self.scene
        assert scene is not None and scene.play_from is not None and scene.play_until is not None
        ms = self.wave.ms_at(client_x)
        return min(max(ms, scene.play_from), scene.play_until)
    def seek(self, ms: float) -> None:
        scene = self.scene
        assert scene is not None and scene.play_from is not None and scene.play_until is not None
        now = js.performance.now()
        paused = self.paused_at is not None
        self.play_clock = now - (ms - scene.play_from)
        if paused: self.paused_at = now
        scene.playhead = ms
        self.spawn(self.restart(ms, scene.play_until, paused))
        self.redraw()
    async def restart(self, ms: float, until: float, paused: bool) -> None:
        await play(self.source_id, ms, until)
        if paused: await pause()
    def pick(self, ms: float) -> None:
        scene = self.scene
        assert scene is not None
        roots = scene.roots
        scene.selection = None
        found = None
        if roots is not None:
            for index, row in enumerate(roots.row):
                lo, hi = bounds(roots.ms(row))
                if lo <= ms < hi: found = index
        self.choose(found)

async def boot() -> Studio:
    sources = await read("source_info")
    factors = await read("dna")
    studio = Studio(sources, factors)
    studio.listen()
    studio.changed()
    return studio
STUDIO = asyncio.ensure_future(boot())
