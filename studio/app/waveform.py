from __future__ import annotations
import asyncio
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Callable, overload
import js # pyright: ignore[reportMissingImports]
import numpy as np
from pyodide.ffi import create_proxy, to_js
from plot import rows
from relay import render_component_list, waveform
from library import Table
from view import View

TICK_PX = 96
LABEL_PX = 10
LABEL_BAND_PX = 22
DOT_PX = 6
HIT_PX = 4
GAIN = 1.0
COMPACT_ALPHA = 0.5
FLOOR = -60.0
PLANES = 2
PROGRESS_PX = 2
COMPONENT_COLUMNS = 2048
COMPONENT_KEYS = 5
QUANTIZATION = 24
HZ_TICK = (1, 2, 5)

@dataclass
class Marker:
    window: int
    field: int

@dataclass
class Look:
    compact: bool
    half:    bool
    decibel: bool
    floor:   float

@dataclass
class Atomization:
    id: int
    x0: float
    x1: float
    floor: float
    index: np.ndarray
    relevance: np.ndarray
    peak: np.ndarray
    onset: np.ndarray
    decay: np.ndarray
    activation: np.ndarray
    @staticmethod
    def parse(atomization_id: int, raw: bytes) -> Atomization:
        header = np.frombuffer(raw, "<f8", count=4)
        count = int(header[3])
        column = int(np.frombuffer(raw, "<u8", count=1, offset=32)[0])
        key = np.frombuffer(raw, "<f8", count=COMPONENT_KEYS*count, offset=40).reshape(count, COMPONENT_KEYS)
        activation = np.frombuffer(raw, "<f4", count=count*column, offset=40 + 8*COMPONENT_KEYS*count).reshape(count, column)
        floor = float(header[2])
        order = np.argsort(-key[:, 1])
        ordered = key[order]
        return Atomization(atomization_id, float(header[0]), float(header[1]), floor, ordered[:, 0].astype(int),
                           ordered[:, 1] / floor, ordered[:, 2], ordered[:, 3], ordered[:, 4], activation[order])
    @staticmethod
    async def fetch(atomization_id: int) -> Atomization:
        raw = await render_component_list(atomization_id, COMPONENT_COLUMNS)
        return Atomization.parse(atomization_id, raw)

@dataclass
class Band:
    offset:  float
    density: float
    db:      np.ndarray

@dataclass
class Stft:
    x0: float
    step:  float
    band:  list[Band]
    @property
    def columns(self) -> int:
        return self.band[0].db.shape[0]
    @property
    def x1(self) -> float:
        return self.x0 + (self.columns - 0.5)*self.step
    @property
    def z0(self) -> float:
        return self.band[0].offset
    @property
    def z1(self) -> float:
        top = self.band[-1]
        return top.offset + top.db.shape[1] / top.density
    def hz(self) -> np.ndarray:
        return np.concatenate([each.offset + np.arange(each.db.shape[1]) / each.density for each in self.band])
    @staticmethod
    def parse(raw: bytes) -> Stft:
        columns = int(np.frombuffer(raw, "<u8", count=1)[0])
        density, x0 = np.frombuffer(raw, "<f8", count=2, offset=8)
        count = int(np.frombuffer(raw, "<u8", count=1, offset=QUANTIZATION)[0])
        at = QUANTIZATION + 8
        axis = []
        for _ in range(count):
            rows = int(np.frombuffer(raw, "<u8", count=1, offset=at)[0])
            hz_density, offset = np.frombuffer(raw, "<f8", count=2, offset=at + 8)
            axis.append((rows, float(hz_density), float(offset)))
            at += QUANTIZATION
        band = []
        for rows, hz_density, offset in axis:
            db = np.frombuffer(raw, "<f4", count=columns*rows, offset=at).reshape(columns, rows)
            band.append(Band(offset, hz_density, db))
            at += 4*columns*rows
        return Stft(float(x0), 1 / float(density), band)

@dataclass
class Overlay:
    x0: float
    x1:   float
    low:   float
    high:  float
    label: str
    layer: Stft | None

@dataclass
class Scene:
    view:       View
    source_id:  int
    waveform:   Waveform
    shown:      tuple[float, float, int]
    roots:      Table       | None
    window:     int         | None
    field:      int         | None
    selection:  list[float] | None
    play_from:  float       | None
    play_until: float       | None
    playhead:   float       | None
    atomization: Atomization | None
    component:  int         | None
    stft:       Stft        | None
    reconstruction: Stft    | None
    residual:   Stft        | None
    overview:   Waveform    | None
    overlay:    Overlay     | None

@dataclass
class Palette:
    background: str
    lane: str
    grid: str
    inside: str
    outside: str
    chosen: str
    marker: str
    strong: str
    field: str
    selection: str
    playhead: str
    progress: str
    component: str
    overlay: str
    axis: str
    text: str
    label: str
    font: str

@dataclass
class Waveform:
    x0: float
    x1: float
    frame_ms: float
    floor: float
    frames_per_pixel: float
    ms: np.ndarray
    y0: np.ndarray
    y1: np.ndarray
    @property
    def channels(self) -> int:
        return len(self.y0)
    def peak(self) -> np.ndarray:
        return np.maximum(np.abs(self.y0), np.abs(self.y1)).max(axis=0)
    @staticmethod
    def parse(raw: bytes) -> Waveform:
        header = np.frombuffer(raw, "<f8", count=5)
        shape = np.frombuffer(raw, "<u8", count=2, offset=40)
        channels = int(shape[0])
        count = int(shape[1])
        ms = np.frombuffer(raw, "<f8", count=count, offset=56)
        values = np.frombuffer(raw, "<f4", count=channels*PLANES*count, offset=56 + 8*count)
        planes = values.reshape(channels, PLANES, count)
        return Waveform(float(header[0]), float(header[1]), float(header[2]), float(header[3]), float(header[4]), ms, planes[:, 0], planes[:, 1])

def millisecond(ms: float) -> str:
    return f"{ms:.1f}"
def token(style: Any, name: str) -> str:
    return style.getPropertyValue(name).strip()
def palette(element: Any) -> Palette:
    style = js.window.getComputedStyle(element)
    return Palette(
        background=token(style, "--wave-background"),
        lane=token(style, "--wave-lane"),
        grid=token(style, "--dark"),
        inside=token(style, "--wave-inside"),
        outside=token(style, "--wave-outside"),
        chosen=token(style, "--wave-chosen"),
        marker=token(style, "--wave-marker"),
        strong=token(style, "--wave-strong"),
        field=token(style, "--bright"),
        selection=token(style, "--wave-selection"),
        playhead=token(style, "--bright"),
        progress=token(style, "--faint"),
        component=token(style, "--wave-component"),
        overlay=token(style, "--wave-overlay"),
        axis=token(style, "--border"),
        text=token(style, "--dim"),
        label=token(style, "--light"),
        font=token(style, "--font-mono"))
@overload
def height(level: float, floor: float) -> float: ...
@overload
def height(level: np.ndarray, floor: float) -> np.ndarray: ...
def height(level: float | np.ndarray, floor: float) -> float | np.ndarray:
    return np.clip((level - floor) / -floor, 0, 1)
def amplitude(value: np.ndarray, decibel: bool, floor: float) -> np.ndarray:
    if not decibel: return value
    with np.errstate(divide="ignore", invalid="ignore"):
        level = 20 * np.log10(np.abs(value))
        scaled = height(level, floor)
    return np.sign(value) * scaled
def tick_step(minimum: float) -> float:
    base = 10**math.floor(math.log10(minimum))
    for multiple in (1, 2, 5):
        if base*multiple >= minimum: return base*multiple
    return base*10
def bounds(ms: list[float]) -> tuple[float, float]:
    return min(ms), max(ms)

class Surface(ABC):
    def __init__(self, canvas: Any) -> None:
        self.canvas = canvas
        self.context = canvas.getContext("2d")
        self.palette = palette(canvas)
        self.scene: Scene | None = None
        self.ratio = 1.0
        self.scheduled = False
        self.frame_proxy = create_proxy(self.frame)
        self.resize_proxy = create_proxy(self.resize)
        self.observer = js.ResizeObserver.new(self.resize_proxy)
        self.observer.observe(canvas)
    def resize(self, entries: object = None, observer: object = None) -> None:
        box = self.canvas.getBoundingClientRect()
        self.ratio = js.window.devicePixelRatio
        self.canvas.width = round(box.width*self.ratio)
        self.canvas.height = round(box.height*self.ratio)
        self.request()
    def request(self) -> None:
        if self.scheduled: return
        self.scheduled = True
        js.window.requestAnimationFrame(self.frame_proxy)
    def frame(self, timestamp: float = 0) -> None:
        self.scheduled = False
        self.draw()
    def active(self) -> Scene:
        scene = self.scene
        assert scene is not None
        return scene
    @overload
    def x(self, ms: float) -> float: ...
    @overload
    def x(self, ms: np.ndarray) -> np.ndarray: ...
    def x(self, ms: float | np.ndarray) -> float | np.ndarray:
        view = self.active().view
        return (ms - view.x0) / view.span * self.canvas.width
    def ms_at(self, client_x: float) -> float:
        box = self.canvas.getBoundingClientRect()
        return self.active().view.ms_at((client_x - box.left) / box.width)
    def ticks(self) -> tuple[int, int, float]:
        view = self.active().view
        step = tick_step(view.span / self.canvas.width*TICK_PX*self.ratio)
        return math.ceil(view.x0 / step), math.floor(view.x1 / step), step
    @abstractmethod
    def draw(self) -> None: ...

class Wave(Surface):
    def __init__(self, canvas: Any, layer: Any) -> None:
        super().__init__(canvas)
        self.layer = layer
        self.look = Look(compact=True, half=True, decibel=False, floor=FLOOR)
        self.on_frame: Callable[[], bool] | None = None
        self.task: asyncio.Task | None = None
        self.loaded: Callable[[], None] | None = None
    def frame(self, timestamp: float = 0) -> None:
        self.scheduled = False
        if self.on_frame and self.on_frame(): self.request()
        self.draw()
    def nearest(self, roots: Table, windows: list[int], pointer: float) -> Marker | None:
        nearest = None
        distance = HIT_PX*self.ratio
        for window in windows:
            for field, ms in enumerate(roots.ms(roots.row[window])):
                gap = abs(self.x(ms) - pointer)
                if gap > distance: continue
                nearest = Marker(window, field)
                distance = gap
        return nearest
    def marker_at(self, client_x: float) -> Marker | None:
        scene = self.scene
        if scene is None or scene.roots is None: return None
        roots = scene.roots
        box = self.canvas.getBoundingClientRect()
        pointer = (client_x - box.left) * self.ratio
        chosen = None if scene.window is None else self.nearest(roots, [scene.window], pointer)
        if chosen: return chosen
        return self.nearest(roots, list(range(len(roots.row))), pointer)
    def top(self) -> float:
        return LABEL_BAND_PX*self.ratio
    def ms_per_css_pixel(self) -> float:
        box = self.canvas.getBoundingClientRect()
        return self.active().view.span / box.width
    def frames_per_pixel(self) -> float:
        return self.active().waveform.frames_per_pixel
    async def fetch(self, scene: Scene, wanted: tuple[float, float, int]) -> None:
        x0, x1, width = wanted
        raw = await waveform(scene.source_id, x0, x1, width)
        scene.waveform = Waveform.parse(raw)
        scene.shown = wanted
        self.task = None
        if self.loaded is not None: self.loaded()
    def draw(self) -> None:
        context = self.context
        width = self.canvas.width
        context.clearRect(0, 0, width, self.canvas.height)
        if self.scene is None:
            self.layer.replaceChildren()
            return
        scene = self.scene
        lane_height = (self.canvas.height - self.top()) / self.lane_count()
        i0, i1, step = self.ticks()
        self.grid(i0, i1, step, lane_height)
        self.windows(scene)
        self.selection(scene)
        self.lanes(scene, width, lane_height)
        self.activation(scene)
        self.markers(scene)
        self.play(scene)
    def grid(self, i0: int, i1: int, step: float, lane_height: float) -> None:
        context = self.context
        width = self.canvas.width
        context.fillStyle = self.palette.grid
        for k in range(i0, i1 + 1):
            context.fillRect(round(self.x(k*step)), 0, self.ratio, self.canvas.height)
        context.fillStyle = self.palette.lane
        for lane in range(self.lane_count()):
            base, scale = self.origin(lane, lane_height)
            context.fillRect(0, round(base), width, self.ratio)
    def lane_count(self) -> int:
        return 1 if self.look.compact else self.active().waveform.channels
    def lane(self, channel: int) -> int:
        return 0 if self.look.compact else channel
    def origin(self, lane: int, lane_height: float) -> tuple[float, float]:
        if self.look.half: return self.top() + lane_height * (lane + 1), lane_height*GAIN
        return self.top() + lane_height * (lane + 0.5), lane_height*0.5*GAIN
    def level(self, value: np.ndarray) -> np.ndarray:
        if self.look.half: value = np.abs(value)
        return amplitude(value, self.look.decibel, self.look.floor)
    def alpha(self) -> float:
        overlaid = self.look.compact and self.active().waveform.channels > 1
        return COMPACT_ALPHA if overlaid else 1.0
    def chosen_span(self, scene: Scene) -> tuple[float, float] | None:
        if scene.window is None or scene.roots is None: return None
        return bounds(scene.roots.ms(scene.roots.row[scene.window]))
    def windows(self, scene: Scene) -> None:
        layer = self.layer
        roots = scene.roots
        if roots is None:
            layer.replaceChildren()
            return
        while layer.childElementCount < len(roots.row):
            layer.append(js.document.createElement("wave-window"))
        while layer.childElementCount > len(roots.row):
            last = layer.lastElementChild
            if last is None: break
            last.remove()
        for index, row in enumerate(roots.row):
            ms = roots.ms(row)
            left = self.x(min(ms))
            right = self.x(max(ms))
            band = layer.children[index]
            band.style.left = f"{left / self.ratio}px"
            band.style.width = f"{(right - left) / self.ratio}px"
            band.style.height = f"{self.canvas.height / self.ratio}px"
            band.style.opacity = "1" if index == scene.window else "0.35"
    def selection(self, scene: Scene) -> None:
        if not scene.selection: return
        left = self.x(min(scene.selection))
        right = self.x(max(scene.selection))
        self.context.fillStyle = self.palette.selection
        self.context.fillRect(left, 0, right - left, self.canvas.height)
    def lanes(self, scene: Scene, width: int, lane_height: float) -> None:
        view = scene.view
        wanted = (view.x0, view.x1, width)
        if wanted != scene.shown and self.task is None:
            self.task = asyncio.ensure_future(self.fetch(scene, wanted))
        if scene.waveform.frames_per_pixel >= 1:
            self.columns(scene, lane_height)
            return
        self.trace(scene, lane_height)
    def columns(self, scene: Scene, lane_height: float) -> None:
        context = self.context
        shown = scene.waveform
        x = self.x(shown.ms) + 0.5
        span = self.chosen_span(scene)
        inside = np.ones(len(shown.ms), bool) if span is None else (shown.ms >= span[0]) & (shown.ms <= span[1])
        context.lineWidth = self.ratio
        context.globalAlpha = self.alpha()
        for c in range(shown.channels):
            base, scale = self.origin(self.lane(c), lane_height)
            if self.look.half:
                peak = np.maximum(np.abs(shown.y0[c]), np.abs(shown.y1[c]))
                top = base - self.level(peak) * scale
                bottom = np.full_like(top, base)
            else:
                top = base - self.level(shown.y1[c]) * scale
                bottom = base - self.level(shown.y0[c]) * scale
            self.bars(x[~inside], top[~inside], bottom[~inside], self.palette.outside)
            self.bars(x[inside], top[inside], bottom[inside], self.palette.inside)
        context.globalAlpha = 1.0
    def activation(self, scene: Scene) -> None:
        atomization = scene.atomization
        if atomization is None or scene.component is None: return
        shape = atomization.activation[scene.component]
        peak = float(shape.max())
        ms = np.linspace(atomization.x0, atomization.x1, len(shape))
        x = self.x(ms)
        y = self.top() + (self.canvas.height - self.top())*(1 - shape / peak)
        line = []
        command = "M"
        for px, py in zip(x, y):
            line.append(f"{command}{px:.2f} {py:.2f}")
            command = "L"
        context = self.context
        context.strokeStyle = self.palette.component
        context.lineWidth = 1.5*self.ratio
        context.stroke(js.Path2D.new("".join(line)))
    def bars(self, x: np.ndarray, top: np.ndarray, bottom: np.ndarray, color: str) -> None:
        if len(x) == 0: return
        bottom = np.maximum(bottom, top + self.ratio)
        segment = []
        for px, y0, y1 in zip(x, top, bottom): segment.append(f"M{px:.2f} {y0:.2f}V{y1:.2f}")
        self.context.strokeStyle = color
        self.context.stroke(js.Path2D.new("".join(segment)))
    def trace(self, scene: Scene, lane_height: float) -> None:
        context = self.context
        shown = scene.waveform
        pixels_per_frame = 1 / shown.frames_per_pixel
        dot = 3*self.ratio
        x = self.x(shown.ms)
        context.strokeStyle = self.palette.inside
        context.fillStyle = self.palette.inside
        context.lineWidth = self.ratio
        context.globalAlpha = self.alpha()
        for c in range(shown.channels):
            base, scale = self.origin(self.lane(c), lane_height)
            y = base - self.level(shown.y0[c]) * scale
            line = []
            command = "M"
            for px, py in zip(x, y):
                line.append(f"{command}{px:.2f} {py:.2f}")
                command = "L"
            context.stroke(js.Path2D.new("".join(line)))
            if pixels_per_frame < DOT_PX*self.ratio: continue
            square = []
            for px, py in zip(x - dot / 2, y - dot / 2):
                square.append(f"M{px:.2f} {py:.2f}h{dot}v{dot}h{-dot}Z")
            context.fill(js.Path2D.new("".join(square)))
        context.globalAlpha = 1.0
    def markers(self, scene: Scene) -> None:
        roots = scene.roots
        if roots is None: return
        context = self.context
        name = roots.ms_name()
        line = (LABEL_PX + 2) * self.ratio
        context.font = f"{LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "top"
        for index, row in enumerate(roots.row):
            chosen = index == scene.window
            place = [round(self.x(ms)) for ms in roots.ms(row)]
            label = [name[field].removesuffix("_x") for field in range(len(place))]
            width = [context.measureText(each).width for each in label]
            stack = rows([(x + 4*self.ratio, x + 4*self.ratio + w) for x, w in zip(place, width)], 4*self.ratio)
            for field, x in enumerate(place):
                active = chosen and field == scene.field
                strong = chosen and label[field] == "peak"
                context.fillStyle = self.palette.field if active else self.palette.strong if strong else self.palette.marker
                context.fillRect(x, 0, self.ratio, self.canvas.height)
                if not chosen: continue
                context.fillStyle = self.palette.field if active else self.palette.label if strong else self.palette.text
                context.fillText(label[field], x + 4*self.ratio, 3*self.ratio + stack[field]*line)
    def play(self, scene: Scene) -> None:
        if scene.playhead is None: return
        context = self.context
        head = round(self.x(scene.playhead))
        context.shadowColor = self.palette.selection
        context.shadowBlur = 6*self.ratio
        context.fillStyle = self.palette.playhead
        context.fillRect(head, 0, self.ratio, self.canvas.height)
        context.shadowBlur = 0

class Ticks(Surface):
    def on_progress(self, client_x: float) -> bool:
        scene = self.scene
        if scene is None or scene.playhead is None or scene.play_from is None or scene.play_until is None: return False
        box = self.canvas.getBoundingClientRect()
        x = (client_x - box.left) * self.ratio
        return self.x(scene.play_from) <= x <= self.x(scene.play_until)
    def draw(self) -> None:
        context = self.context
        context.clearRect(0, 0, self.canvas.width, self.canvas.height)
        if self.scene is None: return
        scene = self.scene
        i0, i1, step = self.ticks()
        decimals = max(0, -math.floor(math.log10(step)))
        context.font = f"{LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "middle"
        middle = self.canvas.height / 2
        for k in range(i0, i1 + 1):
            ms = k*step
            x = round(self.x(ms))
            context.fillStyle = self.palette.axis
            context.fillRect(x, 0, self.ratio, 4*self.ratio)
            context.fillStyle = self.palette.text
            context.fillText(f"{ms:.{decimals}f}", x + 3*self.ratio, middle)
        if scene.play_from is None or scene.play_until is None: return
        bar = PROGRESS_PX*self.ratio
        bottom = self.canvas.height - bar
        position = scene.play_from if scene.playhead is None else scene.playhead
        left = self.x(scene.play_from)
        context.fillStyle = self.palette.progress
        context.fillRect(left, bottom, self.x(scene.play_until) - left, bar)
        context.fillStyle = self.palette.playhead
        context.fillRect(left, bottom, self.x(position) - left, bar)

class Overview(Surface):
    def place(self, ms: float | np.ndarray) -> float | np.ndarray:
        shown = self.active().overview
        assert shown is not None
        return (ms - shown.x0) / (shown.x1 - shown.x0) * self.canvas.width
    def ms_at(self, client_x: float) -> float:
        shown = self.active().overview
        assert shown is not None
        box = self.canvas.getBoundingClientRect()
        return shown.x0 + (client_x - box.left) / box.width * (shown.x1 - shown.x0)
    def draw(self) -> None:
        context = self.context
        width = self.canvas.width
        tall = self.canvas.height
        context.clearRect(0, 0, width, tall)
        scene = self.scene
        if scene is None or scene.overview is None: return
        shown = scene.overview
        peak = shown.peak()
        level = peak / max(float(peak.max()), 1e-12)
        x = self.place(shown.ms)
        view = scene.view
        inside = (shown.ms >= view.x0) & (shown.ms <= view.x1)
        pad = 3*self.ratio
        reach = tall - 2*pad
        middle = tall / 2
        for color, mask in ((self.palette.axis, ~inside), (self.palette.text, inside)):
            segment = []
            for px, share in zip(x[mask], level[mask]):
                half = max(self.ratio / 2, share*reach / 2)
                segment.append(f"M{px:.2f} {middle - half:.2f}V{middle + half:.2f}")
            context.strokeStyle = color
            context.lineWidth = self.ratio
            context.stroke(js.Path2D.new("".join(segment)))
        left = self.place(view.x0)
        right = self.place(view.x1)
        context.fillStyle = self.palette.chosen
        context.fillRect(left, self.ratio, right - left, tall - 2*self.ratio)
        context.strokeStyle = self.palette.text
        context.strokeRect(left + self.ratio / 2, self.ratio / 2, right - left - self.ratio, tall - self.ratio)
        if scene.window is None or scene.roots is None: return
        x0, x1 = bounds(scene.roots.ms(scene.roots.row[scene.window]))
        context.fillStyle = self.palette.field
        context.fillRect(self.place(x0), tall - 2*self.ratio, self.place(x1) - self.place(x0), 2*self.ratio)

def hz_text(hz: float) -> str:
    if hz >= 1: return f"{hz:g}k"
    return f"{hz*1000:g}"
def hz_ticks(z0: float, z1: float) -> list[float]:
    found = []
    for power in range(math.floor(math.log10(z0)), math.ceil(math.log10(z1)) + 1):
        for multiple in HZ_TICK:
            hz = multiple*10**power
            if z0 <= hz <= z1: found.append(hz)
    return found
def rgb(color: str) -> tuple[int, int, int]:
    return int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)

class Spectrum(Surface):
    def __init__(self, canvas: Any) -> None:
        super().__init__(canvas)
        self.floor = FLOOR
        self.layer = "source"
        self.image = js.document.createElement("canvas")
        self.component_image = js.document.createElement("canvas")
        self.painted: tuple[int, int, int, float] | None = None
        self.component_painted: tuple[int, int, int, float] | None = None
    def y(self, hz: np.ndarray | float, stft: Stft) -> np.ndarray | float:
        position = (np.log(hz) - math.log(stft.z0)) / (math.log(stft.z1) - math.log(stft.z0))
        return (1 - position) * self.canvas.height
    def shown(self, scene: Scene) -> Stft | None:
        if self.layer == "reconstruction" and scene.reconstruction is not None: return scene.reconstruction
        if self.layer == "residual" and scene.residual is not None: return scene.residual
        return scene.stft
    def paint(self, stft: Stft, reference: Stft, image: Any, color: str) -> None:
        height = self.canvas.height
        hz = np.concatenate([each.offset + (np.arange(each.db.shape[1]) + 0.5) / each.density for each in stft.band])
        db = np.concatenate([each.db for each in stft.band], axis=1).T
        pixel = np.clip(self.y(hz, stft).astype(int), 0, height - 1)
        centre = np.exp(math.log(stft.z1) - (np.arange(height) + 0.5) / height * (math.log(stft.z1) - math.log(stft.z0)))
        above = np.clip(np.searchsorted(hz, centre), 1, len(hz) - 1)
        closer = np.where(centre - hz[above - 1] < hz[above] - centre, above - 1, above)
        shown = db[closer]
        np.maximum.at(shown, pixel, db)
        peak = max(float(each.db.max()) for each in reference.band)
        level = np.clip((shown - (peak + self.floor)) / -self.floor, 0, 1)
        red, green, blue = rgb(color)
        rgba = np.empty((height, stft.columns, 4), np.uint8)
        rgba[..., 0] = red
        rgba[..., 1] = green
        rgba[..., 2] = blue
        rgba[..., 3] = (level*255).astype(np.uint8)
        pixels = js.Uint8ClampedArray.new(to_js(rgba.tobytes()))
        image.width = stft.columns
        image.height = height
        image.getContext("2d").putImageData(js.ImageData.new(pixels, stft.columns, height), 0, 0)
    def draw(self) -> None:
        context = self.context
        width = self.canvas.width
        height = self.canvas.height
        context.clearRect(0, 0, width, height)
        scene = self.scene
        if scene is None or scene.stft is None: return
        stft = self.shown(scene)
        assert stft is not None
        wanted = (id(stft), id(scene.stft), height, self.floor)
        if self.painted != wanted:
            self.paint(stft, scene.stft, self.image, self.palette.inside)
            self.painted = wanted
        left = self.x(stft.x0 - stft.step / 2)
        right = self.x(stft.x1)
        context.imageSmoothingEnabled = False
        context.drawImage(self.image, left, 0, right - left, height)
        context.fillStyle = self.palette.strong
        context.fillRect(round(left), 0, self.ratio, height)
        context.fillRect(round(right), 0, self.ratio, height)
        self.overlay(scene, stft)
        if scene.playhead is not None:
            context.fillStyle = self.palette.strong
            context.fillRect(round(self.x(scene.playhead)), 0, self.ratio, height)
        self.hz_axis(stft)
    def overlay(self, scene: Scene, stft: Stft) -> None:
        box = scene.overlay
        if box is None: return
        context = self.context
        left = self.x(box.x0)
        right = self.x(box.x1)
        y0 = self.y(min(max(box.high, stft.z0), stft.z1), stft)
        y1 = self.y(min(max(box.low, stft.z0), stft.z1), stft)
        if box.layer is None:
            context.fillStyle = self.palette.overlay
            context.fillRect(left, y0, right - left, y1 - y0)
        else:
            self.component_layer(scene, box.layer)
        context.strokeStyle = self.palette.component
        context.lineWidth = self.ratio
        context.strokeRect(left, y0, right - left, y1 - y0)
        context.font = f"bold {LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "top"
        context.fillStyle = self.palette.component
        context.fillText(box.label, left + 4*self.ratio, y0 + self.ratio)
    def component_layer(self, scene: Scene, layer: Stft) -> None:
        assert scene.stft is not None
        height = self.canvas.height
        wanted = (id(layer), id(scene.stft), height, self.floor)
        if self.component_painted != wanted:
            self.paint(layer, scene.stft, self.component_image, self.palette.component)
            self.component_painted = wanted
        left = self.x(layer.x0 - layer.step / 2)
        right = self.x(layer.x1)
        self.context.drawImage(self.component_image, left, 0, right - left, height)
    def hz_axis(self, stft: Stft) -> None:
        context = self.context
        context.font = f"{LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "bottom"
        for hz in hz_ticks(stft.z0, stft.z1):
            y = round(self.y(hz, stft))
            context.fillStyle = self.palette.axis
            context.fillRect(0, y, 5*self.ratio, self.ratio)
            context.fillStyle = self.palette.text
            label = min(max(y + 0.5*LABEL_PX*self.ratio, (LABEL_PX + 1)*self.ratio), self.canvas.height - self.ratio)
            context.fillText(hz_text(hz), 7*self.ratio, label)
