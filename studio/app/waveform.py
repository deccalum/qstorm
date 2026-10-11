from __future__ import annotations
import asyncio
import math
from dataclasses import dataclass
from typing import Any, Callable, overload
import js # pyright: ignore[reportMissingImports]
import numpy as np
from pyodide.ffi import create_proxy
from relay import waveform
from library import Table
from view import View

AXIS_PX = 22
TICK_PX = 96
LABEL_PX = 11
DOT_PX = 6
PROGRESS_PX = 16
CLOCK_PAD_PX = 4
HIT_PX = 4
GAIN = 1.0
COMPACT_ALPHA = 0.5
FLOOR = -60.0
PLANES = 2

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

@dataclass
class Palette:
    background: str
    lane: str
    grid: str
    wave: str
    chosen: str
    marker: str
    field: str
    selection: str
    playhead: str
    progress: str
    axis: str
    text: str
    font: str

@dataclass
class Waveform:
    start: float
    end: float
    frame_ms: float
    floor: float
    frames_per_pixel: float
    ms: np.ndarray
    lo: np.ndarray
    hi: np.ndarray
    @property
    def channels(self) -> int:
        return len(self.lo)
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
        wave=token(style, "--light"),
        chosen=token(style, "--wave-chosen"),
        marker=token(style, "--wave-marker"),
        field=token(style, "--bright"),
        selection=token(style, "--wave-selection"),
        playhead=token(style, "--bright"),
        progress=token(style, "--dim"),
        axis=token(style, "--border"),
        text=token(style, "--mid"),
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

class Wave:
    def __init__(self, canvas: Any, layer: Any) -> None:
        self.canvas = canvas
        self.layer = layer
        self.context = canvas.getContext("2d")
        self.palette = palette(canvas)
        self.scene: Scene | None = None
        self.look = Look(compact=True, half=True, decibel=False, floor=FLOOR)
        self.ratio = 1.0
        self.scheduled = False
        self.on_frame: Callable[[], bool] | None = None
        self.task: asyncio.Task | None = None
        self.loaded: Callable[[], None] | None = None
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
        if self.on_frame and self.on_frame(): self.request()
        self.draw()
    def active(self) -> Scene:
        scene = self.scene
        assert scene is not None
        return scene
    def ms_at(self, client_x: float) -> float:
        box = self.canvas.getBoundingClientRect()
        return self.active().view.ms_at((client_x - box.left) / box.width)
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
    def on_progress(self, client_x: float, client_y: float) -> bool:
        scene = self.scene
        if scene is None or scene.playhead is None or scene.play_from is None or scene.play_until is None: return False
        box = self.canvas.getBoundingClientRect()
        x = (client_x - box.left) * self.ratio
        y = (client_y - box.top) * self.ratio
        top = self.canvas.height - (PROGRESS_PX + HIT_PX) * self.ratio
        inside_y = top <= y
        inside_x = self.x(scene.play_from) <= x <= self.x(scene.play_until)
        return inside_y and inside_x
    def lane_area(self) -> float:
        return self.canvas.height - (AXIS_PX + PROGRESS_PX) * self.ratio
    def ms_per_css_pixel(self) -> float:
        box = self.canvas.getBoundingClientRect()
        return self.active().view.span / box.width
    def frames_per_pixel(self) -> float:
        return self.active().waveform.frames_per_pixel
    async def fetch(self, scene: Scene, wanted: tuple[float, float, int]) -> None:
        t0, t1, width = wanted
        raw = await waveform(scene.source_id, t0, t1, width)
        scene.waveform = Waveform.parse(raw)
        scene.shown = wanted
        self.task = None
        if self.loaded is not None: self.loaded()
    @overload
    def x(self, ms: float) -> float: ...
    @overload
    def x(self, ms: np.ndarray) -> np.ndarray: ...
    def x(self, ms: float | np.ndarray) -> float | np.ndarray:
        view = self.active().view
        return (ms - view.t0) / view.span * self.canvas.width
    def draw(self) -> None:
        context = self.context
        width = self.canvas.width
        height = self.canvas.height
        context.clearRect(0, 0, width, height)
        if self.scene is None:
            self.layer.replaceChildren()
            return
        scene = self.scene
        lane_area = self.lane_area()
        lane_height = lane_area / self.lane_count()
        step = tick_step(scene.view.span / width*TICK_PX*self.ratio)
        first = math.ceil(scene.view.t0 / step)
        last = math.floor(scene.view.t1 / step)
        self.grid(scene, first, last, step, lane_area, lane_height)
        self.windows(scene, lane_area)
        self.selection(scene, lane_area)
        self.lanes(scene, width, lane_height)
        self.markers(scene, lane_area)
        self.play(scene, lane_area)
        self.axis(first, last, step, lane_area)
    def grid(self, scene: Scene, first: int, last: int, step: float, lane_area: float, lane_height: float) -> None:
        context = self.context
        width = self.canvas.width
        context.fillStyle = self.palette.grid
        for k in range(first, last + 1):
            context.fillRect(round(self.x(k*step)), 0, self.ratio, lane_area)
        context.fillStyle = self.palette.lane
        for lane in range(self.lane_count()):
            base, scale = self.origin(lane, lane_height)
            context.fillRect(0, round(base), width, self.ratio)
            if lane > 0: context.fillRect(0, round(lane_height*lane), width, self.ratio)
    def lane_count(self) -> int:
        return 1 if self.look.compact else self.active().waveform.channels
    def lane(self, channel: int) -> int:
        return 0 if self.look.compact else channel
    def origin(self, lane: int, lane_height: float) -> tuple[float, float]:
        if self.look.half: return lane_height * (lane + 1), lane_height*GAIN
        return lane_height * (lane + 0.5), lane_height*0.5*GAIN
    def level(self, value: np.ndarray) -> np.ndarray:
        if self.look.half: value = np.abs(value)
        return amplitude(value, self.look.decibel, self.look.floor)
    def alpha(self) -> float:
        overlaid = self.look.compact and self.active().waveform.channels > 1
        return COMPACT_ALPHA if overlaid else 1.0
    def windows(self, scene: Scene, lane_area: float) -> None:
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
        context = self.context
        for index, row in enumerate(roots.row):
            ms = roots.ms(row)
            x0 = self.x(min(ms))
            x1 = self.x(max(ms))
            band = layer.children[index]
            band.style.left = f"{x0 / self.ratio}px"
            band.style.width = f"{(x1 - x0) / self.ratio}px"
            band.style.height = f"{lane_area / self.ratio}px"
            if index != scene.window: continue
            context.fillStyle = self.palette.chosen
            context.fillRect(x0, 0, x1 - x0, lane_area)
    def selection(self, scene: Scene, lane_area: float) -> None:
        if not scene.selection: return
        x0 = self.x(min(scene.selection))
        x1 = self.x(max(scene.selection))
        self.context.fillStyle = self.palette.selection
        self.context.fillRect(x0, 0, x1 - x0, lane_area)
    def lanes(self, scene: Scene, width: int, lane_height: float) -> None:
        view = scene.view
        wanted = (view.t0, view.t1, width)
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
        context.lineWidth = self.ratio
        context.globalAlpha = self.alpha()
        for c in range(shown.channels):
            base, scale = self.origin(self.lane(c), lane_height)
            if self.look.half:
                peak = np.maximum(np.abs(shown.lo[c]), np.abs(shown.hi[c]))
                top = base - self.level(peak) * scale
                bottom = np.full_like(top, base)
            else:
                top = base - self.level(shown.hi[c]) * scale
                bottom = base - self.level(shown.lo[c]) * scale
            self.bars(x, top, bottom, self.palette.wave)
        context.globalAlpha = 1.0
    def bars(self, x: np.ndarray, top: np.ndarray, bottom: np.ndarray, color: str) -> None:
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
        context.strokeStyle = self.palette.wave
        context.fillStyle = self.palette.wave
        context.lineWidth = self.ratio
        context.globalAlpha = self.alpha()
        for c in range(shown.channels):
            base, scale = self.origin(self.lane(c), lane_height)
            y = base - self.level(shown.lo[c]) * scale
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
    def markers(self, scene: Scene, lane_area: float) -> None:
        roots = scene.roots
        if roots is None: return
        context = self.context
        name = roots.ms_name()
        line = (LABEL_PX + 3) * self.ratio
        context.font = f"{LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "top"
        for index, row in enumerate(roots.row):
            chosen = index == scene.window
            for field, ms in enumerate(roots.ms(row)):
                x = round(self.x(ms))
                active = chosen and field == scene.field
                context.fillStyle = self.palette.field if active else self.palette.marker
                context.fillRect(x, 0, self.ratio, lane_area)
                if not chosen: continue
                left = x + 3*self.ratio
                top = 2*self.ratio + field*line
                label = context.measureText(name[field]).width
                context.fillStyle = self.palette.background
                context.fillRect(left - 2*self.ratio, top - self.ratio, label + 4*self.ratio, line - 2*self.ratio)
                context.fillStyle = self.palette.field if active else self.palette.text
                context.fillText(name[field], left, top)
    def play(self, scene: Scene, lane_area: float) -> None:
        if scene.play_from is None or scene.play_until is None: return
        context = self.context
        height = PROGRESS_PX*self.ratio
        top = self.canvas.height - height
        position = scene.play_from if scene.playhead is None else scene.playhead
        x0 = self.x(scene.play_from)
        x1 = self.x(scene.play_until)
        head = self.x(position)
        context.fillStyle = self.palette.progress
        context.fillRect(x0, top, x1 - x0, height)
        context.fillStyle = self.palette.playhead
        context.fillRect(x0, top, head - x0, height)
        if scene.playhead is not None: context.fillRect(round(head), 0, self.ratio, lane_area)
        pad = CLOCK_PAD_PX*self.ratio
        middle = top + height / 2
        context.font = f"bold {LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "middle"
        context.fillStyle = self.palette.background
        context.textAlign = "left"
        context.fillText(millisecond(position), x0 + pad, middle)
        context.textAlign = "right"
        context.fillText(f"{millisecond(scene.play_until)} ms", x1 - pad, middle)
        context.textAlign = "left"
    def axis(self, first: int, last: int, step: float, lane_area: float) -> None:
        context = self.context
        decimals = max(0, -math.floor(math.log10(step)))
        context.fillStyle = self.palette.axis
        context.fillRect(0, lane_area, self.canvas.width, self.ratio)
        context.font = f"{LABEL_PX*self.ratio}px {self.palette.font}"
        context.textBaseline = "top"
        for k in range(first, last + 1):
            ms = k*step
            x = round(self.x(ms))
            context.fillStyle = self.palette.axis
            context.fillRect(x, lane_area, self.ratio, 5*self.ratio)
            context.fillStyle = self.palette.text
            context.fillText(f"{ms:.{decimals}f}", x + 3*self.ratio, lane_area + 6*self.ratio)
