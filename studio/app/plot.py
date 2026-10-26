from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any
import js # pyright: ignore[reportMissingImports]
import numpy as np

LABEL_PX = 10
SPARK_PX = 4
SPARK_HEIGHT = 36
FLOOR_DB = -60.0

@dataclass
class Mark:
    fraction: float
    label: str
    strong: bool

@dataclass
class Ink:
    plot: str
    bar: str
    faint: str
    dim: str
    mid: str
    soft: str
    light: str
    bright: str
    amber: str
    mark: str
    strong: str
    font: str

def ink() -> Ink:
    style = js.window.getComputedStyle(js.document.documentElement)
    def token(name: str) -> str: return style.getPropertyValue(name).strip()
    return Ink(token("--plot"), token("--bar"), token("--faint"), token("--dim"), token("--mid"), token("--soft"), token("--light"), token("--bright"), token("--amber"), token("--wave-marker"), token("--wave-strong"), token("--font-mono"))
def fit(canvas: Any) -> tuple[Any, float]:
    ratio = js.window.devicePixelRatio
    box = canvas.getBoundingClientRect()
    canvas.width = max(1, round(box.width*ratio))
    canvas.height = max(1, round(box.height*ratio))
    context = canvas.getContext("2d")
    context.clearRect(0, 0, canvas.width, canvas.height)
    return context, ratio
def peaks(value: np.ndarray, count: int) -> np.ndarray:
    edge = np.linspace(0, len(value), count + 1).astype(int)
    found = np.zeros(count)
    for i in range(count):
        i0 = min(edge[i], len(value) - 1)
        i1 = max(edge[i + 1], i0 + 1)
        found[i] = value[i0:i1].max()
    return found
def octaves(hz: np.ndarray, value: np.ndarray, z0: float, z1: float, count: int) -> np.ndarray:
    position = (np.log(hz) - math.log(z0)) / (math.log(z1) - math.log(z0))
    column = np.clip((position*count).astype(int), 0, count - 1)
    found = np.zeros(count)
    np.maximum.at(found, column, value)
    return found
def decibel(value: np.ndarray, floor: float) -> np.ndarray:
    top = float(value.max())
    with np.errstate(divide="ignore"):
        level = 10*np.log10(value / top)
    return np.clip((level - floor) / -floor, 0, 1)
def bars(context: Any, height: list[float] | np.ndarray, color: list[str] | str, box: tuple[float, float, float, float], ratio: float) -> None:
    left, top, width, tall = box
    count = len(height)
    step = width / count
    gap = ratio if step > 3*ratio else 0
    for i, share in enumerate(height):
        context.fillStyle = color if isinstance(color, str) else color[i]
        bar = max(ratio, float(share)*tall)
        context.fillRect(left + i*step, top + tall - bar, max(ratio / 2, step - gap), bar)
def rows(span: list[tuple[float, float]], pad: float) -> list[int]:
    ends: list[float] = []
    found = [0]*len(span)
    for i in sorted(range(len(span)), key=lambda i: span[i][0]):
        x0, x1 = span[i]
        row = next((r for r, end in enumerate(ends) if end + pad <= x0), len(ends))
        if row == len(ends): ends.append(x1)
        else: ends[row] = x1
        found[i] = row
    return found
def marks(context: Any, found: list[Mark], box: tuple[float, float, float, float], palette: Ink, ratio: float) -> None:
    left, top, width, tall = box
    context.font = f"{LABEL_PX*ratio}px {palette.font}"
    context.textBaseline = "top"
    place = []
    for each in found:
        x = left + round(min(max(each.fraction, 0.0), 1.0)*width)
        label = context.measureText(each.label).width
        start = x - 4*ratio - label if x + 4*ratio + label > left + width else x + 4*ratio
        place.append((x, start, label))
    row = rows([(start, start + label) for x, start, label in place], 4*ratio)
    for each, (x, start, label), r in zip(found, place, row):
        context.fillStyle = palette.strong if each.strong else palette.mark
        context.fillRect(x, top, ratio, tall)
        context.fillStyle = palette.light if each.strong else palette.mid
        context.fillText(each.label, start, top + ratio + r*(LABEL_PX + 2)*ratio)
def spark(height: np.ndarray, color: str) -> Any:
    holder = js.document.createElement("span")
    holder.className = "spark"
    canvas = js.document.createElement("canvas")
    canvas.width = len(height)*SPARK_PX
    canvas.height = SPARK_HEIGHT
    context = canvas.getContext("2d")
    bars(context, height, color, (0, 0, canvas.width, canvas.height), 1)
    holder.append(canvas)
    return holder
