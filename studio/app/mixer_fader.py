from __future__ import annotations
from collections.abc import Callable
from typing import Any
import js # pyright: ignore[reportMissingImports]
from pyodide.ffi import to_js

RANGE  = 'input[type="range"]'
DETENT = 3
NUDGE = 5
FINE = 1
FULL = 499
PAN_RANGE = 100
def position(index: int, count: int) -> float:
    if count == 1: return 0.0
    return 2*index / (count - 1) - 1
def channels(pan: int, count: int) -> list[int]:
    side = pan / PAN_RANGE
    step = []
    for i in range(count):
        gain = min(1.0, max(0.0, 1 + side*position(i, count)))
        step.append(round(FULL*gain))
    return step
def percent(value: float, slider: Any) -> str:
    v0 = float(slider.min)
    v1 = float(slider.max)
    share = (value - v0) / (v1 - v0)
    return f"{round(share*100)}%"
def balance(value: float, slider: Any) -> str:
    if value == 0: return "C"
    side = "L" if value < 0 else "R"
    return f"{side}{abs(round(value))}"
FORMAT = {"percent": percent, "pan": balance}

def paint(slider: Any) -> None:
    v0 = float(slider.min)
    v1 = float(slider.max)
    value = float(slider.value)
    share = (value - v0) / (v1 - v0) * 100
    bipolar = slider.parentElement.classList.contains("bipolar")
    filled_from = min(50, share) if bipolar else 0
    filled_to = max(50, share) if bipolar else share
    slider.style.setProperty("--a", f"{filled_from}%")
    slider.style.setProperty("--b", f"{filled_to}%")
    readout = js.document.querySelector(f'output[for="{slider.id}"]')
    format_name = slider.dataset.format
    if readout is None or format_name is None: return
    readout.textContent = FORMAT[format_name](value, slider)
def commit(slider: Any, value: str) -> None:
    slider.value = value
    bubbling = to_js({"bubbles": True}, dict_converter=js.Object.fromEntries)
    slider.dispatchEvent(js.Event.new("input", bubbling))
    slider.dispatchEvent(js.Event.new("change", bubbling))
def moved(event: Any) -> None:
    slider = event.target
    if not slider.matches(RANGE): return
    detent = slider.dataset.detent
    if slider.hasAttribute("data-detent") and detent is not None and abs(float(slider.value) - float(detent)) <= DETENT:
        slider.value = detent
    paint(slider)
def reset(event: Any) -> None:
    slider = event.target
    if not slider.matches(RANGE): return
    commit(slider, slider.defaultValue)
def nudge(event: Any) -> None:
    fader = event.target.closest(".fader")
    if fader is None: return
    slider = fader.querySelector(RANGE)
    if slider is None: return
    event.preventDefault()
    v0 = float(slider.min)
    v1 = float(slider.max)
    percent_step = FINE if event.shiftKey else NUDGE
    step = percent_step * (v1 - v0) / 100
    direction = -js.Math.sign(event.deltaY)
    value = min(v1, max(v0, float(slider.value) + direction*step))
    commit(slider, str(round(value)))
def install(on: Callable[..., None]) -> None:
    on(js.document, "input", moved, {"capture": True})
    on(js.document, "dblclick", reset)
    on(js.document, "wheel", nudge, {"passive": False})
    for slider in js.document.querySelectorAll(RANGE): paint(slider)
