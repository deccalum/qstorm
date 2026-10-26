from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any
import js # pyright: ignore[reportMissingImports]
import numpy as np
from library import Table, Value, real
from plot import FLOOR_DB, Ink, Mark, bars, decibel, fit, marks, octaves, peaks, spark
from relay import read
from view import definitions
from waveform import Atomization, Stft

SPARK_COLUMNS = 28
DETAIL_COLUMNS = 72
TRACE_COLUMNS = 60

Record = dict[str, Value]

@dataclass
class Component:
    run: Record
    part: dict[int, Record]
    order: list[int]
    spectrum_set: dict[int, int]
    event_set: dict[int, int]
    pair: dict[tuple[int, int], tuple[float, float]]
    related: int

def record(table: Table, row: list[Value]) -> Record:
    return dict(zip(table.name, row))
def latest(table: Table) -> int | None:
    made = [table.key_cells(row)[0] for row in table.row]
    return max(made, default=None)
async def fetch(atomization: Atomization, extraction: int) -> Component:
    runs = await read("atomization", extraction)
    run = next(record(runs, row) for row in runs.row if row[runs.name.index("id")] == atomization.id)
    parts = await read("component", atomization.id)
    part = {int(real(row[parts.name.index("index")])): record(parts, row) for row in parts.row}
    order = [int(k) for k in atomization.index if part[int(k)]["relevant"]]
    spectrum_set: dict[int, int] = {}
    event_set: dict[int, int] = {}
    pair: dict[tuple[int, int], tuple[float, float]] = {}
    related = 0
    groupings = await read("grouping", atomization.id)
    grouping = latest(groupings)
    if grouping is not None:
        groups = await read("component_group", grouping)
        for row in groups.row:
            each = record(groups, row)
            k = int(real(each["component"]))
            spectrum_set[k] = int(real(each["spectrum"]))
            event_set[k] = int(real(each["event"]))
        pairs = await read("component_pair", grouping)
        for row in pairs.row:
            each = record(pairs, row)
            key = (int(real(each["component"])), int(real(each["other"])))
            pair[key] = (real(each["cosine"]), real(each["correlation"]))
            related += 1 if each["related"] else 0
    return Component(run, part, order, spectrum_set, event_set, pair, related)
def hz(value: float) -> str:
    if value >= 1: return f"{value:.2f} kHz"
    return f"{value*1000:.0f} Hz"
def ms(value: float) -> str:
    return f"{value:.1f} ms"
def signif(value: float) -> str:
    if math.isinf(value): return "∞"
    return f"{value:.3g}"
def text(value: Value, shown: Any) -> str:
    return "-" if value is None else shown(real(value))
def spectrum_shape(stft: Stft, spectrum: list[float], count: int) -> np.ndarray:
    frequency = stft.hz()
    width = np.concatenate([np.full(each.db.shape[1], 1 / each.density) for each in stft.band])
    octave = np.log2((frequency + width / 2) / np.maximum(frequency - width / 2, 1e-9))
    density = np.asarray(spectrum) / octave
    return decibel(octaves(frequency, density, stft.z0, stft.z1, count), FLOOR_DB)
def activation_shape(activation: list[float], count: int) -> np.ndarray:
    column = peaks(np.asarray(activation), count)
    return column / max(float(column.max()), 1e-300)
def fraction_hz(stft: Stft, value: float) -> float:
    return (math.log(max(value, stft.z0)) - math.log(stft.z0)) / (math.log(stft.z1) - math.log(stft.z0))
def set_label(report: Component, k: int) -> str:
    found = []
    event = report.event_set.get(k)
    spectrum = report.spectrum_set.get(k)
    if event is not None and sum(1 for each in report.event_set.values() if each == event) > 1: found.append(f"e{event}")
    if spectrum is not None and sum(1 for each in report.spectrum_set.values() if each == spectrum) > 1: found.append(f"s{spectrum}")
    return " ".join(found) if found else "-"
def members(sets: dict[int, int]) -> dict[int, list[int]]:
    found: dict[int, list[int]] = {}
    for k, owner in sets.items(): found.setdefault(owner, []).append(k)
    return {owner: sorted(each) for owner, each in found.items() if len(each) > 1}
def summary(report: Component, stft: Stft, palette: Ink) -> None:
    run = report.run
    share = sum(real(report.part[k]["adsr_energy"]) for k in report.order)
    definitions(js.document.getElementById("run"), [ ("iterations", str(int(real(run["iteration"])))), ("divergence", signif(real(run["divergence"]))), ("components", f"{len(report.order)} of {int(real(run['factor_K']))}"), ("share", f"{share*100:.1f} %")])
    trace = np.asarray(run["trace"], float)
    level = np.log10(peaks(trace, min(TRACE_COLUMNS, len(trace))))
    y0 = float(level.min())
    span = max(float(level.max()) - y0, 1e-12)
    height = 0.06 + 0.94*(level - y0) / span
    canvas = js.document.getElementById("trace")
    context, ratio = fit(canvas)
    color = [palette.bar]*(len(height) - 1) + [palette.bright]
    bars(context, height, color, (0, 0, canvas.width, canvas.height), ratio)
    axis = js.document.getElementById("trace-axis")
    axis.replaceChildren()
    for label in ("1", str(len(trace) - 1)):
        each = js.document.createElement("span")
        each.textContent = label
        axis.append(each)
    residual = [real(each) for each in run["residual"]]
    top = max(residual)
    table = js.document.getElementById("band-residual")
    table.replaceChildren()
    for band, value in zip(stft.band, residual):
        upper = band.offset + band.db.shape[1] / band.density
        name = js.document.createElement("span")
        name.className = "name"
        name.textContent = f"{band.offset*1000:.0f}–{upper*1000:.0f}"
        name.title = f"{hz(band.offset)} – {hz(upper)}"
        track = js.document.createElement("div")
        track.className = "track"
        fill = js.document.createElement("span")
        fill.style.width = f"{value / top*100:.1f}%"
        track.append(fill)
        shown = js.document.createElement("span")
        shown.className = "value"
        shown.textContent = signif(value)
        table.append(name, track, shown)
    definitions(js.document.getElementById("shape-facts"), [ ("event sets", f"{len(members(report.event_set))} · {report.related} related pairs"), ("spectrum sets", str(len(members(report.spectrum_set))))])
def cell(content: str, kind: str = "") -> Any:
    found = js.document.createElement("span")
    found.textContent = content
    if kind: found.className = kind
    return found
def chip(k: int, detail: str, chosen: int | None) -> Any:
    button = js.document.createElement("button")
    button.type = "button"
    button.className = "chip active" if k == chosen else "chip"
    button.dataset.component = str(k)
    button.append(f"k{k} ")
    button.append(cell(detail))
    return button
def components(container: Any, report: Component, stft: Stft, chosen: int | None, palette: Ink) -> None:
    container.replaceChildren()
    head = js.document.createElement("div")
    head.className = "component-head"
    spacer = js.document.createElement("span")
    spacer.className = "spacer"
    grid = js.document.createElement("div")
    grid.className = "component-grid"
    for name in ("k", "λ / floor", "peak", "onset", "decay", "share", "ΔIS", "set", "W", "H"): grid.append(cell(name))
    head.append(spacer, grid)
    container.append(head)
    floor = real(report.run["floor"])
    for k in report.order:
        each = report.part[k]
        active = k == chosen
        color = palette.amber if active else palette.bar
        line = js.document.createElement("div")
        line.className = "component-line active" if active else "component-line"
        play = js.document.createElement("button")
        play.type = "button"
        play.className = "key"
        play.textContent = "▶"
        play.title = "play component"
        play.dataset.component = str(k)
        play.dataset.play = "1"
        row = js.document.createElement("button")
        row.type = "button"
        row.className = "row component-grid"
        row.dataset.component = str(k)
        row.append( cell(f"k{k}", "k"), cell(f"{real(each['lambda']) / floor:.1f}"), cell(text(each["spectral_peak"], hz)), cell(text(each["adsr_x0"], ms)), cell(text(each["adsr_decay"], ms)), cell(f"{real(each['adsr_energy'])*100:.1f}%"), cell("+" + signif(real(each["divergence"]))), cell(set_label(report, k), "set"), spark(spectrum_shape(stft, list(each["spectrum"]), SPARK_COLUMNS), color), spark(activation_shape(list(each["activation"]), SPARK_COLUMNS), color))
        line.append(play, row)
        container.append(line)
def basis(report: Component, group: list[int], index: int, name: str) -> str:
    found = []
    for i, a in enumerate(group):
        for b in group[i + 1:]:
            key = (a, b) if (a, b) in report.pair else (b, a)
            found.append(report.pair[key][index])
    return f"{name} {min(found):.2f}–{max(found):.2f}" if len(found) > 1 else f"{name} {found[0]:.2f}"
def groups(container: Any, report: Component, chosen: int | None) -> None:
    container.replaceChildren()
    if not report.event_set:
        container.append(cell("no grouping run; analyze runs one", "empty"))
        return
    head = js.document.createElement("div")
    head.className = "set-grid heading"
    head.append(cell("set"), cell("components"), cell("basis"), cell("share", "share"))
    container.append(head)
    for kind, sets, index, name in (("e", report.event_set, 1, "H corr"), ("s", report.spectrum_set, 0, "W cos")):
        for owner, group in members(sets).items():
            row = js.document.createElement("div")
            row.className = "set-grid"
            chips = js.document.createElement("div")
            chips.className = "chips"
            for k in group: chips.append(chip(k, text(report.part[k]["spectral_peak"], hz), chosen))
            share = sum(real(report.part[k]["adsr_energy"]) for k in group)
            row.append(cell(f"{kind}{owner}", "label"), chips, cell(basis(report, group, index, name), "basis"), cell(f"{share*100:.1f}%", "share"))
            container.append(row)
def section(title: str, unit: str, pair: list[tuple[str, str]]) -> Any:
    holder = js.document.createElement("div")
    heading = js.document.createElement("h3")
    heading.append(title + " ")
    heading.append(cell(unit))
    listing = js.document.createElement("dl")
    listing.className = "pairs"
    definitions(listing, pair)
    holder.append(heading, listing)
    return holder
def detail(report: Component, k: int | None, atomization: Atomization, stft: Stft, palette: Ink) -> None:
    title = js.document.getElementById("detail-title")
    sub = js.document.getElementById("detail-sub")
    sections = js.document.getElementById("detail-sections")
    solo = js.document.getElementById("solo")
    w_canvas = js.document.getElementById("detail-w")
    h_canvas = js.document.getElementById("detail-h")
    if k is None:
        title.textContent = "component"
        sub.textContent = "pick a component"
        sections.replaceChildren()
        solo.disabled = True
        fit(w_canvas)
        fit(h_canvas)
        return
    each = report.part[k]
    floor = real(report.run["floor"])
    title.textContent = f"component {k}"
    sub.textContent = f"λ / floor {real(each['lambda']) / floor:.1f} · {set_label(report, k)}"
    solo.disabled = False
    context, ratio = fit(w_canvas)
    top = 14*ratio
    box = (2*ratio, top, w_canvas.width - 4*ratio, w_canvas.height - top)
    bars(context, spectrum_shape(stft, list(each["spectrum"]), DETAIL_COLUMNS), palette.amber, box, ratio)
    edge = [("lo", "spectral_low"), ("peak", "spectral_peak"), ("hi", "spectral_high")]
    marks(context, [Mark(fraction_hz(stft, real(each[name])), label, label == "peak") for label, name in edge if each[name] is not None], (0, 0, w_canvas.width, w_canvas.height), palette, ratio)
    js.document.getElementById("detail-w-range").textContent = f"{hz(stft.z0)} – {hz(stft.z1)}"
    context, ratio = fit(h_canvas)
    top = 14*ratio
    box = (2*ratio, top, h_canvas.width - 4*ratio, h_canvas.height - top)
    bars(context, activation_shape(list(each["activation"]), DETAIL_COLUMNS), palette.amber, box, ratio)
    length = atomization.x1 - atomization.x0
    time = [("on", "adsr_x0"), ("pk", "adsr_peak"), ("end", "adsr_x1")]
    marks(context, [Mark((real(each[name]) - atomization.x0) / length, label, label == "pk") for label, name in time if each[name] is not None], (0, 0, h_canvas.width, h_canvas.height), palette, ratio)
    js.document.getElementById("detail-h-range").textContent = f"{atomization.x0:.0f} – {atomization.x1:.0f} ms"
    spectral = [("peak", text(each["spectral_peak"], hz)), ("center", text(each["spectral_center"], hz)), ("spread", text(each["spectral_spread"], hz)), ("low", text(each["spectral_low"], hz)), ("high", text(each["spectral_high"], hz)), ("tilt", text(each["spectral_tilt"], lambda v: f"{v:.1f} dB/oct")), ("flatness", text(each["spectral_flatness"], signif))]
    temporal = [("onset", text(each["adsr_x0"], ms)), ("attack", text(each["adsr_attack"], ms)), ("peak", text(each["adsr_peak"], ms)), ("decay", text(each["adsr_decay"], ms)), ("sustain", text(each["adsr_sustain"], lambda v: f"{v:.1f} dB")), ("release", text(each["adsr_release"], ms)), ("end", text(each["adsr_x1"], ms)), ("slope", text(each["adsr_slope"], lambda v: f"{v:.4f} dB/ms")), ("center", text(each["adsr_center"], ms))]
    relations = [("energy share", f"{real(each['adsr_energy'])*100:.1f} %"), ("ΔIS leave-one-out", "+" + signif(real(each["divergence"])))]
    event = report.event_set.get(k)
    spectrum = report.spectrum_set.get(k)
    if event is not None: relations.append(("event set", " ".join(f"k{m}" for m in sorted(report.event_set) if report.event_set[m] == event and m != k) or "-"))
    if spectrum is not None: relations.append(("spectrum set", " ".join(f"k{m}" for m in sorted(report.spectrum_set) if report.spectrum_set[m] == spectrum and m != k) or "-"))
    sections.replaceChildren(section("spectral", "W", spectral), section("temporal", "H · adsr", temporal), section("relations", "", relations))
