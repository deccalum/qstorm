from __future__ import annotations
from typing import Any
import js # pyright: ignore[reportMissingImports]
from library import Table, Value, real, shortest

class View:
    def __init__(self, src_x0: float, src_x1: float, minimum: float) -> None:
        self.src_x0 = src_x0
        self.src_x1 = src_x1
        self.minimum = minimum
        self.x0 = src_x0
        self.x1 = src_x1
    @property
    def span(self) -> float:
        return self.x1 - self.x0
    def clamp(self) -> None:
        length = self.src_x1 - self.src_x0
        span = min(max(self.span, self.minimum), length)
        x0 = min(max(self.x0, self.src_x0), self.src_x1 - span)
        self.x0 = x0
        self.x1 = x0 + span
    def whole(self) -> None:
        self.x0 = self.src_x0
        self.x1 = self.src_x1
    def pan(self, ms: float) -> None:
        self.x0 += ms
        self.x1 += ms
        self.clamp()
    def zoom(self, factor: float, at: float) -> None:
        length = self.src_x1 - self.src_x0
        span = min(max(self.span*factor, self.minimum), length)
        scale = span/self.span
        x0 = at - (at - self.x0) * scale
        self.x0 = x0
        self.x1 = x0 + span
        self.clamp()
    def focus(self, x0: float, x1: float, margin: float) -> None:
        pad = (x1 - x0) * margin
        self.x0 = x0 - pad
        self.x1 = x1 + pad
        self.clamp()
    def ms_at(self, fraction: float) -> float:
        return self.x0 + fraction*self.span
    def fraction_of(self, ms: float) -> float:
        return (ms - self.x0) / self.span

def definitions(list_element: Any, pair: list[tuple[str, str]]) -> None:
    list_element.replaceChildren()
    for name, text in pair:
        term = js.document.createElement("dt")
        detail = js.document.createElement("dd")
        term.textContent = name
        detail.textContent = text
        list_element.append(term, detail)
def cells(table: Table, index: int, toggle: bool, unit: list[str]) -> list[tuple[str, str]]:
    pair = []
    for i, column in enumerate(table.column):
        text = table.cell(index, i, toggle, False)
        if text is not None and column.unit in unit: pair.append((column.name, text))
    return pair
def number_text(value: float) -> str:
    return f"{value:.0f}"
def fields(container: Any, table: Table, row: list[Value], names: list[str], active: int | None) -> None:
    container.replaceChildren()
    for field, column in enumerate(names):
        label = js.document.createElement("label")
        number = js.document.createElement("input")
        label.textContent = column
        label.htmlFor = f"field-{column}"
        number.id = f"field-{column}"
        number.type = "number"
        number.step = "any"
        index = table.name.index(column)
        value = row[index]
        shown = number_text if table.column[index].unit == "ms" else shortest
        number.value = "" if value is None else shown(real(value))
        number.dataset.field = str(field)
        number.dataset.column = column
        number.classList.toggle("active", field == active)
        container.append(label, number)
def span(text: str, kind: str) -> Any:
    found = js.document.createElement("span")
    found.textContent = text
    found.className = kind
    return found
def markers(container: Any, roots: Table, row: list[Value], active: int | None) -> None:
    container.replaceChildren(span("marker", "heading"), span("x ms", "heading number"),
                              span("y dB", "heading number"), span("+ ms", "heading delta"))
    x0 = min(roots.ms(row))
    for field, column in enumerate(roots.ms_name()):
        value = real(row[roots.name.index(column)])
        number = js.document.createElement("input")
        number.type = "number"
        number.step = "any"
        number.value = number_text(value)
        number.dataset.field = str(field)
        number.dataset.column = column
        number.classList.toggle("active", field == active)
        level = column.removesuffix("_x") + "_y"
        decibel = f"{real(row[roots.name.index(level)]):.1f}" if column.endswith("_x") and level in roots.name else "–"
        container.append(span(column.removesuffix("_x"), "name"), number, span(decibel, "number"),
                         span(f"+{value - x0:.0f}", "delta"))
    if "floor" in roots.name:
        container.append(span("floor", "name"), span("–", "number"),
                         span(f"{real(row[roots.name.index('floor')]):.1f}", "number"), span("", "delta"))
