from __future__ import annotations
from typing import Any
import js # pyright: ignore[reportMissingImports]
from library import Table, Value, real

class View:
    def __init__(self, start: float, end: float, minimum: float) -> None:
        self.start = start
        self.end = end
        self.minimum = minimum
        self.t0 = start
        self.t1 = end
    @property
    def span(self) -> float:
        return self.t1 - self.t0
    def clamp(self) -> None:
        length = self.end - self.start
        span = min(max(self.span, self.minimum), length)
        t0 = min(max(self.t0, self.start), self.end - span)
        self.t0 = t0
        self.t1 = t0 + span
    def whole(self) -> None:
        self.t0 = self.start
        self.t1 = self.end
    def pan(self, ms: float) -> None:
        self.t0 += ms
        self.t1 += ms
        self.clamp()
    def zoom(self, factor: float, at: float) -> None:
        length = self.end - self.start
        span = min(max(self.span*factor, self.minimum), length)
        scale = span/self.span
        t0 = at - (at - self.t0) * scale
        self.t0 = t0
        self.t1 = t0 + span
        self.clamp()
    def focus(self, lo: float, hi: float, margin: float) -> None:
        pad = (hi - lo) * margin
        self.t0 = lo - pad
        self.t1 = hi + pad
        self.clamp()
    def ms_at(self, fraction: float) -> float:
        return self.t0 + fraction*self.span
    def fraction_of(self, ms: float) -> float:
        return (ms - self.t0) / self.span

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
def fields(container: Any, roots: Table, row: list[Value], active: int | None) -> None:
    container.replaceChildren()
    for field, column in enumerate(roots.ms_name()):
        label = js.document.createElement("label")
        number = js.document.createElement("input")
        label.textContent = column
        label.htmlFor = f"field-{column}"
        number.id = f"field-{column}"
        number.type = "number"
        number.step = "any"
        number.value = number_text(real(row[roots.name.index(column)]))
        number.dataset.field = str(field)
        number.dataset.column = column
        number.classList.toggle("active", field == active)
        container.append(label, number)
