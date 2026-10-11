from __future__ import annotations
import struct
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Literal, overload
import js # pyright: ignore[reportMissingImports]

Value = int | float | bool | str | bytes | datetime | list[int] | list[float] | None

INTEGER = ("int4", "int8", "oid")
DATE = ("d", "m")
MS_PER_SECOND = 1000

def integer(value: Value) -> int:
    if isinstance(value, int) and not isinstance(value, bool): return value
    raise TypeError(type(value).__name__)
def real(value: Value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)): raise TypeError(type(value).__name__)
    return float(value)
def line(texts: list[str], kind: str) -> Any:
    row = js.document.createElement("button" if kind != "heading" else "div")
    row.className = f"line {kind}"
    for text in texts:
        cell = js.document.createElement("span")
        cell.textContent = text
        row.append(cell)
    return row
def heading(table: Table) -> Any:
    names = []
    for i, column in enumerate(table.column):
        if table.shown(i, False): names.append(column.name)
    return line(names, "heading")
def entry(table: Table, index: int, second: bool, kind: str) -> Any:
    texts = []
    for i in range(len(table.column)):
        text = table.cell(index, i, False, second)
        if text is not None: texts.append(text)
    return line(texts, kind)
def listing(list_element: Any, sources: Table, windows: dict[int, Table], second: bool, active: int | None, chosen: int | None) -> None:
    top = heading(sources)
    list_element.style.setProperty("--columns", str(top.childElementCount))
    list_element.replaceChildren(top)
    for i in range(len(sources.row)):
        source_id = sources.key_cells(sources.row[i])[0]
        rec = entry(sources, i, second, "rec")
        rec.dataset.source = str(i)
        rec.classList.toggle("active", i == active)
        list_element.append(rec)
        roots = windows.get(source_id)
        if roots is None: continue
        inner = js.document.createElement("div")
        inner.className = "windows"
        names = heading(roots)
        inner.style.setProperty("--columns", str(names.childElementCount))
        inner.append(names)
        rejected = roots.name.index("rejected")
        for j in range(len(roots.row)):
            window = entry(roots, j, second, "window")
            window.classList.toggle("rejected", bool(roots.row[j][rejected]))
            window.dataset.source = str(i)
            window.dataset.window = str(j)
            window.classList.toggle("active", i == active and j == chosen)
            inner.append(window)
        list_element.append(inner)

@dataclass
class Column:
    name: str
    type: str
    key: bool
    reference: bool
    unit: str

@dataclass
class Setting:
    code: list[float]
    token: str
PLAIN = Setting([1], "")

@dataclass
class Appearance:
    visible: bool
    toggle: bool
    round: int | None

def appearance(code: float) -> Appearance:
    whole = int(code)
    amount = round((code - whole) * 10)
    if whole == 0: return Appearance(False, False, None)
    if whole == 1: return Appearance(True, False, None)
    if whole == 2: return Appearance(False, True, None)
    return Appearance(True, False, amount)
def included(code: float, toggle: bool) -> bool:
    look = appearance(code)
    return look.visible or (toggle and look.toggle)
def shortest(value: float) -> str:
    if float(value).is_integer(): return str(int(value))
    return repr(float(value))
def significant(value: float) -> str:
    return shortest(float(f"{value:.5e}"))
def stamp(value: datetime, token: str) -> str:
    if token == "d": return f"{value.day:02}"
    if token == "m": return f"{value.month:02}"
    return f"{value.hour:02}"
def clock(value: datetime, setting: Setting, toggle: bool) -> str:
    text = ""
    previous = ""
    for code, token in zip(setting.code, setting.token):
        if not included(code, toggle) or token == previous: continue
        if previous: text += "." if previous in DATE and token in DATE else " "
        text += stamp(value, token)
        previous = token
    return text
def plain(value: Value) -> str:
    if isinstance(value, bool): return "true" if value else "false"
    if isinstance(value, str): return value
    if isinstance(value, bytes): return f"{len(value)} B"
    if isinstance(value, list):
        each = []
        for number in value: each.append(shortest(number))
        return ",".join(each)
    if isinstance(value, int): return str(value)
    if not isinstance(value, float): raise TypeError(type(value).__name__)
    return shortest(value)
def timed(column: Column, setting: Setting) -> bool:
    return column.unit == "ms" or setting.token == "t"
def number(value: float, column: Column, setting: Setting, second: bool) -> str:
    places = appearance(setting.code[0]).round
    if places is not None: return f"{value:.{places}f}"
    if timed(column, setting) and second: return f"{value / MS_PER_SECOND:.1f}"
    if timed(column, setting): return f"{value:.0f}"
    if column.unit == "db": return f"{value:.1f} dB"
    if column.type == "float8": return significant(value)
    return plain(value)
def present(value: Value, column: Column, setting: Setting, toggle: bool, second: bool) -> str:
    if value is None: return "-"
    if isinstance(value, datetime): return clock(value, setting, toggle)
    if isinstance(value, (bool, str, bytes, list)): return plain(value)
    return number(value, column, setting, second)

class Reader:
    def __init__(self, buffer: bytes) -> None:
        self.buffer = buffer
        self.at = 0
    @overload
    def unpack(self, form: Literal["<B", "<I", "<q"]) -> int: ...
    @overload
    def unpack(self, form: Literal["<d"]) -> float: ...
    def unpack(self, form: str) -> int | float:
        value = struct.unpack_from(form, self.buffer, self.at)[0]
        self.at += struct.calcsize(form)
        return value
    def   u8(self) -> int:   return self.unpack("<B")
    def  u32(self) -> int:   return self.unpack("<I")
    def  i64(self) -> int:   return self.unpack("<q")
    def  f64(self) -> float: return self.unpack("<d")
    def blob(self) -> bytes:
        size = self.u32()
        raw = self.buffer[self.at:self.at + size]
        self.at += size
        return bytes(raw)
    def text(self) -> str:
        return self.blob().decode()
    def value(self, sql_type: str) -> Value:
        if not self.u8(): return None
        if sql_type == "bool": return self.u8() == 1
        if sql_type == "text": return self.text()
        if sql_type == "bytea": return self.blob()
        if sql_type == "timestamptz": return datetime.fromtimestamp(self.i64() / 1e6)
        if sql_type.endswith("[]"):
            count = self.u32()
            each = []
            integer = sql_type[:-2] in INTEGER
            for _ in range(count): each.append(self.i64() if integer else self.f64())
            return each
        if sql_type in INTEGER: return self.i64()
        return self.f64()
    def column(self) -> Column:
        name = self.text()
        sql_type = self.text()
        key = self.u8() == 1
        reference = self.u8() == 1
        unit = self.text()
        return Column(name, sql_type, key, reference, unit)

class Table:
    def __init__(self, column: list[Column], row: list[list[Value]], setting: list[Setting]) -> None:
        self.column = column
        self.row = row
        self.setting = setting
    @staticmethod
    def parse(buffer: bytes, shown: dict[str, Setting]) -> Table:
        reader = Reader(buffer)
        column_count = reader.u32()
        column = []
        for _ in range(column_count): column.append(reader.column())
        row_count = reader.u32()
        row = []
        for _ in range(row_count):
            values = []
            for each in column: values.append(reader.value(each.type))
            row.append(values)
        setting = []
        for each in column: setting.append(shown.get(each.name, PLAIN))
        return Table(column, row, setting)
    @property
    def name(self) -> list[str]:
        names = []
        for each in self.column: names.append(each.name)
        return names
    def key_cells(self, row: list[Value]) -> list[int]:
        keys = []
        for i, value in enumerate(row):
            if self.column[i].key: keys.append(integer(value))
        return keys
    def ms(self, row: list[Value]) -> list[float]:
        times = []
        for i, value in enumerate(row):
            if self.column[i].unit == "ms" and self.shown(i, True): times.append(real(value))
        return times
    def ms_name(self) -> list[str]:
        names = []
        for i, each in enumerate(self.column):
            if each.unit == "ms" and self.shown(i, True): names.append(each.name)
        return names
    def units(self) -> list[str]:
        found = []
        for each in self.column:
            if each.unit not in found: found.append(each.unit)
        return found
    def shown(self, field: int, toggle: bool) -> bool:
        for code in self.setting[field].code:
            if included(code, toggle): return True
        return False
    def cell(self, index: int, field: int, toggle: bool, second: bool) -> str | None:
        if not self.shown(field, toggle): return None
        return present(self.row[index][field], self.column[field], self.setting[field], toggle, second)
