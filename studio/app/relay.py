from __future__ import annotations
import json
import js # pyright: ignore[reportMissingImports]
from pyodide.http import pyfetch
from library import Setting, Table

Argument = str | int | float
async def post(path: str, body: object) -> bytes:
    response = await pyfetch(path, method="POST", body=json.dumps(body))
    return await response.bytes()
async def cli(verb: str, *argument: Argument) -> bytes:
    return await post(f"/cli/{verb}", list(argument))
async def setting(table: str) -> dict[str, Setting]:
    raw = await post("/setting", {"table": table})
    found = {}
    for row in json.loads(raw): found[row["column"]] = Setting(row["code"], row["token"])
    return found
async def read(table: str, *parent: Argument) -> Table:
    raw = await cli("read_table", table, *parent)
    shown = await setting(table)
    return Table.parse(raw, shown)
async def clear_cache(failure: str, table: str, parent: int) -> None:
    if not js.window.confirm(f"{table} for {parent} failed:\n{failure}\n\nclear cache?"): raise RuntimeError(failure)
    await cli("clear_cache", table, parent)
async def waveform(source_id: int, x0: float, x1: float, width: int) -> bytes:
    response = await pyfetch("/cli/render_waveform", method="POST", body=json.dumps([source_id, x0, x1, width]))
    if response.ok: return await response.bytes()
    failure = await response.string()
    await clear_cache(failure, "decoded", source_id)
    await cli("decode_source", source_id)
    return await waveform(source_id, x0, x1, width)
async def read_roots(source_id: int) -> Table:
    extraction = await read("extraction", source_id)
    every = await read("root")
    made = [extraction.key_cells(row)[0] for row in extraction.row]
    newest = max(made, default=None)
    key = every.name.index("extraction")
    row = [each for each in every.row if each[key] == newest]
    return Table(every.column, row, every.setting)
async def stft(extraction: int, index: int, revision: int) -> bytes:
    stored = await read("stft", extraction)
    written = [stored.key_cells(row) for row in stored.row]
    if [extraction, index, revision] not in written: await cli("write_stft", extraction, index)
    response = await pyfetch("/cli/render_stft", method="POST", body=json.dumps([extraction, index]))
    if response.ok: return await response.bytes()
    failure = await response.string()
    await clear_cache(failure, "stft", extraction)
    return await stft(extraction, index, revision)
async def play(audio: int, offset: float) -> None:
    await post("/play", {"audio": audio, "offset": offset})
async def render_component_list(atomization: int, width: int) -> bytes:
    return await cli("render_component_list", atomization, width)
async def mixer(master: int, channel: list[int], mono: bool) -> None:
    await post("/mixer", {"master": master, "channel": channel, "mono": mono})
async def level() -> list[float]:
    raw = await post("/level", {})
    return json.loads(raw)
async def pause() -> None:
    await post("/pause", {})
async def stop() -> None:
    await post("/stop", {})
