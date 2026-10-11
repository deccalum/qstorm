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
    raw = await cli("read", table, *parent)
    shown = await setting(table)
    return Table.parse(raw, shown)
async def waveform(source_id: int, t0: float, t1: float, width: int) -> bytes:
    response = await pyfetch("/cli/view", method="POST", body=json.dumps([source_id, t0, t1, width]))
    if response.ok: return await response.bytes()
    failure = await response.string()
    if not js.window.confirm(f"pyramid for {source_id} failed:\n{failure}\n\nclear cache?"): raise RuntimeError(failure)
    await cli("clear_cache", source_id)
    return await waveform(source_id, t0, t1, width)
async def play(source_id: int, start: float, end: float) -> None:
    await post("/play", {"id": source_id, "start": start, "end": end})
async def mixer(master: int, channel: list[int], mono: bool) -> None:
    await post("/mixer", {"master": master, "channel": channel, "mono": mono})
async def level() -> list[float]:
    raw = await post("/level", {})
    return json.loads(raw)
async def pause() -> None:
    await post("/pause", {})
async def stop() -> None:
    await post("/stop", {})
