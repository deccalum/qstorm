from __future__ import annotations
import json
import os
import subprocess
import threading
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import psycopg2

REPO_ROOT = Path(__file__).resolve().parents[1]
STUDIO = Path(__file__).resolve().parent
BINARY = REPO_ROOT/"build"/"qstorm"

def command(*argument) -> bytes:
    return subprocess.run([str(BINARY), *map(str, argument)], cwd=REPO_ROOT, stdout=subprocess.PIPE,check=True).stdout
def setting(table: str) -> list[dict]:
    with psycopg2.connect(os.environ["DATABASE_URL"]) as connection, connection.cursor() as cursor:
        cursor.execute('SELECT "column", code, token FROM column_setting WHERE "table" = %s', (table,))
        rows = cursor.fetchall()
    found = []
    for column, code, token in rows: found.append({"column": column, "code": code, "token": token})
    return found

class Player:
    def __init__(self) -> None:
        self.process: subprocess.Popen[bytes] | None = None
        self.master = 499
        self.channel: list[int] = []
        self.mono = False
        self.level: list[float] = []
    def play(self, *argument: object) -> None:
        self.stop()
        self.process = subprocess.Popen([str(BINARY), *map(str, argument)], cwd=REPO_ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        reader = threading.Thread(target=self.listen, args=(self.process,), daemon=True)
        reader.start()
        self.send_mixer()
    def listen(self, process: subprocess.Popen[bytes]) -> None:
        stdout = process.stdout
        if stdout is None: return
        for line in stdout:
            word = line.decode().split()
            if word[0] != "level" or process is not self.process: continue
            self.level = [float(each) for each in word[1:]]
    def order(self, word: str) -> None:
        process = self.process
        if process is None or process.poll() is not None: return
        stdin = process.stdin
        if stdin is None: return
        stdin.write(f"{word}\n".encode())
        stdin.flush()
    def send_mixer(self) -> None:
        self.order(f"master {self.master}")
        for c, step in enumerate(self.channel): self.order(f"volume {c} {step}")
        self.order(f"mono {int(self.mono)}")
    def mixer(self, master: int, channel: list[int], mono: bool) -> None:
        self.master = master
        self.channel = channel
        self.mono = mono
        self.send_mixer()
    def pause(self) -> None:
        self.order("pause")
    def stop(self) -> None:
        self.order("stop")
        process = self.process
        if process is not None and process.stdin is not None:
            process.stdin.close()
            process.wait()
        self.process = None
        self.level = []
PLAYER = Player()

class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, directory=str(STUDIO), **kwargs)
    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store")
        super().end_headers()
    def send_bytes(self, body: bytes) -> None:
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
    def send_failure(self, failed: subprocess.CalledProcessError) -> None:
        body = str(failed).encode()
        self.send_response(500)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)
    def do_POST(self) -> None:
        body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        route, _, verb = self.path.removeprefix("/").partition("/")
        if route == "cli":
            try: self.send_bytes(command(verb, *body))
            except subprocess.CalledProcessError as failed: self.send_failure(failed)
        if route == "setting":
            self.send_bytes(json.dumps(setting(body["table"])).encode())
        if route == "play":
            PLAYER.play("play_audio", body["audio"], body["offset"])
            self.send_bytes(b"")
        if route == "mixer":
            PLAYER.mixer(body["master"], body["channel"], body["mono"])
            self.send_bytes(b"")
        if route == "level":
            self.send_bytes(json.dumps(PLAYER.level).encode())
        if route == "pause":
            PLAYER.pause()
            self.send_bytes(b"")
        if route == "stop":
            PLAYER.stop()
            self.send_bytes(b"")

def run() -> None:
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    host, port = server.server_address[:2]
    url = f"http://{host}:{port}/"
    print(url, flush=True)
    webbrowser.open(url)
    try: server.serve_forever()
    finally: PLAYER.stop()
