"""The "MMseqs server": keeps the TED lookup DB hot so each query takes seconds instead of minutes.

MMseqs2 has no query daemon of its own for CPU search; the standard setup (used e.g. by the ColabFold MSA
server) is a pre-computed index that stays in memory between searches:
  CPU:  mmseqs createindex (once) -> mmseqs touchdb (loads it into the page cache) -> searches with
        --db-load-mode 2 reuse the cached index instead of reading ~100+ GB from disk each time.
  GPU:  mmseqs makepaddedseqdb (once) -> `mmseqs gpuserver` holds the DB in GPU memory -> searches with
        --gpu 1 --gpu-server 1 hand their queries to it.
The page cache / GPU memory belong to one node, so the searches have to run on that node. This module runs
a small HTTP service inside the Slurm job that holds the DB (scripts/mmseqs_server.sbatch):

    POST /search  {"sequences": {"name": "SEQ", ...}}  ->  {"hits": {"name": [hit, ...]}}
    GET  /health

and writes its address to db/mmseqs/SERVER, which ted_recreate.classify picks up automatically.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import socket
import sys
import subprocess
import tempfile
import threading
import time
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from . import config
from .classify import mmseqs_search

SERVER_FILE = config.TED_MMSEQS_DB.parent / "SERVER"


def server_url() -> str | None:
    """Address of a running server: $TED_MMSEQS_SERVER or db/mmseqs/SERVER, if it answers /health."""
    url = os.environ.get("TED_MMSEQS_SERVER")
    if not url and SERVER_FILE.exists():
        url = json.loads(SERVER_FILE.read_text()).get("url")
    if not url:
        return None
    try:
        with urllib.request.urlopen(f"{url}/health", timeout=5) as r:
            return url if r.status == 200 else None
    except OSError:
        return None


def remote_search(url: str, seqs: dict, timeout: int = 3600) -> dict:
    req = urllib.request.Request(f"{url}/search", data=json.dumps({"sequences": seqs}).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)["hits"]


class _Handler(BaseHTTPRequestHandler):
    opts: dict = {}
    lock = threading.Semaphore(int(os.environ.get("TED_SERVER_CONCURRENCY", "2")))

    def _send(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/health":
            self._send(200, {"ok": True, "mode": self.opts["mode"], "db": str(self.opts["db"])})
        else:
            self._send(404, {"error": "not found"})

    def do_POST(self):
        if self.path != "/search":
            return self._send(404, {"error": "not found"})
        try:
            seqs = json.loads(self.rfile.read(int(self.headers["Content-Length"])))["sequences"]
            work = Path(tempfile.mkdtemp(prefix="req_", dir=self.opts["tmp"]))
            with self.lock:
                t = time.time()
                hits = mmseqs_search(seqs, work, threads=self.opts["threads"], db=self.opts["db"], db_load_mode=2,
                                     gpu=self.opts["mode"] == "gpu")
            shutil.rmtree(work, ignore_errors=True)
            self._send(200, {"hits": hits, "seconds": round(time.time() - t, 2)})
        except Exception as e:  # report errors to the client instead of dropping the connection
            self._send(500, {"error": repr(e)[:2000]})

    def log_message(self, fmt, *args):
        print(f"[{time.strftime('%H:%M:%S')}] {self.address_string()} {fmt % args}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--port", type=int, default=int(os.environ.get("TED_SERVER_PORT", "8765")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    ap.add_argument("--db", default=str(config.TED_MMSEQS_DB))
    args = ap.parse_args()
    db = Path(args.db)
    tmp = Path(os.environ.get("TMPDIR", tempfile.gettempdir())) / f"ted_mmseqs_server_{os.getpid()}"
    tmp.mkdir(parents=True, exist_ok=True)

    gpu_proc = None
    if args.mode == "cpu":
        idx = Path(str(db) + ".idx")
        target = idx if idx.exists() else db
        print(f"preloading {target} into the page cache (touchdb) ...", flush=True)
        subprocess.run([str(config.MMSEQS), "touchdb", str(target), "--threads", str(args.threads)], check=True)
    else:
        os.environ["TED_MMSEQS_GPU_SERVER"] = "1"
        pad = Path(str(db) + "_pad")
        print(f"starting mmseqs gpuserver on {pad} ...", flush=True)
        gpu_proc = subprocess.Popen([str(config.MMSEQS_GPU), "gpuserver", str(pad), "--max-seqs", "300",
                                     "--db-load-mode", "2"])
        time.sleep(30)

    _Handler.opts = {"mode": args.mode, "db": db, "threads": args.threads, "tmp": tmp}
    httpd = ThreadingHTTPServer(("0.0.0.0", args.port), _Handler)
    url = f"http://{socket.getfqdn()}:{args.port}"
    server_file = db.parent / "SERVER"
    server_file.write_text(json.dumps({"url": url, "mode": args.mode, "slurm_job": os.environ.get("SLURM_JOB_ID"),
                                       "started": time.strftime("%Y-%m-%d %H:%M:%S")}))
    print(f"TED MMseqs server ({args.mode}) listening at {url}", flush=True)
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))   # slurm stops jobs with SIGTERM; clean up
    try:
        httpd.serve_forever()
    finally:
        if server_file.exists() and json.loads(server_file.read_text()).get("url") == url:
            server_file.unlink()
        if gpu_proc:
            gpu_proc.terminate()


if __name__ == "__main__":
    main()
