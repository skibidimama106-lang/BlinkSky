"""
BlinkSky helper server: serves the game and fetches SPHEREx photos for any coordinates.

What it does : A web page is not allowed to call NASA/IPAC IRSA directly (IRSA sends no cross-site permission
               header), so this small server does it on the page's behalf:
                 GET /                       the BlinkSky files (same folder as index.html)
                 GET /api/explore?ra=&dec=   runs blinksky_explorer for that spot (quick mode: 1 photo per detector
                                             per visit), saves it under explore/_live/, answers {"ok", "dir"}
                 GET /api/progress?key=      {"text": "..."} so the page can show what is happening
Inputs       : --port N       (default 8831)
               --root DIR     folder to serve (default: the folder above this file)
               --page NAME    file to open at "/" (default: index.html)
               --cache DIR    cutout cache (default .tmp/blinksky_explorer_cache)
Outputs      : Runs until stopped. Requests are logged to stderr.
Needs        : Python 3.10+ and everything blinksky_explorer.py needs.

Example:
    python blinksky_server.py --port 8831 --page index.v5.claude.html
"""

import argparse
import json
import sys
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

sys.path.insert(0, str(Path(__file__).parent))
from blinksky_explorer import explore   # noqa: E402

PROGRESS: dict[str, str] = {}
LOCKS: dict[str, threading.Lock] = {}
GUARD = threading.Lock()


class Handler(SimpleHTTPRequestHandler):
    root: Path = Path(".")
    page: str = "index.html"
    cache: Path | None = None

    def send_json(self, obj: dict, code: int = 200) -> None:
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if u.path == "/api/progress":
            return self.send_json({"text": PROGRESS.get(q.get("key", [""])[0], "")})
        if u.path == "/api/explore":
            try:
                ra, dec = float(q["ra"][0]), float(q["dec"][0])
                if not (0 <= ra < 360 and -90 <= dec <= 90):
                    raise ValueError
            except (KeyError, ValueError):
                return self.send_json({"ok": False, "error": "Give ra (0-360) and dec (-90 to 90) in degrees."}, 400)
            key = f"{ra:.4f}_{dec:+.4f}".replace(".", "p").replace("+", "n").replace("-", "s")
            out = self.root / "explore" / "_live" / key
            with GUARD:
                lock = LOCKS.setdefault(key, threading.Lock())
            with lock:                                   # two people asking for the same spot share one fetch
                if not (out / "target.json").exists():
                    try:
                        explore(ra, dec, None, 0.1, 1, None, out, self.cache, 4,
                                progress=lambda t: PROGRESS.__setitem__(key, t.strip()))
                    except Exception as e:
                        return self.send_json({"ok": False, "error": str(e)}, 502)
            return self.send_json({"ok": True, "dir": f"explore/_live/{key}/", "key": key})
        if u.path == "/":
            self.path = "/" + self.page
        return super().do_GET()

    def end_headers(self):
        if not self.path.startswith("/api/"):
            self.send_header("Cache-Control", "no-cache")
        super().end_headers()


def main() -> int:
    ap = argparse.ArgumentParser(description="BlinkSky helper server")
    ap.add_argument("--port", type=int, default=8831)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--page", default="index.html")
    ap.add_argument("--cache", default=".tmp/blinksky_explorer_cache")
    a = ap.parse_args()
    Handler.root, Handler.page, Handler.cache = Path(a.root), a.page, Path(a.cache) if a.cache else None
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), partial(Handler, directory=a.root))
    print(f"BlinkSky helper on http://localhost:{a.port}/ (serving {a.root})", file=sys.stderr, flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
