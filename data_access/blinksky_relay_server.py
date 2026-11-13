"""
BlinkSky relay server: the helper server plus the steps between a player's catch and "the gate".

What it does : Everything blinksky_server.py does (serves the game, fetches SPHEREx photos for typed coordinates),
               and adds the review chain for catches that are on no list:
                 POST /api/catch             a player sends a catch (JSON). The server checks it, drops repeats,
                                             limits floods, and puts it with other players' catches of the same change.
                 GET  /api/catch?id=&pid=    where that catch stands: {"status", "agree", "needed", "gate", "answer"}
                 GET  /review                a plain page for the team: every change, who caught it, Approve / Reject
                 POST /api/review            the Approve / Reject buttons on that page
               Status goes: sent -> review (once enough independent players caught the same change) -> approved / rejected.
               A flag only counts as independent when it comes from a different device, from a different place
               (a different network: the first three numbers of the internet address; only a scrambled code of it is
               kept, never the address) and at a different time (--gap-minutes away from every other counted flag).
               Nothing is sent to any outside registry. Approval only marks the draft report as ready.
               The review page and its buttons answer only on this computer, even when --host opens the game to others.
Inputs       : --port N        (default 8833)
               --root DIR      folder to serve (default: the folder above this file)
               --page NAME     file to open at "/" (default: index.html)
               --cache DIR     cutout cache (default .tmp/blinksky_explorer_cache)
               --store DIR     where catches are kept (default .tmp/blinksky_relay)
               --agree N       how many independent players must catch the same change before review (default 3)
               --gap-minutes M counted flags must be at least this many minutes apart (default 5; 0 switches it off)
               --same-place-ok count flags from the same network too (for a demo with phones on one Wi-Fi)
               --host ADDR     address to listen on (default 127.0.0.1; 0.0.0.0 lets phones on the same Wi-Fi play)
               --per-hour N    most catches one player may send in an hour (default 20)
Outputs      : Runs until stopped. Catches are kept in <store>/catches.json. Requests are logged to stderr.
Needs        : Python 3.10+ and everything blinksky_server.py needs.

Example:
    python blinksky_relay_server.py --port 8833 --page index.v6.claude.html
"""

import argparse
import hashlib
import html
import ipaddress
import json
import math
import sys
import threading
import time
from functools import partial
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

sys.path.insert(0, str(Path(__file__).parent))
from blinksky_server import Handler as BaseHandler   # noqa: E402

GATES = {"moved": "Minor Planet Center", "appeared": "Transient Name Server", "vanished": None}
SAME_SPOT_ARCSEC = 15.0          # about 2.5 SPHEREx pixels: two catches this close are the same change
STORE_LOCK = threading.Lock()


def sep_arcsec(a, b) -> float:
    """Angle between two sky positions [ra, dec] in degrees, in arcseconds (small-angle, fine below a degree)."""
    dra = (a[0] - b[0]) * math.cos(math.radians((a[1] + b[1]) / 2))
    return math.hypot(dra, a[1] - b[1]) * 3600


def place_code(ip: str) -> str:
    """A short scrambled code for the network a flag came from. Same Wi-Fi or same provider block = same code."""
    try:
        addr = ipaddress.ip_address(ip)
        net = ipaddress.ip_network(f"{ip}/{24 if addr.version == 4 else 64}", strict=False)
        key = str(net.network_address)
    except ValueError:
        key = "unknown"
    return hashlib.sha256(("blinksky-place:" + key).encode()).hexdigest()[:8]


def clean(text, limit: int) -> str:
    return "".join(ch for ch in str(text or "") if ch.isprintable())[:limit]


class Relay:
    def __init__(self, store: Path, needed: int, per_hour: int, gap_minutes: float = 5, same_place_ok: bool = False):
        self.file, self.needed, self.per_hour = store / "catches.json", needed, per_hour
        self.gap, self.same_place_ok = gap_minutes * 60, same_place_ok
        store.mkdir(parents=True, exist_ok=True)
        self.groups = json.loads(self.file.read_text("utf-8")) if self.file.exists() else []

    def save(self) -> None:
        self.file.write_text(json.dumps(self.groups, indent=1), "utf-8")

    def independent(self, g: dict) -> list[str]:
        """The flags that count: earliest first, each from a new place and far enough in time from every counted one."""
        kept: list[str] = []
        for pid, m in sorted(g["members"].items(), key=lambda kv: kv[1]["ts"]):
            others = [g["members"][k] for k in kept]
            if not self.same_place_ok and any(o.get("place") == m.get("place") for o in others):
                continue
            if any(abs(m["ts"] - o["ts"]) < self.gap for o in others):
                continue
            kept.append(pid)
        return kept

    def state(self, g: dict) -> dict:
        agree = len(self.independent(g))
        status = g["status"] if g["status"] != "open" else ("review" if agree >= self.needed else "sent")
        return {"status": status, "agree": agree, "flagged": len(g["members"]), "needed": self.needed,
                "gate": GATES.get(g["kind"]), "answer": g.get("answer", "")}

    def add(self, c: dict, ip: str = "", now: float | None = None) -> tuple[dict, int]:
        """Second check of a catch, then file it. Returns (reply, http code)."""
        try:
            pid, cid, kind = clean(c.get("pid"), 16), clean(c.get("id"), 12), c.get("kind")
            spots = [[float(c["ra"]), float(c["dec"])]]
            if c.get("ra2") is not None and c.get("dec2") is not None:
                spots.append([float(c["ra2"]), float(c["dec2"])])
            if not pid or not cid or kind not in GATES or not c.get("round"):
                raise ValueError
            if any(not (0 <= ra < 360 and -90 <= dec <= 90) for ra, dec in spots):
                raise ValueError
        except (KeyError, TypeError, ValueError):
            return {"ok": False, "error": "This catch needs a player, a receipt, a kind of change and a sky position."}, 400
        now, rnd = now if now is not None else time.time(), clean(c["round"], 60)
        with STORE_LOCK:
            for g in self.groups:                                   # the same player sending the same change again
                if g["round"] == rnd and pid in g["members"] and any(sep_arcsec(s, p) < SAME_SPOT_ARCSEC for s in spots for p in g["spots"]):
                    return {"ok": True, "repeat": True, "state": self.state(g)}, 200
            recent = sum(1 for g in self.groups for m, v in g["members"].items() if m == pid and now - v["ts"] < 3600)
            if recent >= self.per_hour:
                return {"ok": False, "error": "Too many catches sent this hour. Try again later."}, 429
            member = {"id": cid, "tag": clean(c.get("tag"), 12) or "ANON", "ts": now, "place": place_code(ip)}
            for g in self.groups:                                   # another player already caught this change
                if g["round"] == rnd and g["status"] == "open" and any(sep_arcsec(s, p) < SAME_SPOT_ARCSEC for s in spots for p in g["spots"]):
                    g["members"][pid] = member
                    self.save()
                    return {"ok": True, "state": self.state(g)}, 200
            g = {"gid": f"G{len(self.groups) + 1:04d}", "round": rnd, "kind": kind, "spots": spots, "status": "open",
                 "time": clean(c.get("time"), 40), "time2": clean(c.get("time2"), 40), "band": clean(c.get("band"), 60),
                 "what": clean(c.get("what"), 160), "members": {pid: member}}
            self.groups.append(g)
            self.save()
            return {"ok": True, "state": self.state(g)}, 200

    def find(self, cid: str, pid: str):
        for g in self.groups:
            m = g["members"].get(pid)
            if m and m["id"] == cid:
                return g
        return None

    def decide(self, gid: str, decision: str) -> str:
        with STORE_LOCK:
            g = next((x for x in self.groups if x["gid"] == gid), None)
            if not g or g["status"] != "open":
                return "That change is already decided."
            day = time.strftime("%d %b %Y")
            if decision == "approve":
                if len(self.independent(g)) < self.needed:
                    return f"Needs {self.needed} independent players (different time and place) before it can be approved."
                gate = GATES.get(g["kind"])
                g["status"] = "approved"
                g["answer"] = (f"Passed our review on {day}. Draft report ready for the {gate}. Prototype: not sent." if gate
                               else f"Passed our review on {day}. No registry takes this kind of change yet.")
            else:
                g["status"] = "rejected"
                g["answer"] = f"Reviewed on {day}: not passed on. Most likely a camera glitch."
            self.save()
            return ""

    def page(self, note: str = "") -> str:
        rows = []
        for g in reversed(self.groups):
            st, e = self.state(g), html.escape
            counted = set(self.independent(g))
            who = "<br>".join(f'{"✓" if pid in counted else "–"} {e(m["tag"])} · {time.strftime("%d %b %H:%M", time.localtime(m["ts"]))} · place {e(m.get("place", "?"))}'
                              for pid, m in sorted(g["members"].items(), key=lambda kv: kv[1]["ts"]))
            spot = " and ".join(f"{ra:.5f}, {dec:+.5f}" for ra, dec in g["spots"])
            if g["status"] == "open":
                can = st["agree"] >= self.needed
                act = (f'<form method="post" action="/api/review"><input type="hidden" name="gid" value="{e(g["gid"])}">'
                       f'<button name="decision" value="approve"{"" if can else " disabled"}>Approve</button> '
                       f'<button name="decision" value="reject">Reject</button></form>')
            else:
                act = e(g.get("answer", ""))
            rows.append(f'<tr><td>{e(g["gid"])}</td><td>{e(g["kind"])}<br><small>{e(st["gate"] or "no registry yet")}</small></td>'
                        f'<td>A light {e(g["what"])}<br><small>{e(g["time"])} · {e(g["band"])}<br>RA, Dec: {spot}</small></td>'
                        f'<td>{st["agree"]} of {st["needed"]} count<br><small>{st["flagged"]} flagged<br>{who}</small></td><td>{e(st["status"])}</td><td>{act}</td></tr>')
        body = "".join(rows) or '<tr><td colspan="6">No catches sent yet.</td></tr>'
        return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BlinkSky review</title><style>
body{{background:#0b0b0b;color:#eee;font:15px/1.45 system-ui,sans-serif;margin:0;padding:24px 16px;max-width:1000px;margin:auto}}
h1{{font-weight:500;font-size:22px}}p{{color:#aaa}}.note{{color:#f0c860}}
.wrap{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;min-width:720px}}
td,th{{border-top:1px solid #2c2c2c;padding:10px 8px;text-align:left;vertical-align:top}}th{{color:#999;font-weight:500;font-size:12px;text-transform:uppercase;letter-spacing:.08em}}
small{{color:#999}}button{{background:#eee;color:#000;border:0;border-radius:999px;padding:7px 14px;font:inherit;cursor:pointer}}
button[value=reject]{{background:none;color:#eee;border:1px solid #555}}button:disabled{{opacity:.35;cursor:default}}
</style></head><body><h1>BlinkSky review</h1>
<p>Changes players caught that are on no list. A change can be approved once {self.needed} independent players flagged it:
different devices{'' if self.same_place_ok else ', different places (networks)'}{f', and at least {self.gap / 60:g} minutes apart' if self.gap else ''}.
A tick marks the flags that count. Approving marks the draft report as ready. Nothing is sent to any registry.</p>
{f'<p class="note">{html.escape(note)}</p>' if note else ''}
<div class="wrap"><table><tr><th>No.</th><th>Kind / registry</th><th>What and where</th><th>Players</th><th>Status</th><th>Decision</th></tr>{body}</table></div>
</body></html>"""


class Handler(BaseHandler):
    relay: Relay = None  # type: ignore[assignment]

    def send_html(self, text: str, code: int = 200) -> None:
        body = text.encode()
        self.send_response(code)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def local(self) -> bool:
        try:
            return ipaddress.ip_address(self.client_address[0]).is_loopback
        except ValueError:
            return False

    def do_GET(self):  # noqa: N802
        u = urlparse(self.path)
        q = parse_qs(u.query)
        if u.path == "/api/catch":
            g = self.relay.find(q.get("id", [""])[0], q.get("pid", [""])[0])
            return self.send_json({"ok": True, "state": self.relay.state(g)} if g else {"ok": False, "error": "Not found."}, 200 if g else 404)
        if u.path == "/review":
            if not self.local():
                return self.send_html("<p>The review page opens only on the computer that runs BlinkSky.</p>", 403)
            return self.send_html(self.relay.page(q.get("note", [""])[0]))
        return super().do_GET()

    def do_POST(self):  # noqa: N802
        u = urlparse(self.path)
        size = int(self.headers.get("Content-Length") or 0)
        if size > 20_000:
            return self.send_json({"ok": False, "error": "Too large."}, 413)
        raw = self.rfile.read(size)
        if u.path == "/api/catch":
            try:
                data = json.loads(raw or b"{}")
                if not isinstance(data, dict):
                    raise ValueError
            except ValueError:
                return self.send_json({"ok": False, "error": "Not readable."}, 400)
            reply, code = self.relay.add(data, self.client_address[0])
            return self.send_json(reply, code)
        if u.path == "/api/review":
            if not self.local():
                return self.send_json({"ok": False, "error": "Review works only on the computer that runs BlinkSky."}, 403)
            f = parse_qs(raw.decode("utf-8", "replace"))
            note = self.relay.decide(f.get("gid", [""])[0], f.get("decision", [""])[0])
            self.send_response(303)
            self.send_header("Location", "/review" + (f"?note={note.replace(' ', '+')}" if note else ""))
            self.send_header("Content-Length", "0")
            self.end_headers()
            return None
        return self.send_json({"ok": False, "error": "Unknown address."}, 404)


def main() -> int:
    ap = argparse.ArgumentParser(description="BlinkSky relay server (helper server + review chain)")
    ap.add_argument("--port", type=int, default=8833)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--page", default="index.html")
    ap.add_argument("--cache", default=".tmp/blinksky_explorer_cache")
    ap.add_argument("--store", default=".tmp/blinksky_relay")
    ap.add_argument("--agree", type=int, default=3)
    ap.add_argument("--per-hour", type=int, default=20)
    ap.add_argument("--gap-minutes", type=float, default=5)
    ap.add_argument("--same-place-ok", action="store_true")
    ap.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args()
    if a.agree < 1 or a.per_hour < 1 or a.gap_minutes < 0:
        print("--agree and --per-hour must be 1 or more, and --gap-minutes 0 or more.", file=sys.stderr)
        return 2
    Handler.root, Handler.page, Handler.cache = Path(a.root), a.page, Path(a.cache) if a.cache else None
    Handler.relay = Relay(Path(a.store), a.agree, a.per_hour, a.gap_minutes, a.same_place_ok)
    srv = ThreadingHTTPServer((a.host, a.port), partial(Handler, directory=a.root))
    print(f"BlinkSky relay on http://localhost:{a.port}/  (review page: /review, catches in {a.store})", file=sys.stderr, flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
