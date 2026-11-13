"""
BlinkSky orbit check: where NASA JPL Horizons says every known object in a round should be, seen from SPHEREx.

What it does : Reads a round's round.json and, for each object, asks NASA JPL Horizons for its position at the
               time of each image, as seen from the SPHEREx spacecraft itself (Horizons ID -163182), plus its
               distance from the Sun and from SPHEREx. It works out how far the orbit says the object should
               move between the images, and compares that with how far it moved in the images.
               Writes horizons.json in the same folder.
Inputs       : --round DIR       folder holding round.json (e.g. ../real_round_full)
               --delay SECONDS   pause between requests to be polite to JPL (default 0.3)
Outputs      : DIR/horizons.json, keyed by the object name exactly as it appears in round.json:
               {"predicted_arcsec", "measured_arcsec" (null if not in both images), "sun_au", "spherex_au", ...}.
               JSON summary on stdout; progress on stderr.
Needs        : Python 3.10+, requests.

Example:
    python blinksky_horizons.py --round ../real_round_full

Data credit: NASA JPL Solar System Dynamics, Horizons API (https://ssd-api.jpl.nasa.gov/doc/horizons.html).
"""

import argparse
import json
import math
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
HEADERS = {"User-Agent": "BlinkSky-prototype (NASA Space Apps 2026)"}
SPHEREX = "500@-163182"   # Horizons: the SPHEREx spacecraft as the observer
AU_KM = 149_597_870.7


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def horizons_command(name: str) -> str:
    """'21 Lutetia (A852 VA)' -> '21;'  ·  '(2022 QJ66)' -> 'DES=2022 QJ66;'  (';' = small body only)."""
    m = re.match(r"\s*(\d+)\b", name)
    if m:
        return f"{m.group(1)};"
    m = re.search(r"\(([^)]+)\)", name)
    return f"DES={(m.group(1) if m else name).strip()};"


def to_jd(time_utc: str) -> float:
    """'2026-07-21 07:12 UTC' -> Julian Date."""
    t = datetime.strptime(time_utc.replace(" UTC", ""), "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
    return t.timestamp() / 86400 + 2440587.5


def separation_arcsec(ra1, dec1, ra2, dec2) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    h = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((r2 - r1) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(h))) * 3600


def ephemeris(name: str, jds: list[float]) -> list[dict]:
    """One row per time: ra, dec (deg), sun_au, spherex_au."""
    params = {"format": "json", "COMMAND": f"'{horizons_command(name)}'", "EPHEM_TYPE": "OBSERVER", "CENTER": f"'{SPHEREX}'",
              "TLIST": " ".join(f"'{jd:.6f}'" for jd in jds), "TLIST_TYPE": "JD", "QUANTITIES": "'1,19,20'",
              "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES", "OBJ_DATA": "NO"}
    r = requests.get(HORIZONS, params=params, headers=HEADERS, timeout=60)
    r.raise_for_status()
    text = r.json().get("result", "")
    if "$$SOE" not in text:
        raise ValueError(text.strip().splitlines()[-1][:120] if text.strip() else "empty answer")
    head = text[:text.index("$$SOE")].strip().splitlines()
    cols = [c.strip() for c in next(l for l in reversed(head) if "R.A." in l).split(",")]
    rows = []
    for line in text[text.index("$$SOE") + 5:text.index("$$EOE")].strip().splitlines():
        v = dict(zip(cols, (x.strip() for x in line.split(","))))
        ra = next(v[c] for c in cols if c.startswith("R.A."))
        dec = next(v[c] for c in cols if c.startswith("DEC"))
        rows.append({"ra": float(ra), "dec": float(dec), "sun_au": float(v["r"]), "spherex_au": float(v["delta"])})
    return rows


def load_frames(folder: Path, meta: dict):
    """Photo A and photo B as arrays, each with its difference image (stars cancel, movers stay).
    Only for full-resolution rounds (png_scale 1)."""
    if meta.get("png_scale") != 1:
        return None
    import numpy as np
    from PIL import Image
    a, b = (np.asarray(Image.open(folder / f["src"]).convert("L"), dtype=float) for f in meta["frames"][:2])
    return (a, a - b), (b, b - a)


def seen_at(frame, x: float, y: float, r: int = 4):
    """Centre of the spot near (x, y), measured in the photo itself, or None if the spot does not stand out
    in the difference image (too faint to see). The difference image only decides "visible?": its centre is
    pulled apart where the two spots overlap, so the position comes from the plain photo."""
    import numpy as np
    image, diff = frame
    xi, yi = int(round(x)), int(round(y))
    box = (slice(max(yi - 15, 0), yi + 16), slice(max(xi - 15, 0), xi + 16))
    win = (slice(max(yi - r, 0), yi + r + 1), slice(max(xi - r, 0), xi + r + 1))
    if diff[box].size == 0 or diff[win].size == 0:
        return None
    noise = 1.4826 * np.median(np.abs(diff[box] - np.median(diff[box]))) or 1.0
    if diff[win].max() < 6 * noise:
        return None
    w = (image[win] - np.median(image[box])).clip(min=0)
    ys, xs = np.mgrid[win]
    return float((xs * w).sum() / w.sum()), float((ys * w).sum() / w.sum())


def main() -> int:
    ap = argparse.ArgumentParser(description="Check every object in a BlinkSky round against its NASA JPL Horizons orbit")
    ap.add_argument("--round", required=True)
    ap.add_argument("--delay", type=float, default=0.3)
    a = ap.parse_args()
    folder = Path(a.round)
    try:
        meta = json.loads((folder / "round.json").read_text())
        jds = [to_jd(f["time_utc"]) for f in meta["frames"]]
        scale = meta.get("pixel_arcsec", 6.15) / (meta.get("png_scale") or 6)   # arcsec per PNG pixel (older rounds: 6x)
    except (OSError, KeyError, ValueError) as e:
        log(f"FAILED reading {folder / 'round.json'}: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1
    photos = load_frames(folder, meta)
    out, missing = {}, []
    for n, o in enumerate(meta["objects"], 1):
        try:
            rows = ephemeris(o["name"], jds)
            p = o["positions"]
            measured = round(math.hypot(p[1][0] - p[0][0], p[1][1] - p[0][1]) * scale, 1) if p[0] and p[1] else None
            out[o["name"]] = {"predicted_arcsec": round(separation_arcsec(rows[0]["ra"], rows[0]["dec"], rows[-1]["ra"], rows[-1]["dec"]), 1),
                              "measured_arcsec": measured, "sun_au": round(rows[0]["sun_au"], 3), "spherex_au": round(rows[0]["spherex_au"], 3),
                              "spherex_million_km": round(rows[0]["spherex_au"] * AU_KM / 1e6, 1),
                              "positions": [[round(r["ra"], 6), round(r["dec"], 6)] for r in rows]}
            e = out[o["name"]]
            # measured straight from the image pixels: where the spot really is in photo A and in photo B
            spots = [seen_at(photos[i], *p[i]) for i in (0, 1)] if photos is not None and p[0] and p[1] else [None, None]
            e["seen_arcsec"] = round(math.hypot(spots[1][0] - spots[0][0], spots[1][1] - spots[0][1]) * scale, 1) if all(spots) else None
            log(f"{n}/{len(meta['objects'])} {o['name']}: predicted {e['predicted_arcsec']}″, round file {measured}″, pixels {e['seen_arcsec']}″")
        except (requests.RequestException, ValueError, KeyError, StopIteration) as e:
            missing.append(o["name"])
            log(f"{n}/{len(meta['objects'])} {o['name']}: skipped ({e})")
        time.sleep(a.delay)
    (folder / "horizons.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    both = [v for v in out.values() if v["measured_arcsec"] is not None]
    diffs = sorted(abs(v["predicted_arcsec"] - v["measured_arcsec"]) for v in both)
    print(json.dumps({"ok": True, "written": str(folder / "horizons.json"), "objects": len(meta["objects"]), "with_orbit": len(out),
                      "in_both_images": len(both), "median_difference_arcsec": diffs[len(diffs) // 2] if diffs else None,
                      "seen_in_pixels": {k: [v["predicted_arcsec"], v["seen_arcsec"]] for k, v in out.items() if v.get("seen_arcsec") is not None},
                      "missing": missing}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
