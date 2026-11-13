"""
BlinkSky Gaia overlay: the stars ESA's Gaia mission lists around an Explorer target, ready to draw on the SPHEREx frames.

What it does : For one Explorer folder (made by blinksky_explorer.py) it asks NASA/IPAC IRSA for the Gaia DR3 stars
               in that patch, moves each one from Gaia's date (2016.0) to the wanted year with Gaia's own motion,
               and writes where each star sits relative to the patch centre (arcseconds east and north).
               The app draws these as circles on top of the SPHEREx photos: a light with no circle is not in Gaia.
Inputs       : --folder DIR       an Explorer folder holding target.json (e.g. ../explore/barnard); repeatable
               --gmax MAG         faintest Gaia star to keep (default 18)
               --year Y           year to move the stars to (default 2026.0)
Outputs      : DIR/gaia.json for every folder. JSON summary on stdout.
Needs        : Python 3.10+, requests. Free service, no key.

Example:
    python blinksky_gaia_overlay.py --folder ../explore/barnard --folder ../explore/rwaur

Data credit: Gaia DR3 (ESA/Gaia/DPAC), served by NASA/IPAC IRSA.
"""

import argparse
import csv
import io
import json
import math
import sys
from pathlib import Path

import requests

IRSA_TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
HEADERS = {"User-Agent": "BlinkSky-hackathon-tool"}


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def gaia(ra: float, dec: float, radius: float, gmax: float) -> list[dict]:
    q = (f"SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaia_dr3_source "
         f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.5f},{dec:.5f},{radius:.4f}))=1 AND phot_g_mean_mag < {gmax}")
    r = requests.post(IRSA_TAP, data={"QUERY": q, "FORMAT": "csv"}, headers=HEADERS, timeout=300)
    r.raise_for_status()
    return list(csv.DictReader(io.StringIO(r.text)))


def num(v, default=0.0):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def main() -> int:
    ap = argparse.ArgumentParser(description="Write gaia.json (Gaia DR3 stars as offsets) into BlinkSky Explorer folders")
    ap.add_argument("--folder", action="append", required=True)
    ap.add_argument("--gmax", type=float, default=18)
    ap.add_argument("--year", type=float, default=2026.0)
    a = ap.parse_args()
    done = []
    for folder in map(Path, a.folder):
        try:
            t = json.loads((folder / "target.json").read_text(encoding="utf-8"))
            ra0, dec0 = t["ra"], t["dec"]
            # a little wider than the picture, so fast stars that drifted in since 2016 are not missed
            rows = gaia(ra0, dec0, t["size_deg"] * 0.75 + 0.03, a.gmax)
            years, stars = a.year - 2016.0, []
            half = t["picture_px"] * t["arcsec_per_px"] / 2
            for r in rows:
                e = (num(r["ra"]) - ra0) * math.cos(math.radians(dec0)) * 3600 + num(r["pmra"]) / 1000 * years
                n = (num(r["dec"]) - dec0) * 3600 + num(r["pmdec"]) / 1000 * years
                if abs(e) <= half and abs(n) <= half:
                    moved = math.hypot(num(r["pmra"]), num(r["pmdec"])) / 1000
                    stars.append([round(e, 2), round(n, 2), round(num(r["phot_g_mean_mag"], 99), 1), round(moved, 3), r["source_id"]])
            stars.sort(key=lambda s: s[2])
            (folder / "gaia.json").write_text(json.dumps({
                "source": "Gaia DR3 (ESA/Gaia/DPAC) via NASA/IPAC IRSA", "year": a.year, "g_mag_limit": a.gmax,
                "fields": ["east_arcsec", "north_arcsec", "g_mag", "motion_arcsec_per_year", "gaia_dr3_id"], "stars": stars}), encoding="utf-8")
            log(f"{folder.name}: {len(stars)} Gaia stars in the picture")
            done.append({"folder": str(folder), "stars": len(stars)})
        except (requests.RequestException, OSError, KeyError, ValueError) as e:
            log(f"FAILED {folder}: {e}")
            print(json.dumps({"ok": False, "folder": str(folder), "error": str(e)}))
            return 1
    print(json.dumps({"ok": True, "written": done}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
