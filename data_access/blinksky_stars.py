"""
BlinkSky stars: what are the things that do NOT move? Matches the visible stars in a full-image round to the Gaia catalogue.

What it does : Finds the visible point sources in frame_1.png of a full-image round, turns their pixel positions into
               RA/Dec with the round's saved TAN-SIP keywords, downloads Gaia DR3 stars for the same sky from
               NASA/IPAC IRSA (TAP), and matches each visible source to the brightest Gaia star within 1.5 SPHEREx pixels.
               Writes stars.json: for each matched source its PNG position, Gaia id, Gaia G brightness,
               distance in light years and temperature (when Gaia has them).
Inputs       : --round DIR        a full-image round folder (made by blinksky_full_round.py; needs "wcs" and "crop")
               --max-stars N      how many of the brightest visible sources to keep (default 5000)
               --gaia-gmax MAG    faintest Gaia star to download (default 16)
Outputs      : DIR/stars.json. JSON summary on stdout; progress on stderr.
Needs        : Python 3.10+, requests, numpy, scipy, pillow, and blinksky_full_round.py next to this file.

Example:
    python blinksky_stars.py --round ../real_round_full

Data credit: Gaia DR3 (ESA/Gaia/DPAC) served by NASA/IPAC IRSA; SPHEREx images via NASA/IPAC IRSA.
"""

import argparse
import csv
import io
import json
import math
import sys
from pathlib import Path

import numpy as np
import requests
from PIL import Image

from blinksky_full_round import TanSip

IRSA_TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
HEADERS = {"User-Agent": "BlinkSky-prototype (NASA Space Apps 2026)"}


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def detect(png: Path, keep: int):
    """Local brightness peaks in the (already stretched) image, brightest first."""
    from scipy.ndimage import maximum_filter
    img = np.asarray(Image.open(png).convert("L"), dtype=float)
    peaks = (img == maximum_filter(img, size=5)) & (img > 110)
    ys, xs = np.nonzero(peaks)
    order = np.argsort(-img[ys, xs])[:keep]
    return xs[order], ys[order], img[ys[order], xs[order]], img.shape


def gaia(ra: float, dec: float, radius: float, gmax: float) -> list[dict]:
    q = (f"SELECT source_id, ra, dec, phot_g_mean_mag, parallax, teff_gspphot, distance_gspphot FROM gaia_dr3_source "
         f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.5f},{dec:.5f},{radius:.4f}))=1 AND phot_g_mean_mag < {gmax}")
    r = requests.post(IRSA_TAP, data={"QUERY": q, "FORMAT": "csv"}, headers=HEADERS, timeout=900)
    r.raise_for_status()
    return list(csv.DictReader(io.StringIO(r.text)))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Match the visible stars of a BlinkSky round to Gaia")
    ap.add_argument("--round", required=True)
    ap.add_argument("--max-stars", type=int, default=5000)
    ap.add_argument("--gaia-gmax", type=float, default=16)
    a = ap.parse_args()
    folder = Path(a.round)
    try:
        meta = json.loads((folder / "round.json").read_text())
        wcs, crop = TanSip(meta["wcs"]), meta["crop"]
        xs, ys, vals, (ny, nx) = detect(folder / meta["frames"][0]["src"], a.max_stars)
        log(f"{len(xs)} visible sources found")
        # PNG pixel -> photo A pixel (FITS rows run bottom-up) -> RA/Dec
        ra, dec = wcs.pix2world(xs + crop["x0"], (ny - 1 - ys) + crop["y0"])
        ra_c, dec_c = meta["center"]["ra"], meta["center"]["dec"]
        radius = math.hypot(nx, ny) / 2 * meta["pixel_arcsec"] / 3600 * 1.02
        log(f"downloading Gaia DR3 stars brighter than G {a.gaia_gmax} within {radius:.2f} deg from NASA/IPAC IRSA ...")
        rows = gaia(ra_c, dec_c, radius, a.gaia_gmax)
        log(f"{len(rows)} Gaia stars")
    except (OSError, KeyError, requests.RequestException) as e:
        log(f"FAILED: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1
    from scipy.spatial import cKDTree
    cosd = math.cos(math.radians(dec_c))
    g_ra = np.array([float(r["ra"]) for r in rows]); g_dec = np.array([float(r["dec"]) for r in rows])
    g_mag = np.array([num(r["phot_g_mean_mag"]) or 99 for r in rows])
    tree = cKDTree(np.c_[(g_ra - ra_c) * cosd, g_dec - dec_c])
    match_deg = 1.5 * meta["pixel_arcsec"] / 3600
    stars = []
    for x, y, r_, d_ in zip(xs, ys, ra, dec):
        idx = tree.query_ball_point([(r_ - ra_c) * cosd, d_ - dec_c], match_deg)
        if not idx:
            continue
        best = rows[min(idx, key=lambda i: g_mag[i])]
        dist_pc = num(best["distance_gspphot"])
        if dist_pc is None and (num(best["parallax"]) or 0) > 0.2:
            dist_pc = 1000 / float(best["parallax"])
        teff = num(best["teff_gspphot"])
        stars.append([round(float(x) + .5, 1), round(float(y) + .5, 1), best["source_id"], round(float(best["phot_g_mean_mag"]), 1),
                      round(dist_pc * 3.2616) if dist_pc else None, round(teff) if teff else None])
    (folder / "stars.json").write_text(json.dumps({"source": "Gaia DR3 (ESA/Gaia/DPAC) via NASA/IPAC IRSA", "fields": ["x", "y", "gaia_dr3_id", "g_mag", "distance_ly", "teff_k"], "stars": stars}))
    print(json.dumps({"ok": True, "written": str(folder / "stars.json"), "visible_sources": int(len(xs)), "matched_to_gaia": len(stars),
                      "with_distance": sum(1 for s in stars if s[4]), "with_temperature": sum(1 for s in stars if s[5])}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
