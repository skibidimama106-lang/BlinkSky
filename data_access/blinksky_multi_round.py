"""
BlinkSky four-photo round: several SPHEREx photos of one patch of sky, taken at two moments.

What it does : SPHEREx photographs a patch a few times within minutes, then comes back hours later and does it again.
               This tool takes two photos from each of those two moments (four in all), so that
                 - a light seen in both photos of one moment is real, not a camera glitch, and
                 - a mover has four dated positions instead of two.
               It downloads the FULL images from NASA/IPAC IRSA (~65 MB each, cached in --cache), lines every photo
               up with the first one pixel by pixel, crops all of them to the area they share, writes one PNG per
               photo with the same brightness stretch, asks NASA JPL (sb_ident) which known objects are in each
               photo, and checks every object against its JPL Horizons orbit as seen from SPHEREx.
               round.json has the same layout as blinksky_full_round.py, with more frames; each frame carries
               "moment": 0 or 1.
Inputs       : --ra DEG --dec DEG     a point inside the wanted photos (default: where 21 Lutetia was on 21 Jul 2026)
               --detector D1..D6      SPHEREx detector / colour band (default D5)
               --gap-hours H          wanted time between the two moments (default 3.2)
               --per-moment N         photos to take from each moment (default 2)
               --obs A,B,C,D          pick the photos yourself by observation id (overrides the automatic choice)
               --vmag-lim MAG         only list known objects brighter than this (default 20)
               --cache DIR            where full FITS files are kept between runs
               --out DIR              round folder to write (e.g. ../real_round_4)
               --no-horizons          skip the JPL Horizons check (horizons.json is then not written)
Outputs      : DIR/round.json, DIR/frame_1.png ... frame_N.png, DIR/horizons.json. JSON summary on stdout.
Needs        : Python 3.10+, requests, numpy, scipy, pillow, astropy, and the other blinksky_*.py files next to this one.

Example:
    python blinksky_multi_round.py --cache ../../.tmp/streak_check/fits --out ../real_round_4

Data credit: SPHEREx/JPL/Caltech/NASA via NASA/IPAC IRSA; small-body data from NASA JPL Solar System Dynamics.
"""

import argparse
import json
import math
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).parent))
import blinksky_data as bd            # noqa: E402
import blinksky_full_round as fr      # noqa: E402
import blinksky_horizons as bh        # noqa: E402

MOMENT_GAP_DAYS = 30 / 1440           # photos less than 30 minutes apart belong to the same moment


def pick_photos(imgs: list[dict], gap_hours: float, per_moment: int) -> list[tuple[dict, int]]:
    """Group photos into moments, take the two moments closest to gap_hours apart, spread per_moment picks over each."""
    imgs = sorted(imgs, key=lambda i: i["time_mjd"])
    moments: list[list[dict]] = []
    for i in imgs:
        if moments and i["time_mjd"] - moments[-1][-1]["time_mjd"] < MOMENT_GAP_DAYS:
            moments[-1].append(i)
        else:
            moments.append([i])
    pairs = [(abs((b[0]["time_mjd"] - a[0]["time_mjd"]) * 24 - gap_hours), a, b) for n, a in enumerate(moments) for b in moments[n + 1:]]
    if not pairs:
        raise RuntimeError("SPHEREx photographed this spot at one moment only; pick another spot or detector")
    _, first, second = min(pairs, key=lambda t: t[0])
    out = []
    for m, group in enumerate((first, second)):
        n = min(per_moment, len(group))
        idx = sorted({round(k * (len(group) - 1) / max(n - 1, 1)) for k in range(n)})
        out += [(group[k], m) for k in idx]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Build a four-photo BlinkSky round")
    ap.add_argument("--ra", type=float, default=239.385)
    ap.add_argument("--dec", type=float, default=-20.72)
    ap.add_argument("--detector", default="D5")
    ap.add_argument("--gap-hours", type=float, default=3.2)
    ap.add_argument("--per-moment", type=int, default=2)
    ap.add_argument("--obs", default="")
    ap.add_argument("--vmag-lim", type=float, default=20)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-horizons", action="store_true")
    a = ap.parse_args()
    cache, out = Path(a.cache), Path(a.out)
    cache.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    try:
        imgs = [i for i in bd.find_spherex_images(a.ra, a.dec, 0.05) if i["detector"].endswith(a.detector)]
        if a.obs:
            want = [s.strip() for s in a.obs.split(",") if s.strip()]
            chosen = sorted((i for i in imgs if i["obs_id"] in want), key=lambda i: i["time_mjd"])
            if len(chosen) != len(want):
                raise RuntimeError(f"found {len(chosen)} of the {len(want)} photos asked for on {a.detector}")
            picks, m = [], 0
            for n, i in enumerate(chosen):
                if n and i["time_mjd"] - chosen[n - 1]["time_mjd"] >= MOMENT_GAP_DAYS:
                    m += 1
                picks.append((i, m))
        else:
            picks = pick_photos(imgs, a.gap_hours, a.per_moment)
        if len(picks) < 2:
            raise RuntimeError("fewer than two photos found")
        for img, m in picks:
            t = datetime(1858, 11, 17) + timedelta(days=img["time_mjd"])
            bd.log(f"photo: {img['obs_id']} {t:%Y-%m-%d %H:%M} (moment {m + 1})")

        data = []
        for img, m in picks:
            f = fr.download(img["url"], cache / f"{img['obs_id']}{img['detector'][-2:]}.fits")
            header, pixels = bd.read_fits(f.read_bytes())["IMAGE"]
            cleaned, n_bad = fr.clean_flagged(pixels, fr.read_flags(f))
            bd.log(f"    {img['obs_id']}: cleaned {n_bad:,} flagged pixels")
            data.append((img, m, header, fr.TanSip(header), cleaned))
        ref_img, _, ref_hdr, ref_wcs, ref_pix = data[0]
        layers, inside = [ref_pix], np.ones(ref_pix.shape, dtype=bool)
        for img, _, _, wcs, pix in data[1:]:
            bd.log(f"    lining {img['obs_id']} up with the first photo ...")
            lined, ok = fr.reproject_onto(pix, wcs, ref_wcs, ref_pix.shape)
            layers.append(lined); inside &= ok
        y0, y1, x0, x1 = fr.shared_rectangle(inside)
        bd.log(f"    area all photos share: x {x0}-{x1}, y {y0}-{y1} of {ref_pix.shape[1]}x{ref_pix.shape[0]}")
        crops = [layer[y0:y1, x0:x1] for layer in layers]
        ny, nx = crops[0].shape
        ra_c, dec_c = (float(v) for v in ref_wcs.pix2world(np.array((x0 + x1) / 2), np.array((y0 + y1) / 2)))
        half = max(nx, ny) * ref_img["pixel_arcsec"] / 3600 / 2 * 1.05
        lo, hi = (float(v) for v in np.nanpercentile(crops[0], [1, 99.7]))     # one stretch for every photo

        frames, names = [], {}
        for n, (img, m, *_rest) in enumerate(data):
            fname = f"frame_{n + 1}.png"
            fr.to_png_shared(crops[n], out / fname, lo, hi)
            when = bd.mjd_to_jpl_time(img["time_mjd"])
            bd.log(f"    asking NASA JPL for objects brighter than V{a.vmag_lim} at {when} ...")
            for o in fr.known_bright(ra_c, dec_c, when, round(half, 2), a.vmag_lim, cache):
                x, y = (float(v) for v in ref_wcs.world2pix(np.array(o["ra"]), np.array(o["dec"])))
                pos = [round(x - x0 + .5, 1), round(ny - 1 - (y - y0) + .5, 1)]    # PNG rows run top-down
                e = names.setdefault(o["name"], {"name": o["name"], "v_mag": o["v_mag"], "positions": [None] * len(data)})
                e["positions"][n] = pos if 0 <= pos[0] < nx and 0 <= pos[1] < ny else None
            t = datetime(1858, 11, 17) + timedelta(days=img["time_mjd"])
            frames.append({"src": fname, "time_utc": t.strftime("%Y-%m-%d %H:%M UTC"), "moment": m, "detector": img["detector"],
                           "wavelength_um": img["wavelength_um"], "obs_id": img["obs_id"], "width": nx, "height": ny})
        objects = sorted(names.values(), key=lambda o: (o["v_mag"] is None, o["v_mag"] or 99))
        round_ = {
            "source": "NASA SPHEREx Quick Release images via NASA/IPAC IRSA; known objects from NASA JPL sb_ident",
            "collection": "spherex_qr3", "center": {"ra": round(ra_c, 5), "dec": round(dec_c, 5)}, "full_image": True, "png_scale": 1, "aligned_to": "frame_1",
            "wcs": {k: ref_hdr[k] for k in ref_hdr if k.startswith(("CRPIX", "CRVAL", "CD", "PC", "CDELT", "A_", "B_", "AP_", "BP_")) and not k[-1] in "AW"},
            "crop": {"x0": x0, "y0": y0, "width": nx, "height": ny},
            "vmag_limit": a.vmag_lim, "pixel_arcsec": ref_img["pixel_arcsec"],
            "width": nx, "height": ny, "frames": frames, "objects": objects,
        }
        (out / "round.json").write_text(json.dumps(round_, indent=1))
    except (requests.RequestException, RuntimeError, KeyError) as e:
        bd.log(f"FAILED: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1

    seen = {}
    if not a.no_horizons:
        # The same check as blinksky_horizons.py, but from the first photo to the last one (not photo 1 to photo 2).
        from PIL import Image
        first, last = (np.asarray(Image.open(out / frames[k]["src"]).convert("L"), dtype=float) for k in (0, -1))
        photos = ((first, first - last), (last, last - first))
        jds, scale, hor = [bh.to_jd(f["time_utc"]) for f in frames], round_["pixel_arcsec"], {}
        for n, o in enumerate(objects, 1):
            try:
                rows, p = bh.ephemeris(o["name"], jds), o["positions"]
                both = p[0] and p[-1]
                e = hor[o["name"]] = {
                    "predicted_arcsec": round(bh.separation_arcsec(rows[0]["ra"], rows[0]["dec"], rows[-1]["ra"], rows[-1]["dec"]), 1),
                    "measured_arcsec": round(math.hypot(p[-1][0] - p[0][0], p[-1][1] - p[0][1]) * scale, 1) if both else None,
                    "sun_au": round(rows[0]["sun_au"], 3), "spherex_au": round(rows[0]["spherex_au"], 3),
                    "spherex_million_km": round(rows[0]["spherex_au"] * bh.AU_KM / 1e6, 1),
                    "positions": [[round(r["ra"], 6), round(r["dec"], 6)] for r in rows]}
                spots = [bh.seen_at(photos[0], *p[0]), bh.seen_at(photos[1], *p[-1])] if both else [None, None]
                e["seen_arcsec"] = round(math.hypot(spots[1][0] - spots[0][0], spots[1][1] - spots[0][1]) * scale, 1) if all(spots) else None
                if e["seen_arcsec"] is not None:
                    seen[o["name"]] = [e["predicted_arcsec"], e["seen_arcsec"]]
                bd.log(f"{n}/{len(objects)} {o['name']}: predicted {e['predicted_arcsec']}″, pixels {e['seen_arcsec']}″")
            except (requests.RequestException, ValueError, KeyError, StopIteration) as e:
                bd.log(f"{n}/{len(objects)} {o['name']}: skipped ({e})")
            time.sleep(0.3)
        (out / "horizons.json").write_text(json.dumps(hor, indent=1, ensure_ascii=False), encoding="utf-8")

    print(json.dumps({"ok": True, "out": str(out), "size_px": [nx, ny], "photos": [[f["obs_id"], f["time_utc"], f["moment"]] for f in frames],
                      "objects_listed": len(objects), "in_every_photo": sum(1 for o in objects if all(o["positions"])),
                      "brightest": [(o["name"], o["v_mag"]) for o in objects[:5]], "seen_in_pixels": seen}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
