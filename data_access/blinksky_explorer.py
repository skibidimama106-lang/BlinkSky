"""
BlinkSky Explorer data tool: one patch of sky -> a "little video" of SPHEREx photos + a colour fingerprint.

What it does : For a sky position (RA, Dec) it
               1. finds every SPHEREx image of that spot in NASA/IPAC IRSA (releases DR1 and QR3 together),
               2. groups them into visits (SPHEREx returns to a spot about every 6 months),
               3. downloads small cutouts, a few per visit per detector,
               4. turns each one into a north-up picture on the same grid, so the frames line up like a video,
               5. measures the brightest source near the centre in every frame: where it is, how bright,
                  which of SPHEREx's colours (wavelength) that photo was taken in, and any warning flags,
               6. fits a straight line to the positions (motion per year) and collects the brightness-by-colour
                  points (the "colour fingerprint", i.e. a spectrum built on the fly).
Inputs       : --ra DEG --dec DEG   centre of the patch, degrees (ICRS)
               --name TEXT          label for the target (default "RA .., Dec ..")
               --size DEG           patch width (default 0.1, about 58 SPHEREx pixels)
               --per-visit N        photos per detector per visit (default 2; 1 = quick)
               --detectors LIST     e.g. "2,5" (default all six)
               --out DIR            folder to write (default .tmp/blinksky_explore/<slug>)
               --cache DIR          keep downloaded cutouts here so a re-run is instant (default .tmp/blinksky_explorer_cache)
               --workers N          parallel downloads (default 4)
Outputs      : JSON summary on stdout. In --out: target.json (everything the web page needs) and one PNG per frame.
               Progress goes to stderr. Exit code 0 on success, 1 on failure.
Needs        : Python 3.10+, requests, numpy, scipy, astropy, pillow.

Example:
    python blinksky_explorer.py --ra 269.44627 --dec 4.76859 --name "Barnard's Star" --out ../explore/barnard

Data credit: SPHEREx/JPL/Caltech/NASA via NASA/IPAC IRSA. Releases DR1 (DOI 10.26131/IRSA678) and QR3 (DOI 10.26131/IRSA662); the five ready-made test objects were fetched from QR2 (DOI 10.26131/IRSA652) and QR3.
"""

import argparse
import hashlib
import json
import math
import sys
import time
import warnings
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from io import BytesIO
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from PIL import Image
from scipy.ndimage import map_coordinates

warnings.simplefilter("ignore", FITSFixedWarning)
warnings.simplefilter("ignore", RuntimeWarning)

IRSA_SIA = "https://irsa.ipac.caltech.edu/SIA"
VOT = "{http://www.ivoa.net/xml/VOTable/v1.3}"
HEADERS = {"Connection": "close", "User-Agent": "BlinkSky-explorer (NASA Space Apps 2026)"}
# DR1 (12 Nov 2026) reprocesses QR2 to match QR3; IRSA no longer returns QR2 in searches. A new release = one more line here.
RELEASES = [("spherex_dr1", "10.26131/IRSA678"), ("spherex_qr3", "10.26131/IRSA662")]
OUT_SCALE = 3.075          # arcsec per picture pixel (half a SPHEREx pixel, so re-gridding loses little)
VISIT_GAP_DAYS = 30
# Flag bits that mean "do not trust this pixel" (names from the FLAGS layer's own header).
BAD_BITS = {0: "TRANSIENT", 1: "OVERFLOW", 2: "SUR_ERROR", 6: "NONFUNC", 9: "MISSING_DATA", 10: "HOT",
            15: "NONLINEAR", 17: "PERSIST", 19: "OUTLIER", 26: "BLOOM", 27: "SNOWBALL", 28: "HALO"}


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def mjd_to_iso(mjd: float) -> str:
    return (datetime(1858, 11, 17) + timedelta(days=mjd)).strftime("%Y-%m-%d %H:%M UTC")


# ---------- 1. find the images ----------

def find_images(ra: float, dec: float, radius: float = 0.01) -> list[dict]:
    """Every SPHEREx science image covering the spot, from both public releases, oldest first."""
    out, seen = [], set()
    for collection, doi in RELEASES:
        r = requests.get(IRSA_SIA, params={"COLLECTION": collection, "POS": f"CIRCLE {ra} {dec} {radius}", "MAXREC": 3000},
                         headers=HEADERS, timeout=180)
        r.raise_for_status()
        root = ET.fromstring(r.content)
        fields = [f.get("name") for f in root.iter(VOT + "FIELD")]
        for tr in root.iter(VOT + "TR"):
            x = dict(zip(fields, [(td.text or "") for td in tr]))
            key = (x["obs_id"], x["energy_bandpassname"])     # the same exposure must not come twice from two releases
            if x.get("dataproduct_subtype") != "science" or key in seen:
                continue
            seen.add(key)
            out.append({"obs_id": x["obs_id"], "mjd": float(x["t_min"]), "detector": x["energy_bandpassname"],
                        "band_um": [round(float(x["em_min"]) * 1e6, 3), round(float(x["em_max"]) * 1e6, 3)],
                        "url": x["access_url"], "release": collection, "doi": doi})
    return sorted(out, key=lambda i: i["mjd"])


def split_visits(images: list[dict]) -> list[list[dict]]:
    visits: list[list[dict]] = []
    for im in images:
        if visits and im["mjd"] - visits[-1][-1]["mjd"] <= VISIT_GAP_DAYS:
            visits[-1].append(im)
        else:
            visits.append([im])
    return visits


def pick(visits: list[list[dict]], per_visit: int, detectors: set[str] | None) -> list[dict]:
    """A few photos per detector per visit, spread evenly through the visit."""
    chosen = []
    for v, visit in enumerate(visits):
        by_det: dict[str, list[dict]] = {}
        for im in visit:
            by_det.setdefault(im["detector"], []).append(im)
        for det, ims in sorted(by_det.items()):
            if detectors and det[-1] not in detectors:
                continue
            idx = sorted({round(k * (len(ims) - 1) / max(1, per_visit - 1)) for k in range(per_visit)}) if per_visit > 1 else [len(ims) // 2]
            for i in idx:
                chosen.append({**ims[i], "visit": v})
    return sorted(chosen, key=lambda i: i["mjd"])


# ---------- 2. download a cutout ----------

def get_cutout(im: dict, ra: float, dec: float, size: float, cache: Path | None) -> bytes:
    key = hashlib.sha1(f"{im['url']}|{ra:.6f}|{dec:.6f}|{size}".encode()).hexdigest()[:20]
    f = cache / f"{key}.fits" if cache else None
    if f and f.exists():
        return f.read_bytes()
    last = None
    for attempt in range(3):
        try:
            r = requests.get(im["url"], params={"center": f"{ra},{dec}", "size": str(size)}, headers=HEADERS, timeout=240)
            r.raise_for_status()
            if f:
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_bytes(r.content)
            return r.content
        except requests.RequestException as e:   # IRSA sometimes answers 503 for a while; try again
            last = e
            time.sleep(3 * (attempt + 1))
    raise last


# ---------- 3. measure one photo ----------

def wavelength_at(hdul, x_fits: float, y_fits: float) -> tuple[float, float]:
    """Each SPHEREx pixel sees its own colour. The file carries an 11x11 lookup table over the detector; interpolate it."""
    h = hdul["IMAGE"].header
    t = hdul["WCS-WAVE"].data[0]
    xs, ys, vals = np.asarray(t["X"], float), np.asarray(t["Y"], float), np.asarray(t["VALUES"], float)   # vals[y][x][wave, width]
    xr = x_fits - float(h["CRPIX1W"]) + float(h["CRVAL1W"])      # position on the full detector (raw FITS pixels)
    yr = y_fits - float(h["CRPIX2W"]) + float(h["CRVAL2W"])
    ix = int(np.clip(np.searchsorted(xs, xr) - 1, 0, len(xs) - 2))
    iy = int(np.clip(np.searchsorted(ys, yr) - 1, 0, len(ys) - 2))
    fx = (xr - xs[ix]) / (xs[ix + 1] - xs[ix])
    fy = (yr - ys[iy]) / (ys[iy + 1] - ys[iy])
    v = (vals[iy, ix] * (1 - fx) * (1 - fy) + vals[iy, ix + 1] * fx * (1 - fy)
         + vals[iy + 1, ix] * (1 - fx) * fy + vals[iy + 1, ix + 1] * fx * fy)
    return float(v[0]), float(v[1])


def measure(raw: bytes, ra: float, dec: float, size: float) -> dict | None:
    """Re-grid the photo north-up and measure the brightest source near the centre."""
    with fits.open(BytesIO(raw)) as hdul:
        img = hdul["IMAGE"].data.astype(float)
        hdr = hdul["IMAGE"].header
        flags = hdul["FLAGS"].data.astype(np.int64)
        var = hdul["VARIANCE"].data.astype(float)
        w = WCS(hdr)
        ny, nx = img.shape
        x0, y0 = (float(v) for v in w.all_world2pix(ra, dec, 0))
        if not (6 <= x0 < nx - 6 and 6 <= y0 < ny - 6):
            return None                                   # the spot sits on the edge of this photo
        bad_mask = np.zeros(img.shape, bool)
        for bit in BAD_BITS:
            bad_mask |= (flags >> bit) & 1 == 1
        sky = np.nanmedian(img[~bad_mask]) if (~bad_mask).any() else np.nanmedian(img)
        sig = 1.4826 * np.nanmedian(np.abs(img[~bad_mask] - sky)) if (~bad_mask).any() else 1.0

        # --- the picture: same north-up grid for every photo, so frames line up ---
        n = int(round(size * 3600 / OUT_SCALE))
        j, i = np.mgrid[0:n, 0:n]
        d_e = ((n - 1) / 2 - i) * OUT_SCALE / 3600        # east is left
        d_n = ((n - 1) / 2 - j) * OUT_SCALE / 3600        # north is up (row 0 = top)
        dec_g = dec + d_n
        ra_g = ra + d_e / np.cos(np.radians(dec_g))
        px, py = w.all_world2pix(ra_g, dec_g, 0)
        clean = np.where(np.isfinite(img), img, sky)
        # A star too bright for the detector leaves a hole in its own middle (OVERFLOW). Fill the hole with the
        # brightest good value so the picture looks right.
        # This is for the picture only. Measurements use the untouched pixels: filling moved the centre by up to 10".
        overflow = (flags >> 1) & 1 == 1
        disp = clean.copy()
        if overflow.any() and (~overflow).any():
            disp[overflow] = np.nanmax(clean[~overflow])
        pic = map_coordinates(disp, [py, px], order=1, mode="constant", cval=np.nan)
        pic_bad = map_coordinates(bad_mask.astype(float), [py, px], order=0, mode="constant", cval=0) > 0

        # --- the measurement: brightest source within ~18 arcsec of the centre ---
        cx, cy = int(round(x0)), int(round(y0))
        best, bx, by = -np.inf, cx, cy
        for yy in range(cy - 3, cy + 4):
            for xx in range(cx - 3, cx + 4):
                if math.hypot(xx - x0, yy - y0) > 3.2:
                    continue
                s = np.nansum(clean[yy - 1:yy + 2, xx - 1:xx + 2])
                if s > best:
                    best, bx, by = s, xx, yy
        yy, xx = np.mgrid[0:ny, 0:nx]
        r = np.hypot(xx - bx, yy - by)
        ring = (r >= 6) & (r <= 10) & ~bad_mask & np.isfinite(img)
        bg = np.nanmedian(img[ring]) if ring.sum() > 20 else sky
        box = (slice(max(0, by - 2), by + 3), slice(max(0, bx - 2), bx + 3))
        cut = clean[box] - bg
        flux = float(np.nansum(cut))
        err = float(np.sqrt(np.nansum(np.where(np.isfinite(var[box]), var[box], 0)))) or float(sig * 5)
        pos = np.clip(cut, 0, None)
        if pos.sum() > 0:
            gy, gx = np.mgrid[box]
            mx, my = float((pos * gx).sum() / pos.sum()), float((pos * gy).sum() / pos.sum())
        else:
            mx, my = float(bx), float(by)
        s_ra, s_dec = (float(v) for v in w.all_pix2world(mx, my, 0))
        off_e = (s_ra - ra) * math.cos(math.radians(dec)) * 3600
        off_n = (s_dec - dec) * 3600
        core = flags[max(0, by - 1):by + 2, max(0, bx - 1):bx + 2]
        warn = sorted({name for bit, name in BAD_BITS.items() if ((core >> bit) & 1).any()})
        wave, width = wavelength_at(hdul, mx + 1, my + 1)
        pix_sr = (6.15 / 206265) ** 2             # one SPHEREx pixel on the sky; MJy/sr summed x this = MJy
        return {"pic": pic, "pic_bad": pic_bad, "sky": float(sky), "sigma": float(sig),
                "flux_mjy": round(flux * pix_sr * 1e9, 4), "err_mjy": round(err * pix_sr * 1e9, 4),
                "snr": round(flux / err, 1) if err > 0 else 0.0,
                "off_e": round(off_e, 2), "off_n": round(off_n, 2),
                "wave_um": round(wave, 4), "width_um": round(width, 4), "flags": warn,
                "date_obs": hdr.get("DATE-AVG") or hdr.get("DATE-OBS")}


def save_png(pic: np.ndarray, lo: float, hi: float, path: Path) -> None:
    v = np.where(np.isfinite(pic), pic, lo)
    s = np.arcsinh(np.clip((v - lo) / max(hi - lo, 1e-9), 0, 1) * 60) / np.arcsinh(60)
    Image.fromarray((s * 255).astype(np.uint8)).save(path, optimize=True)


# ---------- 4. put it together ----------

def fit_motion(frames: list[dict]) -> dict | None:
    """Straight-line fit of the measured positions against time, using every clearly detected frame.
    One pass throws out frames that sit far from the line (a glitch, or a neighbour star picked by mistake)."""
    good = [f for f in frames if f["snr"] >= 8]
    if len(good) < 4:
        return None
    if np.ptp([f["mjd"] for f in good]) < 70:
        return None
    for _ in range(2):
        t = np.array([f["mjd"] for f in good]); t_yr = (t - t.mean()) / 365.25
        a = np.vstack([t_yr, np.ones_like(t_yr)]).T
        dist = np.zeros(len(good))
        for key in ("off_e", "off_n"):
            y = np.array([f[key] for f in good])
            coef, *_ = np.linalg.lstsq(a, y, rcond=None)
            dist += (y - a @ coef) ** 2
        dist = np.sqrt(dist)
        keep = dist <= max(3 * np.median(dist) * 1.4826, 1.0)
        if keep.all() or keep.sum() < 4:
            break
        good = [f for f, k in zip(good, keep) if k]
    t = np.array([f["mjd"] for f in good]); t_yr = (t - t.mean()) / 365.25
    res = {}
    for key in ("off_e", "off_n"):
        y = np.array([f[key] for f in good])
        a = np.vstack([t_yr, np.ones_like(t_yr)]).T
        coef, *_ = np.linalg.lstsq(a, y, rcond=None)
        resid = y - a @ coef
        scatter = float(np.sqrt((resid ** 2).sum() / max(1, len(y) - 2)))
        res[key] = (float(coef[0]), scatter / math.sqrt(((t_yr - t_yr.mean()) ** 2).sum()), scatter)
    e, n = res["off_e"], res["off_n"]
    total = math.hypot(e[0], n[0])
    err = math.hypot(e[0] * e[1], n[0] * n[1]) / total if total > 0 else math.hypot(e[1], n[1])
    return {"east_arcsec_per_yr": round(e[0], 2), "north_arcsec_per_yr": round(n[0], 2),
            "total_arcsec_per_yr": round(total, 2), "error_arcsec_per_yr": round(err, 2),
            "scatter_arcsec": round(math.hypot(e[2], n[2]), 2), "frames_used": len(good), "frames_with_warning_flags": sum(1 for f in good if f["flags"]),
            "span_years": round(float(np.ptp(t)) / 365.25, 2)}


def explore(ra: float, dec: float, name: str | None, size: float, per_visit: int, detectors: set[str] | None,
            out: Path, cache: Path | None, workers: int = 4, progress=log) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    progress(f"1/3 Searching NASA/IPAC IRSA for SPHEREx photos of RA {ra:.5f}, Dec {dec:+.5f} ...")
    images = find_images(ra, dec)
    visits = split_visits(images)
    chosen = pick(visits, per_visit, detectors)
    progress(f"    {len(images)} photos in {len(visits)} visits; fetching {len(chosen)} cutouts")
    if not chosen:
        raise RuntimeError("SPHEREx has no public photos of this spot yet")

    done = [0]
    def work(im):
        try:
            raw = get_cutout(im, ra, dec, size, cache)
            m = measure(raw, ra, dec, size)
        except Exception as e:                     # one bad photo must not sink the whole request
            m = None; log(f"    skipped {im['obs_id']}: {e}")
        done[0] += 1
        progress(f"2/3 Photo {done[0]} of {len(chosen)}")
        return im, m
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = [r for r in pool.map(work, chosen) if r[1] is not None]
    if not results:
        raise RuntimeError("no usable SPHEREx photo of this spot")

    progress("3/3 Lining up the frames and measuring ...")
    # One shared brightness scale per detector, so a real change in brightness is visible between frames.
    hi_by_det: dict[str, float] = {}
    for det in {im["detector"] for im, _ in results}:
        tops = [np.nanpercentile(m["pic"] - m["sky"], 99.8) for im, m in results if im["detector"] == det]
        hi_by_det[det] = float(np.median(tops))
    frames = []
    for k, (im, m) in enumerate(results):
        src = f"f{k:03d}.png"
        save_png(m["pic"] - m["sky"], -1.5 * m["sigma"], max(hi_by_det[im["detector"]], 8 * m["sigma"]), out / src)
        frames.append({"src": src, "mjd": round(im["mjd"], 5), "time_utc": mjd_to_iso(im["mjd"]), "visit": im["visit"],
                       "detector": im["detector"][-2:], "obs_id": im["obs_id"], "release": im["release"], "file": im["url"],
                       **{k2: m[k2] for k2 in ("flux_mjy", "err_mjy", "snr", "off_e", "off_n", "wave_um", "width_um", "flags")}})
    by_det: dict[str, list[float]] = {}
    for f in frames:
        by_det.setdefault(f["detector"], []).append(f["snr"])
    best_det = max(by_det, key=lambda d: (len(by_det[d]) >= 3, float(np.median(by_det[d]))))
    n_side = int(round(size * 3600 / OUT_SCALE))
    target = {
        "name": name or f"RA {ra:.4f}, Dec {dec:+.4f}", "ra": ra, "dec": dec, "size_deg": size,
        "picture_px": n_side, "arcsec_per_px": OUT_SCALE,
        "source": "NASA SPHEREx spectral images via NASA/IPAC IRSA",
        "releases": [{"collection": c, "doi": d, "photos": sum(1 for i in images if i["release"] == c)} for c, d in RELEASES],
        "photos_available": len(images),
        "visits": [{"from": mjd_to_iso(v[0]["mjd"])[:10], "to": mjd_to_iso(v[-1]["mjd"])[:10], "photos": len(v)} for v in visits],
        "default_detector": best_det,
        "motion": fit_motion(frames),
        "colours_collected": len({(f["detector"], round(f["wave_um"], 2)) for f in frames if f["snr"] >= 3}),
        "frames": frames,
    }
    (out / "target.json").write_text(json.dumps(target, indent=1))
    return target


def main() -> int:
    ap = argparse.ArgumentParser(description="BlinkSky Explorer data tool")
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--name")
    ap.add_argument("--size", type=float, default=0.1)
    ap.add_argument("--per-visit", type=int, default=2)
    ap.add_argument("--detectors")
    ap.add_argument("--out")
    ap.add_argument("--cache", default=".tmp/blinksky_explorer_cache")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    slug = "".join(ch if ch.isalnum() else "_" for ch in (a.name or f"{a.ra:.4f}_{a.dec:+.4f}")).strip("_").lower()
    out = Path(a.out) if a.out else Path(".tmp/blinksky_explore") / slug
    dets = set(a.detectors.replace("D", "").split(",")) if a.detectors else None
    try:
        t = explore(a.ra, a.dec, a.name, a.size, a.per_visit, dets, out, Path(a.cache) if a.cache else None, a.workers)
    except Exception as e:
        log(f"FAILED: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1
    print(json.dumps({"ok": True, "out": str(out), "name": t["name"], "photos_available": t["photos_available"],
                      "visits": t["visits"], "frames": len(t["frames"]), "default_detector": t["default_detector"],
                      "motion": t["motion"], "colours_collected": t["colours_collected"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
