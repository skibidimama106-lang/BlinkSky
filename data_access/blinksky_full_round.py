"""
BlinkSky full-image round: a whole SPHEREx image pair as one big round, like an IASC image set.

What it does : Finds the SPHEREx blink pair (same detector, about --gap-hours apart) covering --ra/--dec,
               downloads the two FULL images from NASA/IPAC IRSA (~65 MB each, cached in --cache),
               lines photo B up with photo A pixel by pixel (full TAN-SIP positions), crops both to the area they share,
               writes each as a PNG at one PNG pixel per SPHEREx pixel (up to 2040 x 2040),
               and asks NASA JPL (sb_ident) which known objects brighter than --vmag-lim are in it.
               Output uses the same round.json format as blinksky_data.py --export-round, plus "png_scale": 1.
Inputs       : --ra DEG --dec DEG     a point inside the wanted image (default: BlinkSky's ecliptic test patch)
               --gap-hours H          wanted time between the two photos (default 3.17)
               --detector D1..D6      SPHEREx detector / colour band (default D5)
               --vmag-lim MAG         only list known objects brighter than this (default 22; fainter ones are
                                      far below what one SPHEREx exposure can show)
               --cache DIR            where full FITS files are kept between runs
               --out DIR              round folder to write (e.g. ../real_round_full)
Outputs      : DIR/round.json, DIR/frame_1.png, DIR/frame_2.png. JSON summary on stdout; progress on stderr.
Needs        : Python 3.10+, requests, numpy, scipy, pillow, astropy (to unpack the FLAGS layer), and blinksky_data.py next to this file.

Example:
    python blinksky_full_round.py --cache ../.fits_cache --out ../real_round_full

Data credit: SPHEREx/JPL/Caltech/NASA via NASA/IPAC IRSA; small-body data from NASA JPL Solar System Dynamics.
"""

import argparse
import json
import math
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import requests

import blinksky_data as bd


class TanSip:
    """Full TAN-SIP world coordinates for one SPHEREx image (pixel <-> RA/Dec), including lens distortion."""
    def __init__(self, h: dict):
        g = lambda k, d=0.0: float(h.get(k, d))
        self.crpix = (g("CRPIX1"), g("CRPIX2")); self.ra0, self.dec0 = math.radians(g("CRVAL1")), math.radians(g("CRVAL2"))
        if "CD1_1" in h:
            self.cd = np.array([[g("CD1_1"), g("CD1_2")], [g("CD2_1"), g("CD2_2")]])
        else:
            self.cd = np.array([[g("PC1_1", 1) * g("CDELT1", 1), g("PC1_2") * g("CDELT1", 1)],
                                [g("PC2_1") * g("CDELT2", 1), g("PC2_2", 1) * g("CDELT2", 1)]])
        self.cdi = np.linalg.inv(self.cd)

        def poly(prefix):
            order = int(float(h.get(prefix + "ORDER", 0)))
            return [(p, q, g(f"{prefix}{p}_{q}")) for p in range(order + 1) for q in range(order + 1 - p) if f"{prefix}{p}_{q}" in h]
        self.A, self.B, self.AP, self.BP = poly("A_"), poly("B_"), poly("AP_"), poly("BP_")

    @staticmethod
    def _sip(terms, u, v):
        out = np.zeros_like(u, dtype=float)
        for p, q, c in terms:
            out = out + c * u ** p * v ** q
        return out

    def pix2world(self, x, y):
        """0-based pixel arrays -> RA, Dec in degrees."""
        u = x + 1 - self.crpix[0]; v = y + 1 - self.crpix[1]
        u2 = u + self._sip(self.A, u, v); v2 = v + self._sip(self.B, u, v)
        xi = np.radians(self.cd[0, 0] * u2 + self.cd[0, 1] * v2); eta = np.radians(self.cd[1, 0] * u2 + self.cd[1, 1] * v2)
        d = np.cos(self.dec0) - eta * np.sin(self.dec0)
        ra = self.ra0 + np.arctan2(xi, d)
        dec = np.arctan2(np.sin(self.dec0) + eta * np.cos(self.dec0), np.hypot(xi, d))
        return np.degrees(ra) % 360, np.degrees(dec)

    def world2pix(self, ra, dec):
        """RA, Dec in degrees -> 0-based pixel arrays."""
        a, dd = np.radians(ra), np.radians(dec)
        cosc = np.sin(self.dec0) * np.sin(dd) + np.cos(self.dec0) * np.cos(dd) * np.cos(a - self.ra0)
        xi = np.degrees(np.cos(dd) * np.sin(a - self.ra0) / cosc)
        eta = np.degrees((np.cos(self.dec0) * np.sin(dd) - np.sin(self.dec0) * np.cos(dd) * np.cos(a - self.ra0)) / cosc)
        u2 = self.cdi[0, 0] * xi + self.cdi[0, 1] * eta; v2 = self.cdi[1, 0] * xi + self.cdi[1, 1] * eta
        u = u2 + self._sip(self.AP, u2, v2); v = v2 + self._sip(self.BP, u2, v2)
        return u + self.crpix[0] - 1, v + self.crpix[1] - 1


def reproject_onto(src, src_wcs, ref_wcs, shape):
    """Resample image src (with src_wcs) onto the pixel grid of ref_wcs, so the stars line up for the blink."""
    from scipy.ndimage import map_coordinates
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    ra, dec = ref_wcs.pix2world(xx, yy)
    sx, sy = src_wcs.world2pix(ra, dec)
    data = np.nan_to_num(src.astype(float), nan=float(np.nanmedian(src)))
    out = map_coordinates(data, [sy, sx], order=1, mode="constant", cval=np.nan)
    inside = (sx >= 0) & (sx <= src.shape[1] - 1) & (sy >= 0) & (sy <= src.shape[0] - 1)
    out[~inside] = np.nan
    return out, inside


def shared_rectangle(mask):
    """Biggest-ish axis-aligned rectangle fully inside mask: start from its bounding box, trim the worst edge until clean."""
    ys, xs = np.where(mask)
    if not len(ys):
        raise RuntimeError("the two photos do not overlap; pick a closer pair")
    y0, y1, x0, x1 = int(ys.min()), int(ys.max()) + 1, int(xs.min()), int(xs.max()) + 1
    while y1 - y0 > 10 and x1 - x0 > 10:
        sub = mask[y0:y1, x0:x1]
        if sub.all():
            return y0, y1, x0, x1
        bad = {"top": (~sub[0]).mean(), "bottom": (~sub[-1]).mean(), "left": (~sub[:, 0]).mean(), "right": (~sub[:, -1]).mean()}
        side = max(bad, key=bad.get)
        if side == "top": y0 += 1
        elif side == "bottom": y1 -= 1
        elif side == "left": x0 += 1
        else: x1 -= 1
    raise RuntimeError("the shared area is too small; pick a closer pair")


# SPHEREx FLAGS bits that mark pixels which would blink without anything real being there.
# Names and meanings come from the MSKN/MSKD keywords in every SPHEREx image's FLAGS header.
BAD_BITS = {0: "TRANSIENT (cosmic ray)", 1: "OVERFLOW", 2: "SUR_ERROR", 6: "NONFUNC", 9: "MISSING_DATA", 10: "HOT",
            11: "COLD", 17: "PERSIST", 19: "OUTLIER", 25: "STREAK", 26: "BLOOM"}


def clean_flagged(pixels, flags):
    """Replace flagged pixels with the median of their 5x5 neighbourhood. Real stars are left untouched."""
    from scipy.ndimage import median_filter
    mask = np.zeros(flags.shape, dtype=bool)
    for bit in BAD_BITS:
        mask |= ((flags.astype(np.int64) >> bit) & 1).astype(bool)
    data = np.nan_to_num(pixels.astype(float), nan=float(np.nanmedian(pixels)))
    smooth = median_filter(data, size=5)
    data[mask] = smooth[mask]
    return data, int(mask.sum())


def read_flags(path: Path):
    """The FLAGS layer is tile-compressed (RICE), which needs astropy to unpack."""
    from astropy.io import fits
    with fits.open(path) as f:
        return np.asarray(f["FLAGS"].data)


def to_png_shared(image, path: Path, lo: float, hi: float) -> None:
    """Same brightness stretch for both photos, so the sky does not flicker in the blink."""
    from PIL import Image
    v = np.nan_to_num(image.astype(float), nan=lo)
    s = np.arcsinh((np.clip(v, lo, hi) - lo) / max(hi - lo, 1e-9) * 10) / np.arcsinh(10)
    Image.fromarray((s * 255).astype(np.uint8)[::-1]).save(path)   # FITS rows run bottom-up


def download(url: str, path: Path) -> Path:
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    bd.log(f"    downloading {path.name} (~65 MB) ...")
    with requests.get(url, stream=True, headers=bd.HEADERS, timeout=600) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    return path


def _ra_deg(s: str) -> float:
    hh, mm, ss = (float(x) for x in s.split(":"))
    return 15 * (hh + mm / 60 + ss / 3600)


def _dec_deg(s: str) -> float:
    s = s.replace("'", " ").replace('"', " ").split()
    d, m, sec = float(s[0]), float(s[1]), float(s[2])
    sign = -1 if s[0].startswith("-") else 1
    return sign * (abs(d) + m / 60 + sec / 3600)


def known_bright(ra: float, dec: float, when: str, halfwidth: float, vmag_lim: float, cache: Path | None = None) -> list[dict]:
    """Known objects brighter than vmag_lim in a box, from NASA JPL sb_ident.

    Uses JPL's high-precision second pass and reads each object's astrometric RA/Dec directly.
    (The 'distance from centre' columns are rounded to one significant figure, e.g. "4.E3", which is
    far too coarse to place an object on a SPHEREx pixel.) The observer is Earth's centre, which is
    fine for SPHEREx in low Earth orbit (parallax well under one SPHEREx pixel for main-belt objects).
    """
    key = cache / f"jpl2_{when}_{ra:.3f}_{dec:.3f}_{halfwidth}_{vmag_lim}.json".replace(":", "-") if cache else None
    if key and key.exists():
        return json.loads(key.read_text())
    dec_s = bd._dms(dec).replace("-", "M", 1) if dec < 0 else bd._dms(dec)
    params = {"xobs": "0,0,0,0,0,0", "obs-time": when, "fov-ra-center": bd._hms(ra), "fov-dec-center": dec_s,
              "fov-ra-hwidth": halfwidth, "fov-dec-hwidth": halfwidth, "vmag-lim": vmag_lim,
              "two-pass": "true", "suppress-first-pass": "true"}
    r = requests.get(bd.JPL_SB_IDENT, params=params, headers=bd.HEADERS, timeout=600)
    r.raise_for_status()
    j = r.json()
    if "message" in j:
        raise RuntimeError(f"JPL: {j['message']}")
    out = []
    for row in j.get("data_second_pass") or []:
        try:
            vmag = float(row[6])
        except ValueError:
            vmag = None
        out.append({"name": row[0], "v_mag": vmag, "ra": _ra_deg(row[1]), "dec": _dec_deg(row[2])})
    if key:
        key.write_text(json.dumps(out))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Build a full-image BlinkSky round")
    ap.add_argument("--ra", type=float, default=237.82)
    ap.add_argument("--dec", type=float, default=-20.15)
    ap.add_argument("--gap-hours", type=float, default=3.17)
    ap.add_argument("--detector", default="D5")
    ap.add_argument("--vmag-lim", type=float, default=22)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cache, out = Path(a.cache), Path(a.out)
    cache.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    try:
        imgs = [i for i in bd.find_spherex_images(a.ra, a.dec, 0.05) if i["detector"].endswith(a.detector)]
        pairs = [(abs((j["time_mjd"] - i["time_mjd"]) * 24 - a.gap_hours), i, j) for n, i in enumerate(imgs) for j in imgs[n + 1:]]
        if not pairs:
            raise RuntimeError("no image pair found for that detector")
        _, first, second = min(pairs, key=lambda t: t[0])
        bd.log(f"pair: {first['obs_id']} -> {second['obs_id']} ({(second['time_mjd'] - first['time_mjd']) * 24:.2f} h apart)")
        frames, names, data = [], {}, []
        for img in (first, second):
            f = download(img["url"], cache / f"{img['obs_id']}{img['detector'][-2:]}.fits")
            header, pixels = bd.read_fits(f.read_bytes())["IMAGE"]
            cleaned, n_bad = clean_flagged(pixels, read_flags(f))
            bd.log(f"    {img['obs_id']}: cleaned {n_bad:,} flagged pixels (cosmic rays, hot pixels, ...)")
            data.append((img, header, TanSip(header), cleaned))
        (img_a, hdr_a, wcs_a, pix_a), (img_b, hdr_b, wcs_b, pix_b) = data
        bd.log("    lining photo B up with photo A (pixel by pixel, with lens distortion) ...")
        pix_b_on_a, inside = reproject_onto(pix_b, wcs_b, wcs_a, pix_a.shape)
        y0, y1, x0, x1 = shared_rectangle(inside)
        bd.log(f"    shared area: x {x0}-{x1}, y {y0}-{y1} of {pix_a.shape[1]}x{pix_a.shape[0]}")
        crops = [pix_a[y0:y1, x0:x1], pix_b_on_a[y0:y1, x0:x1]]
        ny, nx = crops[0].shape
        ra_c, dec_c = (float(v) for v in wcs_a.pix2world(np.array((x0 + x1) / 2), np.array((y0 + y1) / 2)))
        half = max(nx, ny) * img_a["pixel_arcsec"] / 3600 / 2 * 1.05
        lo, hi = (float(v) for v in np.nanpercentile(crops[0], [1, 99.7]))   # one stretch, taken from photo A
        for n, img in enumerate((img_a, img_b)):
            fname = f"frame_{n + 1}.png"
            to_png_shared(crops[n], out / fname, lo, hi)
            when = bd.mjd_to_jpl_time(img["time_mjd"])
            bd.log(f"    asking NASA JPL for objects brighter than V{a.vmag_lim} in a {2 * half:.1f} deg box at {when} ...")
            for o in known_bright(ra_c, dec_c, when, round(half, 2), a.vmag_lim, cache):
                x, y = (float(v) for v in wcs_a.world2pix(np.array(o["ra"]), np.array(o["dec"])))
                x, y = x - x0, y - y0
                pos = [round(x + .5, 1), round(ny - 1 - y + .5, 1)]   # PNG rows run top-down
                ok = 0 <= pos[0] < nx and 0 <= pos[1] < ny
                e = names.setdefault(o["name"], {"name": o["name"], "v_mag": o["v_mag"], "positions": [None, None]})
                e["positions"][n] = pos if ok else None
            t = datetime(1858, 11, 17) + timedelta(days=img["time_mjd"])
            frames.append({"src": fname, "time_utc": t.strftime("%Y-%m-%d %H:%M UTC"), "detector": img["detector"],
                           "wavelength_um": img["wavelength_um"], "obs_id": img["obs_id"], "width": nx, "height": ny})
        objects = sorted(names.values(), key=lambda o: (o["v_mag"] is None, o["v_mag"] or 99))
        round_ = {
            "source": "NASA SPHEREx Quick Release images via NASA/IPAC IRSA; known objects from NASA JPL sb_ident",
            "collection": "spherex_qr3", "center": {"ra": round(ra_c, 5), "dec": round(dec_c, 5)}, "full_image": True, "png_scale": 1, "aligned_to": "frame_1",
            # photo A's sky-position keywords + where the crop sits, so other tools can turn PNG pixels into RA/Dec
            "wcs": {k: hdr_a[k] for k in hdr_a if k.startswith(("CRPIX", "CRVAL", "CD", "PC", "CDELT", "A_", "B_", "AP_", "BP_")) and not k[-1] in "AW"},
            "crop": {"x0": x0, "y0": y0, "width": nx, "height": ny},
            "vmag_limit": a.vmag_lim, "pixel_arcsec": first["pixel_arcsec"],
            "width": frames[0]["width"], "height": frames[0]["height"], "frames": frames, "objects": objects,
        }
        (out / "round.json").write_text(json.dumps(round_, indent=1))
    except (requests.RequestException, RuntimeError, KeyError) as e:
        bd.log(f"FAILED: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1
    inside_both = [o for o in objects if all(o["positions"])]
    print(json.dumps({"ok": True, "out": str(out), "size_px": [round_["width"], round_["height"]],
                      "objects_listed": len(objects), "in_both_photos": len(inside_both),
                      "brightest": [(o["name"], o["v_mag"]) for o in objects[:8]]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
