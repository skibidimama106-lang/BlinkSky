"""
BlinkSky data access prototype: how BlinkSky gets its NASA data.

What it does : For one patch of sky, it
               1. finds SPHEREx images of that patch in NASA/IPAC IRSA (IVOA Simple Image Access),
               2. downloads small cutouts (not the ~65 MB full files) and turns them into PNG pictures,
               3. asks NASA JPL which known asteroids and comets are in that patch at a given time.
Inputs       : --ra DEG --dec DEG     centre of the patch (degrees)
               --radius DEG           search radius for images (default 0.05)
               --size DEG             cutout width (default 0.1, about 58 SPHEREx pixels)
               --collection NAME      IRSA collection (default spherex_qr3; spherex_qr2 for older data)
               --cutouts N            how many cutouts to download (default 2)
               --gap-hours H          instead, pick two photos about H hours apart (same detector) for a blink pair
               --mark                 circle the known objects JPL predicts on each cutout, and save blink.gif
               --export-round DIR     also write a ready-to-play BlinkSky round (round.json + frame PNGs) into DIR
               --jpl-time T           time for the known-object check, e.g. 2026-07-20_06:00:00
                                      (default: the first image's observation time)
               --jpl-halfwidth DEG    half-width of the known-object search box (default 0.3)
               --out DIR              where to write files (default ./out)
Outputs      : JSON summary on stdout; PNG cutouts, raw FITS cutouts and the JPL answer in --out.
               Progress messages go to stderr.
Needs        : Python 3.10+, requests, numpy, pillow. No astronomy library required.

Example:
    python blinksky_data.py --ra 10.68 --dec 41.27
    python blinksky_data.py --ra 300.0 --dec -20.5 --collection spherex_qr3 --cutouts 4

Data credit: SPHEREx/JPL/Caltech/NASA via NASA/IPAC IRSA; small-body data from NASA JPL Solar System Dynamics.
"""

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import requests
from PIL import Image

IRSA_SIA = "https://irsa.ipac.caltech.edu/SIA"
JPL_SB_IDENT = "https://ssd-api.jpl.nasa.gov/sb_ident.api"
VOT = "{http://www.ivoa.net/xml/VOTable/v1.3}"
# One request per connection: IRSA occasionally garbles a reused keep-alive connection.
HEADERS = {"Connection": "close", "User-Agent": "BlinkSky-prototype (NASA Space Apps 2026)"}


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


# ---------- 1. find SPHEREx images ----------

def find_spherex_images(ra: float, dec: float, radius: float = 0.05, collection: str = "spherex_qr3", limit: int = 500) -> list[dict]:
    """Ask IRSA which SPHEREx images cover a circle on the sky. Returns one dict per image, oldest first."""
    r = requests.get(IRSA_SIA, params={"COLLECTION": collection, "POS": f"CIRCLE {ra} {dec} {radius}", "MAXREC": limit},
                     headers=HEADERS, timeout=180)
    r.raise_for_status()
    root = ET.fromstring(r.content)
    fields = [f.get("name") for f in root.iter(VOT + "FIELD")]
    rows = [dict(zip(fields, [(td.text or "") for td in tr])) for tr in root.iter(VOT + "TR")]
    images = [{
        "obs_id": x["obs_id"],
        "time_mjd": float(x["t_min"]),
        "detector": x["energy_bandpassname"],          # SPHEREx has 6 detectors, each a different colour band
        "wavelength_um": [round(float(x["em_min"]) * 1e6, 3), round(float(x["em_max"]) * 1e6, 3)],
        "pixel_arcsec": float(x["s_pixel_scale"]),
        "full_file_kb": int(float(x["access_estsize"] or 0)),
        "url": x["access_url"],
    } for x in rows if x.get("dataproduct_subtype") == "science"]
    return sorted(images, key=lambda i: i["time_mjd"])


# ---------- 2. download a cutout and make a picture ----------

def get_cutout(url: str, ra: float, dec: float, size: float = 0.1) -> bytes:
    """IRSA cuts a small piece out of the full image on its side, so we download kilobytes, not ~65 MB."""
    r = requests.get(url, params={"center": f"{ra},{dec}", "size": str(size)}, headers=HEADERS, timeout=180)
    r.raise_for_status()
    return r.content


def read_fits(raw: bytes) -> dict:
    """Minimal FITS reader for SPHEREx image files: returns {EXTNAME: (header_dict, array_or_None)}."""
    pos, hdus = 0, {}
    while pos < len(raw):
        cards, done = {}, False
        while not done and pos < len(raw):
            block = raw[pos:pos + 2880]
            pos += 2880
            for i in range(0, 2880, 80):
                card = block[i:i + 80].decode("ascii", "replace")
                key = card[:8].strip()
                if key == "END":
                    done = True
                    break
                if card[8:10] == "= ":
                    cards[key] = card[10:].split("/")[0].strip().strip("'").strip()
        if not cards:
            break
        bitpix, naxis = int(cards.get("BITPIX", 8)), int(cards.get("NAXIS", 0))
        dims = [int(cards[f"NAXIS{i}"]) for i in range(1, naxis + 1)]
        nbytes = (abs(bitpix) // 8) * int(np.prod(dims)) if naxis else 0
        is_table = cards.get("XTENSION", "").startswith("BINTABLE")
        if is_table:
            nbytes += int(cards.get("PCOUNT", 0))
        data = None
        if naxis and not is_table:
            dtype = {8: ">u1", 16: ">i2", 32: ">i4", -32: ">f4", -64: ">f8"}[bitpix]
            data = np.frombuffer(raw[pos:pos + nbytes], dtype=dtype).reshape(dims[::-1])
        hdus[cards.get("EXTNAME", "PRIMARY")] = (cards, data)
        pos += -(-nbytes // 2880) * 2880
    return hdus


def to_png(image: np.ndarray, path: Path, scale: int = 6) -> None:
    """Stretch the brightness so faint things show, then enlarge with crisp pixels."""
    v = image.astype(float)
    v = np.nan_to_num(v, nan=float(np.nanmedian(v)))
    lo, hi = np.percentile(v, [1, 99.7])
    s = np.arcsinh((np.clip(v, lo, hi) - lo) / max(hi - lo, 1e-9) * 10) / np.arcsinh(10)
    pic = Image.fromarray((s * 255).astype(np.uint8)[::-1])  # FITS rows run bottom-up
    pic.resize((pic.width * scale, pic.height * scale), Image.NEAREST).save(path)


def world_to_pixel(header: dict, ra: float, dec: float) -> tuple[float, float]:
    """Gnomonic (TAN) projection using the image's CRVAL/CRPIX/CD keywords. Ignores the small SIP distortion terms."""
    import math
    ra0, dec0 = math.radians(float(header["CRVAL1"])), math.radians(float(header["CRVAL2"]))
    a, d = math.radians(ra), math.radians(dec)
    cosc = math.sin(dec0) * math.sin(d) + math.cos(dec0) * math.cos(d) * math.cos(a - ra0)
    xi = math.degrees(math.cos(d) * math.sin(a - ra0) / cosc)
    eta = math.degrees((math.cos(dec0) * math.sin(d) - math.sin(dec0) * math.cos(d) * math.cos(a - ra0)) / cosc)
    if "CD1_1" in header:
        c11, c12, c21, c22 = (float(header.get(k, 0)) for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2"))
    else:
        d1, d2 = float(header["CDELT1"]), float(header["CDELT2"])
        c11, c12, c21, c22 = (float(header.get(k, dflt)) for k, dflt in (("PC1_1", 1), ("PC1_2", 0), ("PC2_1", 0), ("PC2_2", 1)))
        c11, c12, c21, c22 = c11 * d1, c12 * d1, c21 * d2, c22 * d2
    det = c11 * c22 - c12 * c21
    x = (c22 * xi - c12 * eta) / det
    y = (-c21 * xi + c11 * eta) / det
    return x + float(header["CRPIX1"]) - 1, y + float(header["CRPIX2"]) - 1  # 0-based pixel


def mark_png(png: Path, header: dict, objects: list, centre: tuple, scale: int = 6) -> int:
    """Circle each known object on the enlarged PNG. Returns how many land inside the picture."""
    from PIL import ImageDraw
    import math
    pic = Image.open(png).convert("RGB"); draw = ImageDraw.Draw(pic)
    ny = pic.height // scale; inside = 0
    for o in objects:
        ra, dec = o["ra_deg"], o["dec_deg"]
        x, y = world_to_pixel(header, ra, dec)
        px, py = (x + .5) * scale, (ny - 1 - y + .5) * scale  # the PNG is flipped vertically
        if 0 <= px < pic.width and 0 <= py < pic.height:
            inside += 1
            draw.ellipse([px - 26, py - 26, px + 26, py + 26], outline=(195, 179, 255), width=3)
            draw.text((px + 30, py - 10), f"{o['name']}  V{o['v_mag']}", fill=(255, 210, 63))
    pic.save(png)
    return inside


# ---------- 3. which known objects are there? ----------

def mjd_to_jpl_time(mjd: float) -> str:
    from datetime import datetime, timedelta
    t = datetime(1858, 11, 17) + timedelta(days=mjd)
    return t.strftime("%Y-%m-%d_%H:%M:%S")


def _hms(ra: float) -> str:
    h = ra / 15; hh = int(h); m = (h - hh) * 60; mm = int(m)
    return f"{hh:02d}-{mm:02d}-{(m - mm) * 60:05.2f}"


def _dms(dec: float) -> str:
    sign = "-" if dec < 0 else ""; d = abs(dec); dd = int(d); m = (d - dd) * 60; mm = int(m)
    return f"{sign}{dd:02d}-{mm:02d}-{(m - mm) * 60:04.1f}"


def _ra_deg(s: str) -> float:
    hh, mm, ss = (float(x) for x in s.split(":"))
    return 15 * (hh + mm / 60 + ss / 3600)


def _dec_deg(s: str) -> float:
    s = s.replace("'", " ").replace('"', " ").split()
    sign = -1 if s[0].startswith("-") else 1
    return sign * (abs(float(s[0])) + float(s[1]) / 60 + float(s[2]) / 3600)


def known_objects(ra: float, dec: float, when: str, halfwidth: float = 0.3) -> dict:
    """NASA JPL sb_ident: known asteroids/comets inside a box on the sky at a given time (geocentric view)."""
    params = {"mpc-code": "500", "obs-time": when, "fov-ra-center": _hms(ra), "fov-dec-center": _dms(dec),
              "fov-ra-hwidth": halfwidth, "fov-dec-hwidth": halfwidth}
    r = requests.get(JPL_SB_IDENT, params=params, headers=HEADERS, timeout=300)
    r.raise_for_status()
    j = r.json()
    rows = j.get("data_first_pass") or []
    # The first pass is a rough pre-filter that can return objects far outside the box; keep only those inside.
    box = halfwidth * 3600
    rows = [row for row in rows if abs(float(row[3])) <= box and abs(float(row[4])) <= box]
    objects = []
    for row in rows:
        try:
            vmag = float(row[6])
        except ValueError:
            vmag = None
        # Place objects with the astrometric RA/Dec. The "distance from centre" columns are rounded to one
        # significant figure by JPL (e.g. "4.E3"), which can be off by tens of SPHEREx pixels.
        objects.append({"name": row[0], "ra": row[1], "dec": row[2], "ra_deg": _ra_deg(row[1]), "dec_deg": _dec_deg(row[2]), "v_mag": vmag,
                        "offset_arcsec": [float(row[3]), float(row[4])],
                        "rate_arcsec_per_hour": [float(row[7]), float(row[8])]})
    objects.sort(key=lambda o: (o["v_mag"] is None, o["v_mag"] or 99))
    return {"time": when, "count": len(objects), "checked_by_jpl": j.get("n_first_pass", 0), "objects": objects, "raw": j}


# ---------- 4. hand the data to the game ----------

def export_round(folder: Path, frames: list, a) -> None:
    """Write round.json + frame_N.png: real SPHEREx photos and where NASA JPL says each known object is, in picture pixels."""
    import math
    from datetime import datetime, timedelta
    folder.mkdir(parents=True, exist_ok=True)
    scale = 6
    out_frames, names = [], {}
    for n, (img, header, pixels, ko) in enumerate(frames):
        fname = f"frame_{n + 1}.png"
        to_png(pixels, folder / fname, scale)
        ny, nx = pixels.shape
        t = datetime(1858, 11, 17) + timedelta(days=img["time_mjd"])
        out_frames.append({"src": fname, "time_utc": t.strftime("%Y-%m-%d %H:%M UTC"), "detector": img["detector"],
                           "wavelength_um": img["wavelength_um"], "obs_id": img["obs_id"]})
        for o in ko["objects"]:
            ra, dec = o["ra_deg"], o["dec_deg"]
            x, y = world_to_pixel(header, ra, dec)
            pos = [round((x + .5) * scale, 1), round((ny - 1 - y + .5) * scale, 1)]
            inside = 0 <= pos[0] < nx * scale and 0 <= pos[1] < ny * scale
            e = names.setdefault(o["name"], {"name": o["name"], "v_mag": o["v_mag"], "positions": [None] * len(frames)})
            e["positions"][n] = pos if inside else None
    w, h = frames[0][2].shape[1] * scale, frames[0][2].shape[0] * scale
    round_ = {
        "source": "NASA SPHEREx Quick Release images via NASA/IPAC IRSA; known objects from NASA JPL sb_ident",
        "collection": a.collection, "center": {"ra": a.ra, "dec": a.dec}, "size_deg": a.size,
        "pixel_arcsec": frames[0][0]["pixel_arcsec"], "width": w, "height": h,
        "frames": out_frames,
        "objects": sorted(names.values(), key=lambda o: (o["v_mag"] is None, o["v_mag"] or 99)),
    }
    (folder / "round.json").write_text(json.dumps(round_, indent=1))
    log(f"    exported a playable round to {folder} ({len(out_frames)} photos, {len(names)} known objects)")


# ---------- demo ----------

def main() -> int:
    ap = argparse.ArgumentParser(description="BlinkSky data access prototype")
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--radius", type=float, default=0.05)
    ap.add_argument("--size", type=float, default=0.1)
    ap.add_argument("--collection", default="spherex_qr3")
    ap.add_argument("--cutouts", type=int, default=2)
    ap.add_argument("--gap-hours", type=float)
    ap.add_argument("--mark", action="store_true")
    ap.add_argument("--export-round")
    ap.add_argument("--jpl-time")
    ap.add_argument("--jpl-halfwidth", type=float, default=0.3)
    ap.add_argument("--out", default="out")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    try:
        log(f"1/3 Searching IRSA {a.collection} around RA {a.ra}, Dec {a.dec} ...")
        images = find_spherex_images(a.ra, a.dec, a.radius, a.collection)
        log(f"    {len(images)} SPHEREx images found")
        saved, frames_for_round = [], []
        # Pick images far apart in time: asteroids move ~5 pixels per hour, so minutes apart shows nothing.
        picks = images[:: max(1, len(images) // max(1, a.cutouts))][: a.cutouts] if images else []
        if a.gap_hours and images:
            pairs = [(abs((j["time_mjd"] - i["time_mjd"]) * 24 - a.gap_hours), i, j)
                     for n, i in enumerate(images) for j in images[n + 1:] if i["detector"] == j["detector"]]
            if pairs:
                _, i, j = min(pairs, key=lambda t: t[0]); picks = [i, j]
                log(f"    blink pair: {i['detector']}, {(j['time_mjd'] - i['time_mjd']) * 24:.2f} hours apart")
        for n, img in enumerate(picks, 1):
            log(f"2/3 Cutout {n}/{len(picks)}: {img['detector']} at MJD {img['time_mjd']:.4f} ...")
            raw = get_cutout(img["url"], a.ra, a.dec, a.size)
            base = f"cutout_{n}_{img['detector']}_mjd{img['time_mjd']:.4f}".replace(".", "p")
            fits_path, png_path = out / (base + ".fits"), out / (base + ".png")
            fits_path.write_bytes(raw)
            hdus = read_fits(raw)
            header, pixels = hdus["IMAGE"]
            to_png(pixels, png_path)
            marked, ko = None, None
            if a.export_round:
                ko = known_objects(a.ra, a.dec, mjd_to_jpl_time(img["time_mjd"]), a.size / 2)
                frames_for_round.append((img, header, pixels, ko))
            if a.mark:
                when_i = mjd_to_jpl_time(img["time_mjd"])
                log(f"    asking NASA JPL for known objects at {when_i} ...")
                ko = ko or known_objects(a.ra, a.dec, when_i, a.size / 2)
                marked = mark_png(png_path, header, ko["objects"], (a.ra, a.dec))
                log(f"    circled {marked} known objects on the picture")
            saved.append({"known_objects_circled": marked, "png": str(png_path), "fits_kb": len(raw) // 1024, "shape": list(pixels.shape),
                          "date_obs": header.get("DATE-OBS"), "detector": img["detector"], "wavelength_um": img["wavelength_um"],
                          "extensions": list(hdus)})
        if a.mark and len(saved) >= 2:
            frames = [Image.open(x["png"]).convert("P", palette=Image.ADAPTIVE) for x in saved]
            frames[0].save(out / "blink.gif", save_all=True, append_images=frames[1:], duration=700, loop=0)
            log(f"    saved {out / 'blink.gif'}")
        if a.export_round and frames_for_round:
            export_round(Path(a.export_round), frames_for_round, a)
        when = a.jpl_time or (mjd_to_jpl_time(images[0]["time_mjd"]) if images else None)
        jpl = None
        if when:
            log(f"3/3 Asking NASA JPL for known objects at {when} ...")
            jpl = known_objects(a.ra, a.dec, when, a.jpl_halfwidth)
            (out / "jpl_known_objects.json").write_text(json.dumps(jpl["raw"], indent=1))
            log(f"    {jpl['count']} known objects inside the box (JPL pre-checked {jpl['checked_by_jpl']})")
    except requests.RequestException as e:
        log(f"FAILED talking to a NASA service: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1

    summary = {
        "ok": True,
        "patch": {"ra": a.ra, "dec": a.dec, "radius_deg": a.radius},
        "spherex_images_found": len(images),
        "detectors_seen": sorted({i["detector"] for i in images}),
        "time_span_days": round(images[-1]["time_mjd"] - images[0]["time_mjd"], 2) if images else 0,
        "cutouts": saved,
        "known_objects": None if not jpl else {
            "time": jpl["time"], "count": jpl["count"],
            "brighter_than_v18": sum(1 for o in jpl["objects"] if o["v_mag"] is not None and o["v_mag"] <= 18),
            "brightest": jpl["objects"][:5],
        },
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
