"""
BlinkSky colour fingerprint for a game round: what colour each photo saw each catchable object in.

What it does : SPHEREx sees the sky in 102 colours, and every photo catches a given spot in just one of them.
               For each object in a round that is bright enough to see, this tool goes back to the two SPHEREx
               photos of the round, cuts out the spot where NASA JPL Horizons says the object was, and measures
               (a) which colour (wavelength) that photo saw it in and (b) how bright it was.
               One photo = one slice of the object's 102-colour fingerprint. Writes fingerprint.json in the round folder.
Inputs       : --round DIR      folder holding round.json and horizons.json (run blinksky_horizons.py first)
               --max-mag M      only objects at least this bright (default 14; fainter ones do not show in one photo)
               --min-snr S      keep a slice only if the object stands out this clearly (default 5)
               --cache DIR      cutout cache (default .tmp/blinksky_explorer_cache)
Outputs      : DIR/fingerprint.json  {name: {"slices": [{frame, detector, wave_um, width_um, flux_mjy, snr, time_utc}]}}
               JSON summary on stdout; progress on stderr. Exit code 0 on success, 1 on failure.
Needs        : Python 3.10+, requests, numpy, scipy, astropy, pillow (same as blinksky_explorer.py).

Example:
    python blinksky_fingerprint.py --round ../real_round_full

Data credit: SPHEREx/JPL/Caltech/NASA via NASA/IPAC IRSA; positions from NASA JPL Horizons.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from blinksky_explorer import find_images, get_cutout, measure, log   # noqa: E402

CUT_DEG = 0.06


def main() -> int:
    ap = argparse.ArgumentParser(description="BlinkSky colour fingerprint for a round")
    ap.add_argument("--round", required=True)
    ap.add_argument("--max-mag", type=float, default=14.0)
    ap.add_argument("--min-snr", type=float, default=5.0)
    ap.add_argument("--cache", default=".tmp/blinksky_explorer_cache")
    a = ap.parse_args()
    folder = Path(a.round)
    try:
        rnd = json.loads((folder / "round.json").read_text())
        hor = json.loads((folder / "horizons.json").read_text())
    except OSError as e:
        log(f"FAILED: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1
    cache = Path(a.cache) if a.cache else None
    out, tried = {}, 0
    for o in rnd["objects"]:
        if o.get("v_mag") is None or o["v_mag"] > a.max_mag or o["name"] not in hor:
            continue
        slices = []
        for n, (frame, pos) in enumerate(zip(rnd["frames"], hor[o["name"]]["positions"])):
            if not pos:
                continue
            ra, dec = pos
            tried += 1
            try:
                match = [im for im in find_images(ra, dec)
                         if im["obs_id"] == frame["obs_id"] and im["detector"] == frame["detector"]]
                if not match:
                    log(f"    {o['name']}: photo {n + 1} not found in IRSA at that spot")
                    continue
                m = measure(get_cutout(match[0], ra, dec, CUT_DEG, cache), ra, dec, CUT_DEG)
            except Exception as e:
                log(f"    {o['name']}: photo {n + 1} failed: {e}")
                continue
            if not m or m["snr"] < a.min_snr:
                log(f"    {o['name']}: photo {n + 1} too faint to measure")
                continue
            slices.append({"frame": n, "detector": frame["detector"][-2:], "time_utc": frame["time_utc"],
                           "wave_um": m["wave_um"], "width_um": m["width_um"], "flux_mjy": m["flux_mjy"],
                           "snr": m["snr"], "flags": m["flags"]})
            log(f"    {o['name']}: photo {n + 1} saw it at {m['wave_um']} um, {m['flux_mjy']} mJy (S/N {m['snr']})")
        if slices:
            out[o["name"]] = {"slices": slices}
    (folder / "fingerprint.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({"ok": True, "round": str(folder), "objects_with_slices": len(out), "photos_tried": tried,
                      "slices": sum(len(v["slices"]) for v in out.values())}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
