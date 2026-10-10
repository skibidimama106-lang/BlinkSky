"""
BlinkSky object facts: real facts about every known object in a BlinkSky round, from NASA JPL.

What it does : Reads a round's round.json, asks the NASA JPL Small-Body Database (SBDB API) about each object,
               and writes a short, player-friendly summary to facts.json in the same folder:
               orbit type, size, how long one orbit takes, when and by whom it was discovered,
               and who or what it is named after (when it has a name).
Inputs       : --round DIR       folder holding round.json (e.g. ../real_round_wide)
               --delay SECONDS   pause between requests to be polite to JPL (default 0.3)
Outputs      : DIR/facts.json, keyed by the object name exactly as it appears in round.json.
               JSON summary on stdout; progress on stderr.
Needs        : Python 3.10+, requests.

Example:
    python blinksky_facts.py --round ../real_round_wide

Data credit: NASA JPL Solar System Dynamics, Small-Body Database API (https://ssd-api.jpl.nasa.gov/doc/sbdb.html).
"""

import argparse
import json
import re
import sys
import time
from pathlib import Path

import requests

SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"
HEADERS = {"User-Agent": "BlinkSky-prototype (NASA Space Apps 2026)"}

# JPL orbit class codes -> BlinkSky's kinds of object
NEAR_EARTH = {"APO", "ATE", "AMO", "IEO"}
COMETS = {"COM", "JFc", "JFC", "HTC", "ETc", "CTc", "PAR", "HYP"}


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def lookup_key(name: str) -> str:
    """'57742 (2001 UA163)' -> '57742'; '(2022 QJ66)' -> '2022 QJ66'; '25133 Douglin (1998 SU14)' -> '25133'."""
    m = re.match(r"\s*(\d+)\b", name)
    if m:
        return m.group(1)
    m = re.search(r"\(([^)]+)\)", name)
    return m.group(1) if m else name.strip()


def first_sentence(text: str, limit: int = 160) -> str:
    # end of sentence = . ! ? after a real word (not an initial like "N."), followed by a capital letter
    s = re.split(r"(?<=[a-z\d)]{2}[.!?])\s+(?=[A-Z])", text.strip(), maxsplit=1)[0]
    return s if len(s) <= limit else s[: limit - 1].rsplit(" ", 1)[0] + "…"


def summarise(d: dict) -> dict:
    obj, disc = d.get("object", {}), d.get("discovery") or {}
    code = (obj.get("orbit_class") or {}).get("code", "")
    kind = "uncommon" if code in NEAR_EARTH else "rare" if code in COMETS or obj.get("kind", "").startswith("c") else "common"
    phys = {p["name"]: p for p in d.get("phys_par", [])}
    elements = {e["name"]: e for e in (d.get("orbit") or {}).get("elements", [])}
    size = None
    if "diameter" in phys and phys["diameter"].get("value"):
        km = float(phys["diameter"]["value"])
        size = f"{km:.1f} km across" if km >= 1 else f"{km * 1000:.0f} m across"
    elif "H" in phys and phys["H"].get("value"):
        # Not measured: estimate from brightness (absolute magnitude H) with a typical albedo of 0.14,
        # D = 1329 km / sqrt(albedo) * 10^(-H/5). Labelled as an estimate on the card.
        km = 1329 / 0.14 ** 0.5 * 10 ** (-float(phys["H"]["value"]) / 5)
        size = (f"about {km:.0f} km (estimate)" if km >= 1.5 else f"about {km * 1000:.0f} m (estimate)") if km >= 0.05 else None
    period = None
    if "per" in elements and elements["per"].get("value"):
        yrs = float(elements["per"]["value"]) / 365.25
        period = f"{yrs:.1f} years" if yrs < 100 else f"{yrs:,.0f} years"
    named = first_sentence(disc["citation"]) if disc.get("citation") else None
    year = (disc.get("date") or "")[:4] or None
    who = disc.get("who") or None
    if who and who.count(",") >= 2:          # "van Houten, C. J., van Houten-Groeneveld, I., ..." -> "van Houten et al."
        who = who.split(",")[0].strip() + " et al."
    return {
        "kind": kind,
        "orbit_class": (obj.get("orbit_class") or {}).get("name"),
        "neo": bool(obj.get("neo")), "pha": bool(obj.get("pha")),
        "size": size, "period": period,
        "discovered": " · ".join(x for x in (year, who) if x) or None,
        "named_after": named,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch NASA JPL facts for every object in a BlinkSky round")
    ap.add_argument("--round", required=True)
    ap.add_argument("--delay", type=float, default=0.3)
    a = ap.parse_args()
    folder = Path(a.round)
    try:
        objects = json.loads((folder / "round.json").read_text())["objects"]
    except (OSError, KeyError, json.JSONDecodeError) as e:
        log(f"FAILED reading {folder / 'round.json'}: {e}")
        print(json.dumps({"ok": False, "error": str(e)}))
        return 1
    facts, missing = {}, []
    for n, o in enumerate(objects, 1):
        key = lookup_key(o["name"])
        try:
            r = requests.get(SBDB, params={"sstr": key, "phys-par": 1, "discovery": 1}, headers=HEADERS, timeout=60)
            r.raise_for_status()
            d = r.json()
            if "object" not in d:
                raise ValueError(d.get("message", "no match"))
            facts[o["name"]] = summarise(d)
            log(f"{n}/{len(objects)} {o['name']}: ok")
        except (requests.RequestException, ValueError) as e:
            missing.append(o["name"])
            log(f"{n}/{len(objects)} {o['name']}: skipped ({e})")
        time.sleep(a.delay)
    (folder / "facts.json").write_text(json.dumps(facts, indent=1, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"ok": True, "written": str(folder / "facts.json"), "objects": len(objects), "with_facts": len(facts),
                      "named": sum(1 for f in facts.values() if f["named_after"]), "with_size": sum(1 for f in facts.values() if f["size"]),
                      "missing": missing}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
