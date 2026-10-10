# BlinkSky

**Team Earth Intruders (Tahia Tabassum, Mustain Islam) · NASA Space Apps Challenge 2026 · Dhaka, Bangladesh**
**Challenge 12: Planet X and SPHEREx**

BlinkSky is a game where you look for asteroids in photos of the sky above your city.
The photos flip back and forth. The stars stay where they are. You tap whatever moved, and it goes into your collection.

> ## ⚠️ A note on the hackathon rules
>
> ***We are aware that NASA Space Apps does not allow teams to start building their project before the hackathon.
> This prototype was made only because the Space Apps Bangladesh organisers asked for a prototype with the
> prescreening submission, to show the concept and how the NASA data can be accessed. It is here to exhibit
> the concept to the local judges.***
>
> ***Nothing in this repository will be reused on the event day. We will build everything from scratch
> during the hackathon, based on the full challenge statement, which NASA has not published yet (due 28 October 2026).***

## Status: early prototype for the Bangladesh prescreening round

This repository is the prototype the Space Apps Bangladesh organisers asked for with the prescreening video:
a prototype showing how BlinkSky gets NASA data and how that data is used in the code.
The full game will be built during the hackathon (13–14 November 2026).

- **Real NASA data:** `data_access/blinksky_data.py` fetches real SPHEREx images from NASA/IPAC IRSA and asks NASA JPL which known objects are in the photo. Its output is in `real_round/`, and the game's **Real NASA sky** mode plays it.
- **Sample data:** the game's other rounds use **generated sample** sky photos and **sample** object facts, and every such card is labelled "Sample data".

Prescreening video: https://youtu.be/enedMO1WhUw

## How to try it

**Online:** https://skibidimama106-lang.github.io/BlinkSky/

**On your computer:** serve the folder, then visit http://localhost:8000:

```
python -m http.server 8000
```

(Double-clicking `index.html` also opens the game, but browsers block it from loading the real-data round, so **Real NASA sky** only works when the folder is served or opened online.)

1. Press **Hunt now**, then **Start** on the SPHEREx images. (A practice sky with generated images is there too, as a warm-up.)
2. Watch the images flip. Drag to look around, pinch or scroll to zoom, and tap the dot that moves.
   The mover looks exactly like the stars around it: nothing glows or points at it. Stuck? **Hint** circles an area that contains it (never the dot itself), and the card then says "with a hint". In the SPHEREx images, **Show NASA's list** circles every known object NASA JPL lists there.
3. Your catch arrives face up, with a piece of the image you found it in, a ring where you tapped, and the symbol of its type in the corner. Tap the card for its details: real facts from NASA JPL, like who an asteroid is named after.
4. Open **Deck** (bottom bar): your cards sit in one stack per type. Sort the stacks by type, rarity or source, and tap a stack to swipe through it.
5. **Search** (bottom bar) lists every known object NASA JPL found in the SPHEREx images. Tap one and the images open with that object circled. Looking never adds a card; cards only come from catching.

The first few rounds show each kind of object once. After that, rarer objects show up less often.

## The game in short

Five kinds of object to find. Each has its own colour, and every card shows its rarity as a plain word (Common, Uncommon, Rare, Very rare, Unique) plus the real number behind it. No gems, foil or levels.

| Kind | How often (rough share of known objects) |
|---|---|
| Main-belt asteroid | most movers you will see |
| Near-Earth asteroid | about 1 in 40 known asteroids |
| Comet | about 1 in 350 known objects |
| Interstellar visitor | 3 ever seen |
| Not on any list | very rare; passed on for experts to check |

### Look and feel (October 2026 update)

- Black screens, thin type (Hanken Grotesk), one colour per kind of object. The only bright thing on screen is the object.
- Cards: black metal fronts with your own photo, name, kind and three facts. One uniform card back for every card: flat purple with line drawings of things you can catch, the name top and bottom like a playing card.
- Removed from the earlier version: face-down "reveal" cards, rarity gems, rainbow foil, XP and levels, and the glowing mover. They made the game feel like a gambling app, and the glow gave the answer away. Practice rounds never show the answer; only **Real NASA sky** has "Show NASA's list", because its real asteroids are too faint to see by eye.

## How BlinkSky gets its NASA data

The folder `data_access/` contains a working script, `blinksky_data.py`, that shows the data path the game will use:

```
patch of sky (RA, Dec)
   │
   ├─▶ 1. NASA/IPAC IRSA: which SPHEREx images cover this patch?      (IVOA Simple Image Access, collection spherex_qr3)
   │
   ├─▶ 2. IRSA cutout: download only a small piece of each image     (?center=RA,Dec&size=0.1 → ~4 MB instead of ~65 MB)
   │        └─ read the FITS file and turn it into a picture for the blink
   │
   └─▶ 3. NASA JPL: which known asteroids and comets are in this patch at the time of the photo?   (sb_ident API)
            └─ a tap is checked against these positions; a mover that matches nothing is flagged for experts
```

Run it:

```
cd data_access
pip install requests numpy pillow
python blinksky_data.py --ra 10.68 --dec 41.27
```

It prints a JSON summary and saves PNG cutouts, the raw FITS cutouts, and JPL's answer into `out/`.

### What we learned from the real data (checked October 2026)

- SPHEREx pixels are **6.15 arcseconds** wide, so movers show as small dots.
- Near the Andromeda galaxy, IRSA's Quick Release 3 has **58 images within 0.05°**, spread over 18 days.
- Each exposure is recorded by several detectors at once, **each in a different infrared colour band** (for example 1.63–2.44 µm and 4.41–5.01 µm).
- Every image file also carries a **FLAGS** layer that marks bad pixels, and a **WCS-WAVE** layer giving the wavelength at each pixel.
- In three SPHEREx-covered 1° × 1° fields near the ecliptic, NASA JPL lists **116 to 208 known asteroids** each, and **none brighter than magnitude 16.5**.
- A typical asteroid moves about **30 arcseconds per hour**: 0.16 pixels in 2 minutes, about 5 pixels in an hour, about 117 pixels in a day. So the blink has to use photos taken **hours apart**, not minutes.
- Fields far from the ecliptic (like Andromeda's) have almost no known asteroids, so the game should favour patches near the ecliptic.
- SPHEREx revisits a spot in **bursts of 4 exposures about 2.4 minutes apart, then again 3 to 24 hours later**. A pair 3 hours apart gives roughly 15 pixels of motion, which is a good blink.
- A test on a 0.2° patch on the ecliptic (`--ra 237.82 --dec -20.15 --size 0.2 --gap-hours 3 --mark`): JPL predicted 6 known asteroids, all fainter than magnitude 19, and **none were visible** in single exposures. Bright asteroids are rare, so the game needs **large fields** (a full SPHEREx image is about 3.5° wide) and still has to be tested for how faint a mover SPHEREx can show.
- Single bright pixels that appear in one photo and not the next are **camera glitches** (very likely cosmic-ray hits), not asteroids. The FLAGS layer is the first filter for these.

Run the ecliptic blink test with marked known objects:

```
python blinksky_data.py --ra 237.82 --dec -20.15 --size 0.2 --gap-hours 3 --mark --jpl-halfwidth 0.1
```

This saves the two photos with JPL's predicted positions circled, plus `blink.gif`.

### Full SPHEREx image round (October 2026)

The game now plays a **whole SPHEREx image pair** (`real_round_full/`, about 2,034 × 1,809 SPHEREx pixels, roughly 3.5° of sky), in the spirit of the image sets IASC citizen scientists search: a big image you scan piece by piece, with no promise of how many movers it holds. It is shown at search sharpness (about 2–3 screen pixels per SPHEREx pixel) and zoom is capped both ways, so movers stay visible and the image never turns to blur.

`data_access/blinksky_full_round.py` builds it:

- downloads the two full images (D5, 3.8–4.4 µm, 2026-07-21 07:12 and 10:24 UTC) from NASA/IPAC IRSA,
- **cleans flagged pixels** (cosmic rays, hot pixels, persistence, streaks…, from each image's FLAGS layer) so they do not blink like movers,
- **lines photo B up with photo A pixel by pixel** using the full TAN-SIP positions in the files, and crops both to the area they share,
- uses **one brightness stretch** for both photos so the sky does not flicker,
- asks NASA JPL (sb_ident, high-precision second pass) for known objects brighter than V20, and places each one from its **precise RA/Dec**.

What we found: **21 Lutetia (V 11.2) is clearly visible and moves about 4 SPHEREx pixels in 3 hours**; known asteroids around V 17 are not visible in a single exposure in this band. So in band D5 the visibility limit lies somewhere between magnitude 11 and ~16; finding it exactly, and testing the shorter-wavelength detectors, is next.

**A bug we fixed:** JPL's "distance from centre" columns are rounded to one significant figure (e.g. `4.E3` arcseconds), which can be off by tens of SPHEREx pixels. Object positions are now computed from JPL's astrometric RA/Dec instead (in `blinksky_data.py` and `blinksky_full_round.py`). The older `real_round/` and `real_round_wide/` folders were made before this fix; the game uses them only as a fallback.

```
python blinksky_full_round.py --vmag-lim 20 --cache ../.fits_cache --out ../real_round_full
python blinksky_facts.py --round ../real_round_full
```

### Things that stay put: stars

Tapping a star is not a wrong answer. The game checks the tapped spot in both images: a star is bright in the same place in both, a mover is in one and not the other. For stars it offers **"What is it?"**: distance in light years, temperature compared with the Sun, colour and brightness, from the **Gaia DR3** catalogue served by **NASA/IPAC IRSA**. Every star you look up is kept in a **Stars seen** log in the Deck, like a Pokédex records what you have seen as well as what you have caught.

`data_access/blinksky_stars.py` finds the visible sources in the full round, turns their pixels into RA/Dec with the round's saved positions, downloads Gaia stars for the same sky (IRSA TAP) and matches them: 3,089 of the 5,000 brightest visible sources match a Gaia star.

```
python blinksky_stars.py --round ../real_round_full
```

If something changed between the images and NASA JPL has nothing there, the game flags it as an **Unlisted change** for expert review (most such flags turn out to be glitches; experts check each one).

### Real facts on every card

`data_access/blinksky_facts.py` asks the **NASA JPL Small-Body Database API** (`https://ssd-api.jpl.nasa.gov/sbdb.api`) about every object in a round and saves a short summary to `facts.json` next to the photos:
type of orbit, size, how long one trip around the Sun takes, when and by whom it was discovered, and who it is named after (when it has a name).
Where JPL has no measured size, the size is estimated from the object's brightness with a typical reflectivity and is labelled "(estimate)" on the card.
Cards show only three of these facts; the rest sit in the card's details, so players are not buried in numbers.

```
python blinksky_facts.py --round ../real_round_wide
```

### How the data gets into the game

The script can export a ready-to-play round straight into the game folder:

```
python blinksky_data.py --ra 237.82 --dec -20.15 --size 0.2 --gap-hours 3 --jpl-halfwidth 0.1 --export-round ../real_round
```

This writes `real_round/round.json` (photo times, detector, wavelength band, and every known object's position in picture pixels for each photo) plus `frame_1.png` and `frame_2.png`.

In the game, **Real NASA sky** loads that round:
- the two real SPHEREx photos blink full-screen, labelled with their real observation times,
- **Show NASA's list** circles and names every known asteroid NASA JPL says is in the photos (a data view, only in this mode; a catch made with it is marked "with NASA's list"),
- **Hint** circles an area that contains a known asteroid, without marking it,
- tapping one opens a card built from the real data (JPL name, brightness, which photos it is in, the SPHEREx band) with a crop of the real photo, marked **✓ NASA JPL** in your deck.

### How the real photos are shown

A SPHEREx pixel covers 6.15 arcseconds of sky, so the 0.2° patch in this repository is only 117 × 117 real pixels.
Enlarged to fill a phone screen as hard squares, it looks blocky. The game draws it the way survey images are usually displayed:
the empty sky is set to near-black, brightness is shown in grey, the image is enlarged smoothly, and bright stars get a soft glow.
No detail is added: every star stays exactly where SPHEREx recorded it.

The game first looks for a wider round in `real_round_wide/` (made with `--size 0.6`, so less zoom is needed) and falls back to `real_round/`.

**Wider round (`real_round_wide/`, added October 2026):** a 0.6° patch (9× the area) around the same spot, from the same morning (2026-07-21, 07:12 and 10:24 UTC, detector D5). IRSA's cutout service was down, so it was cut locally from the two full SPHEREx images. NASA JPL lists 61 known objects in it (brightest magnitude 18.6). The game loads this one first, so the real sky looks like a real star field.

The round included in this repository is the 0.2° ecliptic patch above. Its known asteroids are fainter than magnitude 19, so they are hard to see without the circles; the full game will use larger fields with brighter movers.

Data credit: SPHEREx/JPL/Caltech/NASA, via NASA/IPAC IRSA. Small-body data: NASA JPL Solar System Dynamics.

## Built with

Plain HTML, CSS and JavaScript. No libraries. All illustrations and cards are drawn with vector shapes and CSS in code.

## Use of AI

Code and illustrations were written with AI assistance (Claude, by Anthropic), directed and reviewed by the team.
