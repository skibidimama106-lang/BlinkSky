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
- **Sample data:** the game's other rounds use **generated sample** sky photos and **sample** object facts, and every such card is labelled "Sample data · Prototype".

Prescreening video: https://youtu.be/enedMO1WhUw

## How to try it

**Online:** https://skibidimama106-lang.github.io/BlinkSky/

**On your computer:** serve the folder, then visit http://localhost:8000:

```
python -m http.server 8000
```

(Double-clicking `index.html` also opens the game, but browsers block it from loading the real-data round, so **Real NASA sky** only works when the folder is served or opened online.)

1. Press **Start hunt**.
2. Watch the photos flip. Drag to look around, pinch or scroll to zoom, and tap the dot that moves.
3. Your catch arrives as a face-down card. Tap to turn it over, tap again to read its field notes.
4. Open **Deck** to look through your cards: swipe the top card away to see the next one.

The first few rounds show each kind of object once. After that, rarer objects show up less often.

## The game in short

| Rarity | Object |
|---|---|
| Common | Asteroid (main belt) |
| Uncommon | Near-Earth asteroid |
| Rare | Comet |
| Legendary | Interstellar visitor |
| One of a kind | A mover that isn't on any list, passed on for experts to check |

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

### How the data gets into the game

The script can export a ready-to-play round straight into the game folder:

```
python blinksky_data.py --ra 237.82 --dec -20.15 --size 0.2 --gap-hours 3 --jpl-halfwidth 0.1 --export-round ../real_round
```

This writes `real_round/round.json` (photo times, detector, wavelength band, and every known object's position in picture pixels for each photo) plus `frame_1.png` and `frame_2.png`.

In the game, **Real NASA sky** loads that round:
- the two real SPHEREx photos blink full-screen, labelled with their real observation times,
- **Show known objects** circles where NASA JPL says each known asteroid is,
- tapping one opens a card built from the real data (JPL name, brightness, which photos it is in, the SPHEREx band), marked **Real NASA data** in your deck.

### How the real photos are shown

A SPHEREx pixel covers 6.15 arcseconds of sky, so the 0.2° patch in this repository is only 117 × 117 real pixels.
Enlarged to fill a phone screen as hard squares, it looks blocky. The game draws it the way survey images are usually displayed:
the empty sky is set to near-black, brightness is mapped to a colour scale, the image is enlarged smoothly, and bright stars get a soft glow.
No detail is added: every star stays exactly where SPHEREx recorded it.

The game first looks for a wider round in `real_round_wide/` (made with `--size 0.6`, so less zoom is needed) and falls back to `real_round/`.

The round included in this repository is the 0.2° ecliptic patch above. Its known asteroids are fainter than magnitude 19, so they are hard to see without the circles; the full game will use larger fields with brighter movers.

Data credit: SPHEREx/JPL/Caltech/NASA, via NASA/IPAC IRSA. Small-body data: NASA JPL Solar System Dynamics.

## Built with

Plain HTML, CSS and JavaScript. No libraries. All illustrations and cards are drawn with vector shapes and CSS in code.

## Use of AI

Code and illustrations were written with AI assistance (Claude, by Anthropic), directed and reviewed by the team.
