# Astronomy — SkySafari observing-list toolkit

A small Python toolkit built around a SkySafari observing list (`.skylist`). It
computes each object's position in the sky, filters the list down to what is
actually observable from a given site and time window, exports the result to a
telescope-friendly CSV, and generates beginner-friendly field notes.

Everything here is configured for a session at **Adin, CA**, on the night of
**Saturday 11 July 2026**.

## Pipeline

1. **`positions.py`** — prints the altitude/azimuth (and RA/Dec) of every object
   in the list at a given time.
2. **`filter_skylist.py`** — keeps only objects that rise at least 10° above the
   horizon between 10 pm and midnight → `2026-best-objects-filtered.skylist`.
3. **`to_csv.py`** — exports the filtered list to an `Objects4.csv`-style CSV,
   filling in object type and magnitude from external catalogs →
   `2026-best-objects-filtered.csv`.
4. **`descriptions.py`** — resolves each object to its Wikipedia article and
   assembles a notes file.

`skylib.py` is the shared library used by all of the above: skylist parsing,
coordinate/altitude resolution, and the various catalog lookups.

## Data sources

- Stars & deep-sky positions — `astropy` name resolution (CDS/Sesame)
- Planets — `astropy` built-in ephemeris
- Moons, asteroids, comets — JPL Horizons (via `astroquery`)
- Deep-sky type & magnitude — [OpenNGC](https://github.com/mattiaverga/OpenNGC)
- Star type & magnitude — SIMBAD (via `astroquery`)
- Descriptions — Wikipedia REST API

`skylib.py` also bootstraps a TLS trust bundle so the JPL Horizons endpoint
verifies correctly (it ships an incomplete certificate chain).

## Files

| File | What it is |
|------|-----------|
| `2026-best-objects.skylist` | Input observing list (88 objects) |
| `2026-best-objects-filtered.skylist` | Filtered list (62 observable objects) |
| `2026-best-objects-filtered.csv` | Filtered list in `Objects4.csv` point-to format |
| `field-notes.md` | Beginner-friendly write-up of each object |
| `Objects4.csv` | Reference file defining the CSV export format |

## Requirements

Python 3.9+, with `astropy`, `astroquery`, and `numpy`. An internet connection
is needed for name resolution, JPL Horizons, OpenNGC, and Wikipedia.

```bash
python positions.py         # positions at 11 pm
python filter_skylist.py    # write the filtered .skylist
python to_csv.py            # write the filtered .csv
```
