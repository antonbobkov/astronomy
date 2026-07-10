"""
Compute sky positions (altitude / azimuth) for every object in a SkySafari
observing list, for a given place and time.

Configured for: Adin, CA -- 11 pm PDT on Sat Jul 11 2026 (= 2026-07-12 06:00 UTC).
Shared logic lives in skylib.py.
"""

import sys
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

import skylib

# Windows consoles default to cp1252; ensure Greek Bayer letters print cleanly.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

SKYLIST_PATH = "2026-best-objects.skylist"

# 11 pm PDT on Sat Jul 11 2026
TIME = skylib.local_pdt(2026, 7, 11, 23, 0)
LOCAL_LABEL = "11:00 pm PDT, Sat Jul 11 2026"


def fmt_radec(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    return c.to_string("hmsdms", sep=":", precision=0)


def main():
    location = skylib.make_location()
    _, objects = skylib.parse_skylist(SKYLIST_PATH)
    objects.sort(key=lambda o: (o["index"] if o["index"] is not None else 1e9))

    print("Sky positions -- %s" % skylib.SITE_NAME)
    print("Time: %s  (%s UTC)" % (LOCAL_LABEL, TIME.iso))
    print("Site: lat %+.4f, lon %+.4f, %.0f m" % (
        skylib.LAT_DEG, skylib.LON_DEG, skylib.HEIGHT_M))
    print()
    header = "%-3s %-22s %-24s %6s %6s %-4s %-6s" % (
        "#", "Object", "RA / Dec (J2000)", "Alt", "Az", "Dir", "Vis")
    print(header)
    print("-" * len(header))

    for obj in objects:
        name = skylib.display_name(obj)
        if name == "Earth":
            print("%-3d %-22s  -- skipped (observing site) --" % (obj["index"], name))
            continue
        try:
            r = skylib.object_altaz(obj, TIME, location)
            alt, az = float(r["alt"][0]), float(r["az"][0])
            ra, dec = float(r["ra"][0]), float(r["dec"][0])
            vis = "up" if alt >= 0 else "below"
            print("%-3d %-22s %-24s %6.1f %6.1f %-4s %-6s" % (
                obj["index"], name[:22], fmt_radec(ra, dec), alt, az,
                skylib.compass(az), vis))
        except Exception as e:
            print("%-3d %-22s  !! could not resolve: %s" % (
                obj["index"], name[:22], repr(e)[:70]))


if __name__ == "__main__":
    main()
