"""
Filter a SkySafari observing list down to objects that are observable during a
given window -- i.e. reach at least a minimum altitude at some point in it.

Configured for: Adin, CA, 10:00 pm -> midnight PDT on Jul 11->12, 2026.
Keeps an object if its maximum altitude across the window is >= MIN_ALT_DEG.
Writes a new .skylist in the original format & order, with DefaultIndex
renumbered contiguously from 0.

Shared logic lives in skylib.py.
"""

import sys
import numpy as np
from astropy.time import TimeDelta

import skylib

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

SRC = "2026-best-objects.skylist"
DST = "2026-best-objects-filtered.skylist"
MIN_ALT_DEG = 10.0
WINDOW_LABEL = "10:00 pm - midnight PDT, Jul 11->12 2026"


def window_times():
    """13 sample times, every 10 min, 10:00 pm -> 12:00 am PDT (inclusive)."""
    start = skylib.local_pdt(2026, 7, 11, 22, 0)          # 10 pm PDT
    minutes = np.arange(0, 121, 10)                        # 0..120 inclusive
    return start + TimeDelta(minutes * 60.0, format="sec")


def set_default_index(raw_lines, new_index):
    """Return the block's raw lines with the DefaultIndex value replaced."""
    out = []
    for line in raw_lines:
        if line.strip().startswith("DefaultIndex="):
            prefix = line[:line.index("DefaultIndex=")]      # preserve leading tab
            out.append("%sDefaultIndex=%d" % (prefix, new_index))
        else:
            out.append(line)
    return out


def main():
    location = skylib.make_location()
    times = window_times()
    skylib.ensure_horizons_ca_bundle()

    header_lines, objects = skylib.parse_skylist(SRC)
    objects.sort(key=lambda o: (o["index"] if o["index"] is not None else 1e9))

    print("Filtering %s" % SRC)
    print("Window: %s  (%s -> %s UTC)" % (WINDOW_LABEL, times[0].iso, times[-1].iso))
    print("Keep rule: peak altitude in window >= %.0f deg\n" % MIN_ALT_DEG)

    kept, dropped = [], []
    for obj in objects:
        name = skylib.display_name(obj)
        if name == "Earth":
            dropped.append((name, None, "not observable"))
            continue
        try:
            peak = float(np.max(skylib.object_altaz(obj, times, location)["alt"]))
        except Exception as e:
            # Fail open: keep anything we couldn't evaluate rather than silently drop it.
            print("  !! %-22s could not resolve, KEEPING: %s" % (name[:22], repr(e)[:60]))
            kept.append(obj)
            continue
        if peak >= MIN_ALT_DEG:
            kept.append(obj)
        else:
            dropped.append((name, peak, "peak %.1f deg" % peak))

    # Write filtered list, preserving format; renumber DefaultIndex 0..N-1.
    with open(DST, "w", encoding="utf-8", newline="\n") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for new_index, obj in enumerate(kept):
            for line in set_default_index(obj["raw_lines"], new_index):
                fh.write(line + "\n")

    print("Kept %d of %d objects; dropped %d." % (
        len(kept), len(objects), len(dropped)))
    print("Wrote %s\n" % DST)
    print("Dropped:")
    for name, peak, reason in dropped:
        print("  - %-22s %s" % (name[:22], reason))


if __name__ == "__main__":
    main()
