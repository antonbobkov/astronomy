"""
Export a SkySafari observing list to the Objects4.csv point-to format:

    <title line>
    N_RH_RM_DD_DM
    <id> <Type> <Mag>,RH,RM,DD,DM
    ...

where RH/RM are RA hours / minutes (RM to 1 decimal) and DD/DM are Dec
degrees / minutes (signed).  Type and magnitude are looked up:
  * deep-sky objects   -> OpenNGC (Type, integrated V-mag, Hubble morphology)
  * stars              -> SIMBAD (object type, V-mag)
  * solar-system bodies-> JPL Horizons (apparent magnitude) + skylib type map

Shared logic lives in skylib.py.
"""

import sys
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

import skylib

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

SRC = "2026-best-objects-filtered.skylist"
DST = "2026-best-objects-filtered.csv"
TITLE = "2026 Best Objects (Adin, 11 Jul 2026) - by constellation, west first"
SCHEMA = "N_RH_RM_DD_DM"

# Reference epoch (window midpoint) for solar-system positions.
REF_TIME = skylib.local_pdt(2026, 7, 11, 23, 0)

# --- OpenNGC type code -> Objects4 vocabulary ------------------------------
OPENNGC_TYPE = {
    "OCl": "Open Cluster", "GCl": "Globular Cluster", "Cl+N": "Open Cluster",
    "*Ass": "Open Cluster",
    "PN": "Planetary Nebula", "SNR": "Supernova Remnant", "DrkN": "Dark Nebula",
    "Neb": "Bright Nebula", "HII": "Bright Nebula", "RfN": "Bright Nebula",
    "EmN": "Bright Nebula", "SFR": "Bright Nebula",
    "G": "Galaxy", "GPair": "Galaxy", "GGroup": "Galaxy", "GTrpl": "Galaxy",
    "GClstr": "Galaxy",
    "**": "Double Star", "*": "Star",
}


def galaxy_word(hubble):
    """Refine a galaxy's type from its Hubble class (e.g. 'SBb', 'E2', 'Irr')."""
    h = (hubble or "").strip()
    if not h:
        return "Galaxy"
    if h[0] == "E" or h.startswith("S0"):
        return "Elliptical Galaxy"
    if h[0] == "S":
        return "Spiral Galaxy"
    if h[0] == "I":
        return "Irregular Galaxy"
    return "Galaxy"


# --- SIMBAD otype -> word (for stars, and DSO fallback) --------------------
DOUBLE_OTYPES = {"**", "EB*", "SB*"}
VARIABLE_OTYPES = {"V*", "Al*", "Ce*", "cC*", "Mi*", "sr*", "RR*", "LP*", "RS*",
                   "BY*", "SX*", "Pu*", "Ir*", "bC*", "a2*", "Ro*", "gD*"}
SIMBAD_DSO_TYPE = {
    "OpC": "Open Cluster", "Cl*": "Open Cluster", "As*": "Open Cluster",
    "GlC": "Globular Cluster", "PN": "Planetary Nebula", "SNR": "Supernova Remnant",
    "DNe": "Dark Nebula", "HII": "Bright Nebula", "RNe": "Bright Nebula",
    "Cld": "Bright Nebula", "GNe": "Bright Nebula", "ISM": "Bright Nebula",
    "EmO": "Bright Nebula", "MoC": "Bright Nebula",
}


def star_word(otype):
    if otype in DOUBLE_OTYPES:
        return "Double Star"
    if otype in VARIABLE_OTYPES:
        return "Variable Star"
    return "Star"


# --- Row label -------------------------------------------------------------
LABEL_PREFERENCE = ["M ", "NGC ", "IC ", "C ", "Cr ", "Mel ", "HR ", "HD ",
                    "HIP ", "SAO "]


def primary_label(obj):
    if obj["type"] == 1:                       # solar-system: use the name
        return obj["common_names"][0] if obj["common_names"] else obj["catalog_numbers"][0]
    cats = obj["catalog_numbers"]
    # Bayer designation ("$a Lyr") is the nicest star id -> decode to Greek
    bayer = [c for c in cats if c.startswith("$")]
    flamsteed = [c for c in cats if c[:1].isdigit()]
    for pref in LABEL_PREFERENCE[:6]:          # DSO catalogs first
        for c in cats:
            if c.startswith(pref):
                return c
    if bayer:
        return skylib.decode_bayer_ascii(bayer[0])
    if flamsteed:
        return flamsteed[0]
    for pref in LABEL_PREFERENCE[6:]:          # HR/HD/HIP/SAO
        for c in cats:
            if c.startswith(pref):
                return c
    if obj["common_names"]:
        return obj["common_names"][0]
    return cats[0] if cats else "?"


def fmt_mag(mag):
    if mag is None or (isinstance(mag, float) and not np.isfinite(mag)):
        return "--"
    return "%.1f" % mag


def type_and_mag(obj):
    """Return (type_word, magnitude_or_None) for any object."""
    if obj["type"] == 1:
        return skylib.solar_type_mag(obj, REF_TIME)

    if obj["type"] == 4:                        # deep-sky: OpenNGC first
        row = skylib.opengc_row(obj)
        if row is not None:
            code = (row.get("Type") or "").strip()
            vmag = (row.get("V-Mag") or "").strip()
            bmag = (row.get("B-Mag") or "").strip()
            mag = float(vmag) if vmag else (float(bmag) if bmag else None)
            if code == "G":
                word = galaxy_word(row.get("Hubble"))
            else:
                word = OPENNGC_TYPE.get(code, "Deep Sky")
            return word, mag

    # stars, and DSOs not in OpenNGC -> SIMBAD
    otype, morph, vmag = skylib.simbad_type_mag(obj)
    if otype is None:
        return ("Star" if obj["type"] == 2 else "Deep Sky"), None
    if obj["type"] == 2:
        return star_word(otype), vmag
    # DSO fallback via SIMBAD
    if otype in ("G", "GiG", "GiC", "SyG", "Sy1", "Sy2", "AGN", "LIN", "rG",
                 "H2G", "EmG", "GiP", "IG", "LSB", "BiC", "GrG"):
        return galaxy_word(morph), vmag
    return SIMBAD_DSO_TYPE.get(otype, "Deep Sky"), vmag


def radec_fields(ra_deg, dec_deg):
    """Return (RH, RM_str, DD_signed_str, DM_str) for the CSV."""
    # RA -> hours / minutes (1 decimal), with 60-carry
    hours = (ra_deg % 360.0) / 15.0
    rh = int(hours)
    rm = round((hours - rh) * 60.0, 1)
    if rm >= 60.0:
        rm -= 60.0
        rh += 1
    rh %= 24
    rm_str = ("%.1f" % rm).rstrip("0").rstrip(".")   # drop trailing .0 like Objects4

    # Dec -> signed degrees / minutes (integer), with 60-carry
    sign = "-" if dec_deg < 0 else ""
    ad = abs(dec_deg)
    dd = int(ad)
    dm = int(round((ad - dd) * 60.0))
    if dm >= 60:
        dm -= 60
        dd += 1
    return rh, rm_str, "%s%d" % (sign, dd), "%d" % dm


def order_by_proximity(coords, start):
    """Order a group of sky positions into a short 'slew' path.

    Nearest-neighbor open path seeded at index `start`, then 2-opt cleanup that
    keeps `start` first. `coords` is a SkyCoord array; returns a list of indices.
    Costs are great-circle angular separations (degrees).
    """
    n = len(coords)
    if n <= 2:
        return [start] + [i for i in range(n) if i != start]

    dist = np.vstack([coords[i].separation(coords).deg for i in range(n)])

    def path_len(order):
        return float(sum(dist[order[i], order[i + 1]] for i in range(len(order) - 1)))

    # nearest-neighbor
    unvisited = set(range(n))
    order = [start]
    unvisited.discard(start)
    while unvisited:
        last = order[-1]
        order.append(min(unvisited, key=lambda j: dist[last, j]))
        unvisited.discard(order[-1])

    # 2-opt, keeping order[0] (the westernmost seed) fixed
    improved = True
    while improved:
        improved = False
        for i in range(1, n - 1):
            for k in range(i + 1, n):
                cand = order[:i] + order[i:k + 1][::-1] + order[k + 1:]
                if path_len(cand) + 1e-9 < path_len(order):
                    order, improved = cand, True
    return order


def main():
    location = skylib.make_location()
    skylib.ensure_horizons_ca_bundle()
    _, objects = skylib.parse_skylist(SRC)

    # Pass 1: resolve each object's label, type/mag, and position.
    items, notes = [], []
    for obj in objects:
        label = primary_label(obj)
        tword, mag = type_and_mag(obj)
        pos = skylib.object_altaz(obj, REF_TIME, location)
        ra, dec = float(pos["ra"][0]), float(pos["dec"][0])
        # Prefer a colourful common name over the dry type word; fall back to the
        # type for catalog-only objects. (Solar-system label already IS the name.)
        nickname = obj["common_names"][0] if (obj["type"] != 1 and obj["common_names"]) else None
        descriptor = nickname if (nickname and nickname != label) else tword
        items.append({"label": label, "descriptor": descriptor, "mag": mag,
                      "ra": ra, "dec": dec, "fields": radec_fields(ra, dec)})
        if tword == "Deep Sky" or mag is None:
            notes.append("  %-24s -> %s %s" % (label, tword, fmt_mag(mag)))

    ras = np.array([it["ra"] for it in items])
    decs = np.array([it["dec"] for it in items])
    cons = skylib.constellations(ras, decs)
    has = skylib.hour_angle_hours(ras, REF_TIME)
    for it, con, ha in zip(items, cons, has):
        it["con"], it["ha"] = str(con), float(ha)

    # Pass 2: group by constellation; order groups west-first (descending mean
    # hour angle) and each group into a short slew path.
    groups = {}
    for idx, it in enumerate(items):
        groups.setdefault(it["con"], []).append(idx)

    rows = []
    for con in sorted(groups, key=lambda c: -np.mean([items[i]["ha"] for i in groups[c]])):
        idxs = np.array(groups[con])
        coords = SkyCoord(ra=ras[idxs] * u.deg, dec=decs[idxs] * u.deg)
        start = int(np.argmax([items[i]["ha"] for i in idxs]))   # westernmost first
        for local in order_by_proximity(coords, start):
            it = items[idxs[local]]
            rh, rm, dd, dm = it["fields"]
            n_field = "%s %s (%s) %s" % (it["label"], it["descriptor"], it["con"], fmt_mag(it["mag"]))
            rows.append("%s,%d,%s,%s,%s" % (n_field, rh, rm, dd, dm))

    with open(DST, "w", encoding="ascii", newline="") as fh:
        fh.write(TITLE + "\r\n")
        fh.write(SCHEMA + "\r\n")
        for r in rows:
            fh.write(r + "\r\n")

    print("Wrote %s  (%d objects, %d constellations)" % (DST, len(rows), len(groups)))
    if notes:
        print("\nObjects with no magnitude or a coarse type (worth a manual check):")
        print("\n".join(notes))


if __name__ == "__main__":
    main()
