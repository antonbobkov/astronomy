"""
Shared library for working with SkySafari observing lists (.skylist):

  * parsing the file format (preserving raw blocks for re-emission)
  * resolving each object's coordinates from its catalog identifiers
  * computing altitude / azimuth at one or many times for a site
  * a TLS trust-bundle bootstrap so JPL Horizons queries verify correctly

Data sources:
  * stars & deep-sky objects   -> SkyCoord.from_name (CDS/Sesame, needs internet)
  * major planets              -> astropy builtin ephemeris
  * moons / asteroids / comets -> JPL Horizons via astroquery
"""

import warnings
warnings.filterwarnings("ignore")

import os
import ssl
import socket
import urllib.request

import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_body

# ---------------------------------------------------------------------------
# Observing site
# ---------------------------------------------------------------------------
SITE_NAME = "Adin, CA"
LAT_DEG = 41.1969
LON_DEG = -120.9472
HEIGHT_M = 1290.0

PDT_OFFSET_H = -7  # Pacific Daylight Time = UTC-7 (July)


def make_location():
    return EarthLocation(lat=LAT_DEG * u.deg, lon=LON_DEG * u.deg, height=HEIGHT_M * u.m)


def local_pdt(year, month, day, hour, minute=0):
    """A PDT wall-clock time returned as a UTC astropy Time (PDT = UTC-7)."""
    return Time("%04d-%02d-%02d %02d:%02d:00" % (year, month, day, hour, minute),
                scale="utc") - TimeDelta(PDT_OFFSET_H * 3600, format="sec")


def constellations(ra_deg, dec_deg):
    """3-letter IAU constellation abbreviation(s) for the given coordinates.

    Accepts scalars or arrays; returns a numpy array of short names.
    """
    ra = np.atleast_1d(np.asarray(ra_deg, dtype=float))
    dec = np.atleast_1d(np.asarray(dec_deg, dtype=float))
    coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    return np.atleast_1d(coords.get_constellation(short_name=True))


def hour_angle_hours(ra_deg, time):
    """Local hour angle in hours, wrapped to (-12, +12].

    HA = local sidereal time - RA, at the site's longitude. Positive = west of
    the meridian (past transit, setting); negative = east (rising).
    """
    lst = time.sidereal_time("apparent", longitude=LON_DEG * u.deg).hour
    ra_h = np.atleast_1d(np.asarray(ra_deg, dtype=float)) / 15.0
    ha = (lst - ra_h) % 24.0
    ha = np.where(ha > 12.0, ha - 24.0, ha)
    return ha


# ---------------------------------------------------------------------------
# JPL Horizons ids for the solar-system minor bodies (keyed by CommonName).
# ---------------------------------------------------------------------------
MINOR_BODIES = {
    "Titan":      ("606", "majorbody"),
    "Pallas":     ("2",   "smallbody"),
    "Vesta":      ("4",   "smallbody"),
    "Juno":       ("3",   "smallbody"),
    "Flora":      ("8",   "smallbody"),
    "Melpomene":  ("18",  "smallbody"),
    "Bamberga":   ("324", "smallbody"),
    # Comet: a bare "10P" is ambiguous across apparitions; DES=..;CAP selects the
    # closest apparition. id_type must be None for this designation form.
    "Tempel 2":   ("DES=10P;CAP", None),
}

MAJOR_PLANETS = {"Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"}

CA_BUNDLE_PATH = "ca-bundle.pem"       # cached trust bundle for JPL Horizons
HORIZONS_HOST = "ssd.jpl.nasa.gov"


# ---------------------------------------------------------------------------
# SSL trust bundle for JPL Horizons.
# ssd.jpl.nasa.gov serves a valid Entrust/Sectigo chain but omits the
# intermediate CAs, and `requests` (unlike browsers/Windows) does not follow
# the AIA "caIssuers" pointers to fetch them. We fetch the missing chain over
# plain HTTP (standard AIA behaviour) and append it to certifi's roots, so TLS
# verification stays fully ON.  Cached to CA_BUNDLE_PATH after the first build.
# ---------------------------------------------------------------------------
def ensure_horizons_ca_bundle():
    if os.path.exists(CA_BUNDLE_PATH):
        os.environ["REQUESTS_CA_BUNDLE"] = os.path.abspath(CA_BUNDLE_PATH)
        return
    import certifi
    from cryptography import x509
    from cryptography.hazmat.primitives.serialization import Encoding, pkcs7

    ctx = ssl._create_unverified_context()
    with socket.create_connection((HORIZONS_HOST, 443), timeout=20) as s:
        with ctx.wrap_socket(s, server_hostname=HORIZONS_HOST) as ss:
            leaf = x509.load_der_x509_certificate(ss.getpeercert(binary_form=True))

    def load_any(data):
        for loader in (x509.load_der_x509_certificate, x509.load_pem_x509_certificate):
            try:
                return [loader(data)]
            except Exception:
                pass
        for loader in (pkcs7.load_der_pkcs7_certificates, pkcs7.load_pem_pkcs7_certificates):
            try:
                return list(loader(data))
            except Exception:
                pass
        return []

    def ca_issuers(cert):
        try:
            aia = cert.extensions.get_extension_for_class(
                x509.AuthorityInformationAccess).value
            return [d.access_location.value for d in aia
                    if d.access_method ==
                    x509.oid.AuthorityInformationAccessOID.CA_ISSUERS]
        except Exception:
            return []

    pems, seen_certs, seen_urls = [], set(), set()
    queue = ca_issuers(leaf)
    while queue:
        url = queue.pop(0)
        if url in seen_urls:
            continue
        seen_urls.add(url)
        try:
            data = urllib.request.urlopen(url, timeout=20).read()
        except Exception:
            continue
        for cert in load_any(data):
            fp = cert.fingerprint(cert.signature_hash_algorithm)
            if fp in seen_certs:
                continue
            seen_certs.add(fp)
            pems.append(cert.public_bytes(Encoding.PEM).decode())
            if cert.subject != cert.issuer:      # keep climbing until self-signed root
                queue.extend(ca_issuers(cert))

    with open(certifi.where()) as fh:
        base = fh.read()
    with open(CA_BUNDLE_PATH, "w") as fh:
        fh.write(base + "\n" + "\n".join(pems))
    os.environ["REQUESTS_CA_BUNDLE"] = os.path.abspath(CA_BUNDLE_PATH)


# ---------------------------------------------------------------------------
# Parse the skylist (capturing raw blocks so they can be re-emitted verbatim)
# ---------------------------------------------------------------------------
def parse_skylist(path):
    """Return (header_lines, objects).

    header_lines : list[str]  -- lines before the first SkyObject block (no newlines).
    objects      : list[dict] with keys type, common_names, catalog_numbers, index,
                   and raw_lines (the block's lines incl. Begin/End, no newlines).
    """
    header_lines = []
    objects = []
    cur = None
    seen_first_block = False
    with open(path, "r", encoding="utf-8", newline="") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            stripped = line.strip()
            if stripped == "SkyObject=BeginObject":
                seen_first_block = True
                cur = {"type": None, "common_names": [], "catalog_numbers": [],
                       "index": None, "raw_lines": [line]}
            elif stripped == "EndObject=SkyObject":
                if cur is not None:
                    cur["raw_lines"].append(line)
                    objects.append(cur)
                cur = None
            elif cur is not None:
                cur["raw_lines"].append(line)
                if "=" in stripped:
                    key, val = stripped.split("=", 1)
                    if key == "ObjectID":
                        cur["type"] = int(val.split(",", 1)[0])
                    elif key == "CommonName":
                        cur["common_names"].append(val)
                    elif key == "CatalogNumber":
                        cur["catalog_numbers"].append(val)
                    elif key == "DefaultIndex":
                        cur["index"] = int(val)
            elif not seen_first_block:
                header_lines.append(line)
    return header_lines, objects


# ---------------------------------------------------------------------------
# Display names (decode SkySafari's "$<x>" Bayer-Greek encoding)
# ---------------------------------------------------------------------------
GREEK = {"a": "α", "b": "β", "g": "γ", "d": "δ", "e": "ε",
         "z": "ζ", "h": "η", "q": "θ", "i": "ι", "k": "κ",
         "l": "λ", "m": "μ", "n": "ν", "c": "ξ", "o": "ο",
         "p": "π", "r": "ρ", "s": "σ", "t": "τ", "u": "υ",
         "f": "φ", "x": "χ", "y": "ψ", "w": "ω"}


def decode_bayer(name):
    if name.startswith("$") and len(name) > 1:
        letter = GREEK.get(name[1])
        if letter:
            return letter + name[2:]
    return name


# Plain-ASCII Bayer abbreviations (for output that must avoid Unicode).
GREEK_ABBREV = {"a": "alf", "b": "bet", "g": "gam", "d": "del", "e": "eps",
                "z": "zet", "h": "eta", "q": "the", "i": "iot", "k": "kap",
                "l": "lam", "m": "mu", "n": "nu", "c": "xi", "o": "omi",
                "p": "pi", "r": "rho", "s": "sig", "t": "tau", "u": "ups",
                "f": "phi", "x": "chi", "y": "psi", "w": "ome"}


def decode_bayer_ascii(name):
    if name.startswith("$") and len(name) > 1:
        abbr = GREEK_ABBREV.get(name[1])
        if abbr:
            return abbr + name[2:]
    return name


def display_name(obj):
    if obj["common_names"]:
        return obj["common_names"][0]
    if obj["catalog_numbers"]:
        return decode_bayer(obj["catalog_numbers"][0])
    return "(unnamed)"


# ---------------------------------------------------------------------------
# Coordinate resolution
# ---------------------------------------------------------------------------
# Preference order of catalog prefixes for name resolution (most reliable first).
CATALOG_PREFERENCE = ["M ", "NGC ", "IC ", "HIP ", "HD ", "HR ", "SAO ", "Mel ", "Cr ", "C "]

# Some abbreviations SkySafari uses aren't understood by the CDS name resolver;
# expand them to the forms Sesame accepts.
PREFIX_EXPANSION = {"Cr ": "Collinder ", "Mel ": "Melotte "}


def compass(az_deg):
    COMPASS_16 = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                  "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return COMPASS_16[int((az_deg % 360) / 22.5 + 0.5) % 16]


def candidate_identifiers(obj):
    """Ordered, de-duplicated list of names to try against the resolver."""
    cats = obj["catalog_numbers"]
    ordered = []
    for pref in CATALOG_PREFERENCE:
        for c in cats:
            if c.startswith(pref):
                ordered.append(c)
                if pref in PREFIX_EXPANSION:
                    ordered.append(c.replace(pref, PREFIX_EXPANSION[pref], 1))
    ordered.extend(obj["common_names"])          # e.g. "Hyades", "Coathanger"
    for c in cats:                                # anything left (skip Bayer "$..")
        if not c.startswith("$"):
            ordered.append(c)
    seen, out = set(), []
    for name in ordered:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def _skycoord_for(obj):
    """Resolve a star/DSO to a fixed SkyCoord; returns (coord, resolved_via)."""
    last_err = None
    for ident in candidate_identifiers(obj):
        try:
            return SkyCoord.from_name(ident), ident
        except Exception as e:
            last_err = e
    raise last_err if last_err else ValueError("no identifiers to resolve")


def _minor_body_name(obj):
    for n in obj["common_names"]:
        if n in MINOR_BODIES:
            return n
    return None


def _major_planet_name(obj):
    for n in obj["common_names"]:
        if n in MAJOR_PLANETS:
            return n
    return None


def _as_float_array(quantity):
    """Return a 1-D numpy array of degrees for a scalar-or-array Angle/quantity."""
    return np.atleast_1d(np.asarray(quantity.to_value(u.deg), dtype=float))


def object_altaz(obj, times, location):
    """Compute horizon coordinates for one object at one or many times.

    `times` is an astropy Time (scalar or array). Returns a dict with keys
    alt, az, ra, dec (each a 1-D numpy array in degrees, aligned to `times`)
    and `via` (the identifier / body name used).
    """
    t = obj["type"]

    if t == 1:  # solar-system body
        minor = _minor_body_name(obj)
        if minor is not None:
            ensure_horizons_ca_bundle()
            from astroquery.jplhorizons import Horizons
            hid, id_type = MINOR_BODIES[minor]
            jd = list(np.atleast_1d(times.jd))
            hz = Horizons(id=hid, id_type=id_type,
                          location={"lon": LON_DEG, "lat": LAT_DEG,
                                    "elevation": HEIGHT_M / 1000.0},
                          epochs=jd)
            eph = hz.ephemerides(quantities="4,1")  # 4 = az/el, 1 = astrometric RA/Dec
            return {
                "alt": np.asarray(eph["EL"], dtype=float),
                "az": np.asarray(eph["AZ"], dtype=float),
                "ra": np.asarray(eph["RA"], dtype=float),
                "dec": np.asarray(eph["DEC"], dtype=float),
                "via": minor,
            }
        planet = _major_planet_name(obj)
        if planet is not None:
            body = get_body(planet.lower(), times, location)
            aa = body.transform_to(AltAz(obstime=times, location=location))
            icrs = body.icrs
            return {"alt": _as_float_array(aa.alt), "az": _as_float_array(aa.az),
                    "ra": _as_float_array(icrs.ra), "dec": _as_float_array(icrs.dec),
                    "via": planet}
        raise ValueError("unknown solar-system body: %s" % obj["common_names"])

    # star (2) or deep-sky (4): fixed coordinate, transform at each time
    coord, via = _skycoord_for(obj)
    aa = coord.transform_to(AltAz(obstime=times, location=location))
    icrs = coord.icrs
    n = np.atleast_1d(times.jd).size
    return {
        "alt": _as_float_array(aa.alt),
        "az": _as_float_array(aa.az),
        "ra": np.full(n, float(icrs.ra.to_value(u.deg))),
        "dec": np.full(n, float(icrs.dec.to_value(u.deg))),
        "via": via,
    }


# ---------------------------------------------------------------------------
# Type / magnitude lookups (used for the Objects4-style CSV export)
# ---------------------------------------------------------------------------
_SIMBAD = None


def _get_simbad():
    global _SIMBAD
    if _SIMBAD is None:
        ensure_horizons_ca_bundle()          # shared trust bundle also covers CDS/SIMBAD
        from astroquery.simbad import Simbad
        s = Simbad()
        for field in ("otype", "V", "morphtype"):
            try:
                s.add_votable_fields(field)
            except Exception:
                pass
        _SIMBAD = s
    return _SIMBAD


def simbad_type_mag(obj):
    """Return (otype_code, morph_type, v_mag) for a star/DSO, trying each identifier.

    Any element may be None if SIMBAD has no value / the object doesn't resolve.
    """
    s = _get_simbad()
    for ident in candidate_identifiers(obj):
        try:
            r = s.query_object(ident)
        except Exception:
            r = None
        if r is None or len(r) == 0:
            continue

        def cell(*names):
            for col in names:
                if col in r.colnames:
                    val = r[0][col]
                    if val is None:
                        continue
                    try:
                        if np.ma.is_masked(val):
                            continue
                    except Exception:
                        pass
                    return val
            return None

        otype = cell("otype", "OTYPE")
        morph = cell("morph_type", "morphtype", "MORPH_TYPE")
        vmag = cell("V", "FLUX_V")
        otype = str(otype).strip() if otype is not None else None
        morph = str(morph).strip() if morph is not None else None
        try:
            vmag = float(vmag) if vmag is not None else None
        except Exception:
            vmag = None
        return otype, morph, vmag
    return None, None, None


# ---------------------------------------------------------------------------
# OpenNGC: authoritative Type / V-magnitude / morphology for NGC/IC/Messier
# objects. Downloaded once and cached to OPENNGC_PATH.
# ---------------------------------------------------------------------------
OPENNGC_PATH = "OpenNGC.csv"
OPENNGC_URL = ("https://raw.githubusercontent.com/mattiaverga/OpenNGC/"
               "master/database_files/NGC.csv")
_OPENNGC = None  # indices: by_name, by_m, by_ngc, by_ic


def _ensure_opengc_file():
    if not os.path.exists(OPENNGC_PATH):
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
        data = urllib.request.urlopen(OPENNGC_URL, context=ctx, timeout=60).read()
        with open(OPENNGC_PATH, "wb") as fh:
            fh.write(data)


def _load_opengc():
    global _OPENNGC
    if _OPENNGC is not None:
        return _OPENNGC
    _ensure_opengc_file()
    import csv
    by_name, by_m, by_ngc, by_ic = {}, {}, {}, {}
    with open(OPENNGC_PATH, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            by_name[row["Name"].strip()] = row
            for col, idx in (("M", by_m), ("NGC", by_ngc), ("IC", by_ic)):
                val = (row.get(col) or "").strip()
                if val.isdigit():
                    idx.setdefault(int(val), row)
    _OPENNGC = (by_name, by_m, by_ngc, by_ic)
    return _OPENNGC


def _catalog_number(cats, prefix):
    """Return the integer N from a 'PREFIX N' catalog id (e.g. 'NGC 6960' -> 6960)."""
    for c in cats:
        if c.startswith(prefix):
            digits = "".join(ch for ch in c[len(prefix):] if ch.isdigit())
            if digits:
                return int(digits)
    return None


def opengc_row(obj):
    """Return the OpenNGC row (dict) for a deep-sky object, or None."""
    by_name, by_m, by_ngc, by_ic = _load_opengc()
    cats = obj["catalog_numbers"]
    ngc = _catalog_number(cats, "NGC ")
    ic = _catalog_number(cats, "IC ")
    m = _catalog_number(cats, "M ")
    if ngc is not None and ("NGC%04d" % ngc) in by_name:
        return by_name["NGC%04d" % ngc]
    if ic is not None and ("IC%04d" % ic) in by_name:
        return by_name["IC%04d" % ic]
    if ngc is not None and ngc in by_ngc:
        return by_ngc[ngc]
    if ic is not None and ic in by_ic:
        return by_ic[ic]
    if m is not None and m in by_m:
        return by_m[m]
    return None


# Horizons body id for the major planets (planet centre, avoids barycentre ambiguity).
PLANET_HORIZONS_ID = {"Mercury": "199", "Venus": "299", "Mars": "499",
                      "Jupiter": "599", "Saturn": "699", "Uranus": "799",
                      "Neptune": "899"}


def solar_type_mag(obj, time):
    """Return (type_word, magnitude_or_None) for a solar-system (type-1) object."""
    minor = _minor_body_name(obj)
    planet = _major_planet_name(obj)
    if minor is not None:
        hid, id_type = MINOR_BODIES[minor]
        if hid.startswith("DES=") or "P" in hid and hid[0].isdigit():
            type_word = "Comet"
        elif id_type == "majorbody":
            type_word = "Moon"
        else:
            type_word = "Asteroid"
    elif planet is not None:
        hid, id_type = PLANET_HORIZONS_ID[planet], "majorbody"
        type_word = "Planet"
    else:
        return "Solar System", None

    ensure_horizons_ca_bundle()
    from astroquery.jplhorizons import Horizons
    try:
        hz = Horizons(id=hid, id_type=id_type,
                      location={"lon": LON_DEG, "lat": LAT_DEG,
                                "elevation": HEIGHT_M / 1000.0},
                      epochs=float(np.atleast_1d(time.jd)[0]))
        eph = hz.ephemerides(quantities="9")
        mag = None
        for col in ("V", "Tmag", "Nmag", "APmag"):
            if col in eph.colnames:
                val = eph[0][col]
                try:
                    if not np.ma.is_masked(val):
                        mag = float(val)
                        break
                except Exception:
                    pass
    except Exception:
        mag = None
    return type_word, mag
