"""
Build a Markdown field-notes file describing every object in
2026-best-objects-filtered.csv, one short paragraph each, sourced from
Wikipedia's REST summary API.

Output: 2026-best-objects-filtered.md  (UTF-8 - view rendered or as UTF-8).
"""

import sys
import csv
import json
import time
import urllib.parse
import urllib.request
import ssl
import certifi

import skylib
import to_csv

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

CSV_PATH = "2026-best-objects-filtered.skylist"
CSV_OUT = "2026-best-objects-filtered.csv"
MD_OUT = "2026-best-objects-filtered.md"
UA = ("astro-observing-list/1.0 (personal astronomy project; "
      "contact antongml@gmail.com)")
SUMMARY = "https://en.wikipedia.org/api/rest_v1/page/summary/"
_CTX = ssl.create_default_context(cafile=certifi.where())

# --- Wikipedia title candidates --------------------------------------------
GREEK_FULL = {"a": "Alpha", "b": "Beta", "g": "Gamma", "d": "Delta", "e": "Epsilon",
              "z": "Zeta", "h": "Eta", "q": "Theta", "i": "Iota", "k": "Kappa",
              "l": "Lambda", "m": "Mu", "n": "Nu", "c": "Xi", "o": "Omicron",
              "p": "Pi", "r": "Rho", "s": "Sigma", "t": "Tau", "u": "Upsilon",
              "f": "Phi", "x": "Chi", "y": "Psi", "w": "Omega"}
GENITIVE = {"And": "Andromedae", "Aql": "Aquilae", "Aqr": "Aquarii", "Boo": "Bootis",
            "Cap": "Capricorni", "Cas": "Cassiopeiae", "Cep": "Cephei", "Cet": "Ceti",
            "CrB": "Coronae Borealis", "Cyg": "Cygni", "Del": "Delphini",
            "Her": "Herculis", "Lyr": "Lyrae", "Oph": "Ophiuchi", "Per": "Persei",
            "Psc": "Piscium", "Ser": "Serpentis", "Tau": "Tauri", "UMa": "Ursae Majoris",
            "UMi": "Ursae Minoris", "Vul": "Vulpeculae", "Aur": "Aurigae",
            "PsA": "Piscis Austrini", "Sge": "Sagittae"}
SOLAR_WIKI = {"Melpomene": "18 Melpomene", "Flora": "8 Flora", "Juno": "3 Juno",
              "Bamberga": "324 Bamberga", "Tempel 2": "10P/Tempel"}


def _num(cats, prefix):
    return skylib._catalog_number(cats, prefix)


def wiki_candidates(obj):
    t, cats, commons = obj["type"], obj["catalog_numbers"], obj["common_names"]
    cand = []
    if t == 1:
        name = commons[0] if commons else cats[0]
        cand.append(SOLAR_WIKI.get(name, name))
        return cand
    if t == 4:
        m, ngc, ic, cr = (_num(cats, "M "), _num(cats, "NGC "),
                          _num(cats, "IC "), _num(cats, "Cr "))
        if m:
            cand.append("Messier %d" % m)
        for cn in commons:
            cand.append(cn)
        if ngc:
            cand.append("NGC %d" % ngc)
        if ic:
            cand.append("IC %d" % ic)
        if cr:
            cand.append("Collinder %d" % cr)
        return cand
    # stars: proper names, then Bayer, then Flamsteed, then HR/HD
    for cn in commons:
        cand.append(cn)
    for c in cats:
        if c.startswith("$"):
            letter = GREEK_FULL.get(c[1])
            parts = c[2:].strip().split()
            num, con = "", None
            if parts and parts[0].isdigit():
                num, con = parts[0], (parts[1] if len(parts) > 1 else None)
            elif parts:
                con = parts[0]
            gen = GENITIVE.get(con)
            if letter and gen:
                if num:
                    cand.append("%s%s %s" % (letter, num, gen))
                cand.append("%s %s" % (letter, gen))
    for c in cats:
        p = c.split()
        if len(p) == 2 and p[0].isdigit() and p[1] in GENITIVE:
            cand.append("%s %s" % (p[0], GENITIVE[p[1]]))
    for c in cats:
        if c.startswith(("HR ", "HD ", "HIP ")):
            cand.append(c)
    return cand


def fetch_summary(title):
    url = SUMMARY + urllib.parse.quote(title, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        raw = urllib.request.urlopen(req, context=_CTX, timeout=25).read()
    except Exception:
        return None
    d = json.loads(raw)
    if d.get("type") == "disambiguation":
        return None
    extract = (d.get("extract") or "").strip()
    if not extract:
        return None
    page = d.get("content_urls", {}).get("desktop", {}).get("page")
    return {"extract": extract, "url": page, "title": d.get("title")}


def describe(obj):
    for title in wiki_candidates(obj):
        res = fetch_summary(title)
        time.sleep(0.05)
        if res:
            return res
    return None


def heading(obj):
    label = to_csv.primary_label(obj)                     # ASCII catalog/star id
    nickname = obj["common_names"][0] if (obj["type"] != 1 and obj["common_names"]) else None
    if nickname and nickname != label:
        return "%s (%s)" % (nickname, label)
    return label


def main():
    _, objects = skylib.parse_skylist(CSV_PATH)

    # magnitude + coordinates straight from the generated CSV (same order)
    with open(CSV_OUT, "r", encoding="utf-8", newline="") as fh:
        rows = [r for r in csv.reader(fh) if r]
    data_rows = rows[2:]                                   # skip title + schema

    out = []
    out.append("# 2026 Best Objects - Field Notes")
    out.append("")
    out.append("Observing list for Adin, CA, night of Sat 11 Jul 2026 "
               "(filtered to objects above 10 deg between 10 pm and midnight PDT). "
               "Descriptions from Wikipedia.")
    out.append("")

    missing = []
    for i, (obj, row) in enumerate(zip(objects, data_rows), start=1):
        n_field = row[0]
        mag = n_field.rsplit(" ", 1)[-1]
        rh, rm, dd, dm = row[1], row[2], row[3], row[4]
        res = describe(obj)
        cats = " / ".join(obj["catalog_numbers"][:4])

        out.append("## %d. %s" % (i, heading(obj)))
        meta = "*%s | mag %s | RA %sh %sm | Dec %s° %s′*" % (cats, mag, rh, rm, dd, dm)
        out.append(meta)
        out.append("")
        if res:
            out.append(res["extract"])
            if res.get("url"):
                out.append("")
                out.append("[Wikipedia](%s)" % res["url"])
        else:
            out.append("_No description found._")
            missing.append(heading(obj))
        out.append("")

    with open(MD_OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))

    print("Wrote %s  (%d objects)" % (MD_OUT, len(data_rows)))
    if missing:
        print("No Wikipedia summary for:", ", ".join(missing))


if __name__ == "__main__":
    main()
