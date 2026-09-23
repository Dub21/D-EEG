#!/usr/bin/env python3
"""Construit manifest.json a partir des JSON de parametres presents.

Le manifeste ne liste que les modeles dont le JSON existe, donc restreindre le
site a un sous-ensemble revient a ne generer que les JSON voulus.

Variables d'environnement :
  OUT    dossier des JSON            (defaut ~/mnt/D-EEG/static/data/normative)
  STATS  model_stats.csv, optionnel  (diagnostics affiches sous le graphe)
"""
import os, re, json, csv, sys

OUT   = os.path.expanduser(os.environ.get("OUT", "~/mnt/D-EEG/static/data/normative"))
STATS = os.path.expanduser(os.environ["STATS"]) if os.environ.get("STATS") else None
ONLY  = re.compile(os.environ["ONLY"]) if os.environ.get("ONLY") else None

BANDS    = ["delta", "theta", "alpha", "beta", "gamma", "all"]
# types sans bande : la bande fait partie du nom du parametre lui-meme
NO_BAND  = ["cf_alpha", "pw_alpha", "bw_alpha", "exponent", "offset", "knee"]


def parse(marker):
    """marker -> (scope, region, hemi, type, band)."""
    if marker.endswith("_mean"):
        stem = marker[:-len("_mean")]
        for t in NO_BAND:
            if stem == t:
                return ("global", None, None, t, None)
        for b in BANDS:
            if stem.endswith("_" + b):
                return ("global", None, None, stem[:-len(b) - 1], b)
        return ("global", None, None, stem, None)

    # regional : <region>_<hemi>_<band>_<type>
    m = re.match(r"^(.+)_(lh|rh)_(" + "|".join(BANDS) + r")_(.+)$", marker)
    if m:
        return ("regional", m.group(1), m.group(2), m.group(4), m.group(3))
    m = re.match(r"^(.+)_(lh|rh)_(.+)$", marker)
    if m:
        return ("regional", m.group(1), m.group(2), m.group(3), None)
    return (None, None, None, None, None)


def load_stats(path):
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = row.get("marker") or row.get("") or row.get("Unnamed: 0")
            if not key:
                continue
            def num(k, nd=2):
                v = row.get(k)
                try:
                    return round(float(v), nd)
                except (TypeError, ValueError):
                    return None
            out[key] = {"shapiro_p": num("Shapiro_p", 3), "smse": num("SMSE"),
                        "msll": num("MSLL"), "skew": num("Skewness"),
                        "kurt": num("Kurtosis")}
    return out


def main():
    stats = load_stats(STATS)
    models, regions, types, bands = [], set(), set(), set()

    for f in sorted(os.listdir(OUT)):
        if not f.endswith(".json") or f == "manifest.json":
            continue
        marker = f[:-5]
        if marker.endswith("_qc"):
            continue                      # variante, signalee par un drapeau
        if ONLY and not ONLY.search(marker):
            continue
        scope, region, hemi, typ, band = parse(marker)
        if scope is None:
            print("  ignore (nom non reconnu) :", marker)
            continue
        d = json.load(open(os.path.join(OUT, f)))
        entry = {"scope": scope, "region": region, "hemi": hemi,
                 "type": typ, "band": band, "marker": marker,
                 "n": d.get("n"), "conv": d.get("converged"),
                 "raw": d.get("qc_ref") is None,
                 "qc": os.path.exists(os.path.join(OUT, marker + "_qc.json"))}
        entry.update(stats.get(marker, {"shapiro_p": None, "smse": None,
                                        "msll": None, "skew": None, "kurt": None}))
        models.append(entry)
        types.add(typ)
        if band:
            bands.add(band)
        if region:
            regions.add(region)

    man = {"n_models": len(models),
           "regions": sorted(regions),
           "types": sorted(types),
           "bands": [b for b in BANDS if b in bands],
           "models": models}
    with open(os.path.join(OUT, "manifest.json"), "w") as fh:
        json.dump(man, fh, separators=(",", ":"))
    print(f"manifest : {len(models)} modeles, {len(types)} types, "
          f"{len(man['bands'])} bandes, {len(man['regions'])} regions")
    print("types :", ", ".join(man["types"]))


if __name__ == "__main__":
    main()
