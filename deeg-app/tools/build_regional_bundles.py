#!/usr/bin/env python3
"""Regroupe les JSON regionaux en un bundle par famille de marqueur.

Un fichier par famille (periodic_alpha, connectivity_beta, exponent, offset...)
contenant les 68 regions dans les deux variantes, charge a la demande par la
page quand l'utilisateur selectionne cette famille. Evite un bundle unique.

Variables d'environnement :
  SRC   dossier des JSON regionaux  (defaut ~/mnt/Desktop/regional_json)
  OUT   dossier de sortie           (defaut ~/mnt/D-EEG/static/data/normative/regional)
  NAGE  points de grille            (defaut 120)
  NORM  csv marker,mean,std         (optionnel, constantes de standardisation)
"""
import os, json, csv, re, collections
import numpy as np

SRC  = os.path.expanduser(os.environ.get("SRC", "~/mnt/Desktop/regional_json"))
OUT  = os.path.expanduser(os.environ.get("OUT", "~/mnt/D-EEG/static/data/normative/regional"))
NAGE = int(os.environ.get("NAGE", "120"))
NORM = os.path.expanduser(os.environ["NORM"]) if os.environ.get("NORM") else None

BANDS = ["delta", "theta", "alpha", "beta", "gamma", "all"]
RE_B  = re.compile(r"^(.+)_(lh|rh)_(" + "|".join(BANDS) + r")_(.+)$")
RE_N  = re.compile(r"^(.+)_(lh|rh)_(.+)$")


def split(marker):
    """marker -> (famille, cle de region) ; famille = <feature>[_<band>]."""
    m = RE_B.match(marker)
    if m:
        region, hemi, band, feat = m.groups()
        return f"{feat}_{band}", f"{region}_{hemi}"
    m = RE_N.match(marker)
    if m:
        region, hemi, feat = m.groups()
        return feat, f"{region}_{hemi}"
    return None, None


def resample(d):
    age = np.asarray(d["age"], float)
    a0, a1 = float(age[0]), float(age[-1])
    grid = np.exp(np.linspace(np.log(a0), np.log(a1), NAGE))
    rec = {"a0": round(a0, 4), "a1": round(a1, 4),
           "nu": d["nu"], "tau": d["tau"], "site_sd": d.get("site_sd", 0.0)}
    for k in ("qc_ref", "qc_coef_mu", "qc_coef_sigma"):
        if d.get(k) is not None:
            rec[k] = d[k]
    for sex in ("m", "f"):
        mu = np.interp(grid, age, np.asarray(d[sex]["mu"], float))
        sg = np.interp(grid, age, np.asarray(d[sex]["sigma"], float))
        rec[sex] = [[round(float(v), 5) for v in mu],
                    [round(float(v), 5) for v in sg]]
    return rec


def main():
    norm = {}
    if NORM and os.path.exists(NORM):
        for r in csv.DictReader(open(NORM)):
            try:
                norm[r["marker"]] = (float(r["mean"]), float(r["std"]))
            except (TypeError, ValueError):
                pass

    os.makedirs(OUT, exist_ok=True)
    fams = collections.defaultdict(dict)
    regions, n, skipped = set(), 0, []
    for f in sorted(os.listdir(SRC)):
        if not f.endswith(".json"):
            continue
        stem = f[:-5]
        suffix = "_qc" if stem.endswith("_qc") else ""
        base = stem[:-3] if suffix else stem
        fam, key = split(base)
        if not fam:
            print("  ignore :", base); continue
        path = os.path.join(SRC, f)
        try:
            d = json.load(open(path))
        except (ValueError, OSError):
            skipped.append(f)          # fichier vide ou tronque
            continue
        rec = resample(d)
        if base in norm:
            rec["mean"], rec["std"] = norm[base]
        fams[fam][key + suffix] = rec
        regions.add(key); n += 1

    # une famille incomplete au point de ne couvrir qu'une poignee de regions
    # est un residu d'un run interrompu, pas une famille : on l'ecarte
    strays = [f for f, m in fams.items() if len(m) < 20]
    for f in strays:
        del fams[f]

    index, total = {}, 0
    for fam, markers in sorted(fams.items()):
        dst = os.path.join(OUT, fam + ".json")
        with open(dst, "w") as fh:
            json.dump({"n_age": NAGE, "family": fam, "markers": markers},
                      fh, separators=(",", ":"))
        kb = os.path.getsize(dst) / 1024
        total += kb
        nq = sum(1 for k in markers if k.endswith("_qc"))
        index[fam] = {"n": len(markers) - nq, "n_qc": nq,
                      "kb": round(kb), "with_units": any("mean" in r for r in markers.values())}
        print(f"  {fam:24s} {len(markers)-nq:3d} + {nq:3d} qc   {kb:7.0f} Ko")

    with open(os.path.join(OUT, "index.json"), "w") as fh:
        json.dump({"n_age": NAGE, "regions": sorted(regions),
                   "families": index}, fh, separators=(",", ":"))
    print(f"\n{len(fams)} familles, {n} modeles, {len(regions)} regions, {total/1024:.1f} Mo au total")
    if strays:
        print("familles trop incompletes, ecartees :", ", ".join(sorted(strays)))
    if skipped:
        print(f"{len(skipped)} fichier(s) illisible(s), ignore(s) :", ", ".join(skipped[:5]))
    if not norm:
        print("ATTENTION : pas de constantes de standardisation, les bundles sont en unites standardisees")


if __name__ == "__main__":
    main()
