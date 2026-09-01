#!/usr/bin/env python3
"""Extrait chaque modele gamlss SHASHo2 en un petit JSON de parametres.
Aucune donnee individuelle en sortie : les ages observes servent a interpoler
le terme lisse, puis sont remplaces par une grille reguliere."""
import os, sys, json, re, time, warnings, csv
warnings.filterwarnings("ignore")
import numpy as np, rdata
from multiprocessing import Pool

D    = os.path.expanduser("~/mnt/Desktop/models")
OUT  = os.path.expanduser("~/mnt/D-EEG/static/data/normative")
NAGE = 200
FORCE = os.environ.get("FORCE") == "1"
SEX, SITE_PFX, QC, AGE = "as.factor(female_bin)1", "as.factor(Unique_Site_ID)", "ratio_ch_good", "ps(age)"
EXCL = re.compile(r"(cf_alpha|pw_alpha|bw_alpha|knee)")
def keep(m): return not EXCL.search(m)

def parts(o, par):
    X = np.asarray(o[f"{par}.x"], dtype=float)
    cols = [str(c) for c in o[f"{par}.x"].coords[list(o[f"{par}.x"].coords)[-1]].values]
    coef = dict(zip(cols, np.asarray(o[f"{par}.coefficients"], dtype=float)))
    sm = o.get(f"{par}.s")
    if sm is None:
        s = np.zeros(X.shape[0]); extra = {}
    else:
        arr = np.asarray(sm, dtype=float)
        if arr.ndim == 1: arr = arr[:, None]
        # avec random(site) la matrice des lisseurs a 2 colonnes : ps(age) et le
        # terme aleatoire. On isole la colonne de l'age par son nom si dispo,
        # sinon par position (l'ordre suit la formule, ps(age) en premier).
        names = []
        try: names = [str(c) for c in sm.coords[list(sm.coords)[-1]].values]
        except Exception: pass
        iage = next((k for k, n in enumerate(names) if "age" in n), 0)
        s = arr[:, iage]
        extra = {n: arr[:, k] for k, n in enumerate(names) if k != iage}
    return X[:, cols.index(AGE)], s, coef, cols, X

def site_mean(coef, cols):
    """Moyenne des offsets de site, en incluant le site de reference (offset 0).
    Recentrer dessus rend la courbe neutre vis-a-vis du site, au lieu de la coller
    a un site arbitraire (le premier par ordre alphabetique)."""
    offs = [coef.get(c, 0.) for c in cols if c.startswith(SITE_PFX)]
    return (sum(offs) + 0.0) / (len(offs) + 1) if offs else 0.0

def eta(o, par, grid, female, qc):
    age, s, coef, cols, X = parts(o, par)
    ua, inv = np.unique(age, return_inverse=True)
    us = np.bincount(inv, weights=s) / np.bincount(inv)
    e = coef.get("(Intercept)", 0.) + coef.get(AGE, 0.)*grid + np.interp(grid, ua, us)
    if female: e = e + coef.get(SEX, 0.)
    return e + coef.get(QC, 0.)*qc + site_mean(coef, cols)

def one(f):
    marker = f[:-4]
    dst = os.path.join(OUT, marker + ".json")
    if os.path.exists(dst) and not FORCE: return (marker, "skip", None)
    try:
        o = rdata.conversion.convert(rdata.parser.parse_file(os.path.join(D, f)))
        age, _, coef, cols, X = parts(o, "mu")
        qc = float(np.median(X[:, cols.index(QC)]))
        grid = np.linspace(float(age.min()), float(age.max()), NAGE)
        out = {
            "marker": marker,
            "n": int(np.asarray(o["N"]).ravel()[0]),
            "family": str(np.asarray(o["family"]).ravel()[0]),
            "converged": bool(np.asarray(o["converged"]).ravel()[0]),
            "qc_ref": round(qc, 5),
            "nu": round(float(np.asarray(o["nu.coefficients"]).ravel()[0]), 6),
            "tau": round(float(np.exp(np.asarray(o["tau.coefficients"]).ravel()[0])), 6),
            "age": [round(float(x), 4) for x in grid],
        }
        for lab, fem in (("m", 0), ("f", 1)):
            out[lab] = {
                "mu":    [round(float(x), 5) for x in eta(o, "mu", grid, fem, qc)],
                "sigma": [round(float(x), 5) for x in np.exp(eta(o, "sigma", grid, fem, qc))],
            }
        sites = sorted({c[len(SITE_PFX):] for c in cols if c.startswith(SITE_PFX)})
        _, _, cs, ccols, _ = parts(o, "sigma")
        out["site_ref"] = "moyenne des sites"
        out["site_offsets_mu"]    = {s: round(coef.get(SITE_PFX+s, 0.) - site_mean(coef, cols), 5) for s in sites}
        out["site_offsets_sigma"] = {s: round(cs.get(SITE_PFX+s, 0.) - site_mean(cs, ccols), 5) for s in sites}
        with open(dst, "w") as fh: json.dump(out, fh, separators=(",", ":"))
        return (marker, "ok", None)
    except Exception as e:
        return (marker, "err", type(e).__name__)

if __name__ == "__main__":
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 95.
    os.makedirs(OUT, exist_ok=True)
    files = sorted(f for f in os.listdir(D) if f.endswith(".rds") and keep(f[:-4]))
    todo = files if FORCE else [f for f in files if not os.path.exists(os.path.join(OUT, f[:-4] + ".json"))]
    if not todo:
        print(f"complet : {len(files)} marqueurs"); sys.exit(0)
    t0, n, errs = time.time(), 0, []
    with Pool(3) as p:
        for m, st, e in p.imap_unordered(one, todo, chunksize=2):
            n += 1
            if st == "err": errs.append((m, e))
            if time.time() - t0 > budget: p.terminate(); break
    print(f"+{n} traites, {len(errs)} erreurs, restant ~{len(todo)-n}")
    for m, e in errs[:5]: print("  ERR", m, e)
