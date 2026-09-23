#!/usr/bin/env python3
"""Extrait chaque modele gamlss SHASHo2 en un petit JSON de parametres.

Aucune donnee individuelle en sortie : les valeurs observees du predicteur
d'age servent a interpoler le terme lisse, puis sont remplacees par une
grille reguliere.

Les noms de termes sont detectes dans la matrice de design, donc le script
suit la specification du modele au lieu de la supposer : ps(age), pb(age),
pb(log_age), avec ou sans covariable de qualite.

Variables d'environnement :
  MODELS   dossier des .rds            (defaut ~/mnt/Desktop/models)
  OUT      dossier de sortie           (defaut ~/mnt/D-EEG/static/data/normative)
  ONLY     regex de selection          (ex '_mean$' pour les 20 marqueurs globaux)
  FORCE=1  reecrit les JSON existants
"""
import os, sys, json, re, time, warnings
warnings.filterwarnings("ignore")
import numpy as np, rdata
from multiprocessing import Pool

D    = os.path.expanduser(os.environ.get("MODELS", "~/mnt/Desktop/models"))
OUT  = os.path.expanduser(os.environ.get("OUT", "~/mnt/D-EEG/static/data/normative"))
NAGE = 200
FORCE = os.environ.get("FORCE") == "1"
ONLY  = re.compile(os.environ["ONLY"]) if os.environ.get("ONLY") else None
SUFFIX = os.environ.get("SUFFIX", "")   # ex "_qc" pour la version ajustee sur la qualite

SITE_KEY = "Unique_Site_ID"
SEX_KEY  = "female_bin"
QC_KEY   = "ratio_ch_good"
# marqueurs de pic alpha et de knee : valides seulement en global, pas par region
EXCL = re.compile(r"(cf_alpha|pw_alpha|bw_alpha|knee)")


def keep(m):
    if ONLY and not ONLY.search(m):
        return False
    return m.endswith("_mean") or not EXCL.search(m)


def find_col(cols, key, exclude=()):
    for c in cols:
        if key in c and not any(e in c for e in exclude):
            return c
    return None


def design(o, par):
    """Colonnes, coefficients et matrice de design d'un parametre."""
    X = np.asarray(o[f"{par}.x"], dtype=float)
    last = list(o[f"{par}.x"].coords)[-1]
    cols = [str(c) for c in o[f"{par}.x"].coords[last].values]
    coef = dict(zip(cols, np.asarray(o[f"{par}.coefficients"], dtype=float)))
    return X, cols, coef


def smooth_of(o, par, cols):
    """Colonne lissee liee a l'age, sommee sur les sujets. Les autres colonnes
    de lisseurs (effet aleatoire) sont ignorees."""
    sm = o.get(f"{par}.s")
    if sm is None:
        return None
    arr = np.asarray(sm, dtype=float)
    if arr.ndim == 1:
        arr = arr[:, None]
    names = []
    try:
        names = [str(c) for c in sm.coords[list(sm.coords)[-1]].values]
    except Exception:
        pass
    iage = next((k for k, n in enumerate(names) if "age" in n), 0)
    return arr[:, iage]


def terms(o, par):
    """Renvoie tout ce qu'il faut pour reconstruire le predicteur lineaire."""
    X, cols, coef = design(o, par)
    age_col  = find_col(cols, "age")
    sex_col  = find_col(cols, SEX_KEY)
    qc_col   = find_col(cols, QC_KEY)
    site_cols = [c for c in cols if SITE_KEY in c
                 and not c.startswith(('re(', 'random('))]
    x = X[:, cols.index(age_col)] if age_col else np.zeros(X.shape[0])
    return dict(X=X, cols=cols, coef=coef, x=x, age_col=age_col, sex_col=sex_col,
                qc_col=qc_col, site_cols=site_cols, s=smooth_of(o, par, cols))


def site_mean(coef, site_cols):
    """Moyenne des offsets de site en incluant le niveau de reference (0), pour
    que la courbe soit neutre vis-a-vis du site plutot que collee au premier."""
    if not site_cols:
        return 0.0
    offs = [cf(coef, c) for c in site_cols]
    return (sum(offs) + 0.0) / (len(offs) + 1)


def cf(coef, key, default=0.0):
    v = coef.get(key, default)
    return default if v is None or not np.isfinite(v) else float(v)



def blup(o, par):
    """Offsets de site quand le terme est un effet aleatoire random() : les BLUP
    vivent dans coefSmo, pas dans les coefficients parametriques. Renvoie un
    dict niveau -> offset, recentre sur la moyenne des niveaux."""
    cs = o.get(f"{par}.coefSmo")
    if cs is None:
        return {}
    try:
        items = list(cs)
    except TypeError:
        return {}
    for el in items:
        if not hasattr(el, "keys") or "factor" not in el or "coef" not in el:
            continue
        c = el["coef"]
        try:
            names = [str(x) for x in c.coords[list(c.coords)[-1]].values]
        except Exception:
            continue
        vals = np.asarray(c, dtype=float).ravel()
        if len(names) != len(vals):
            continue
        m = float(np.nanmean(vals))
        return {n: float(v - m) for n, v in zip(names, vals)}
    return {}


def eta(t, grid_x, female, qc):
    coef = t["coef"]
    e = np.full(len(grid_x), cf(coef, "(Intercept)"))
    if t["age_col"]:
        e = e + cf(coef, t["age_col"]) * grid_x
    if t["s"] is not None:
        ux, inv = np.unique(t["x"], return_inverse=True)
        us = np.bincount(inv, weights=t["s"]) / np.bincount(inv)
        e = e + np.interp(grid_x, ux, us)
    if female and t["sex_col"]:
        e = e + cf(coef, t["sex_col"])
    if t["qc_col"] and qc is not None:
        e = e + cf(coef, t["qc_col"]) * qc
    return e + site_mean(coef, t["site_cols"])


def one(f):
    marker = f[:-4]
    dst = os.path.join(OUT, marker + SUFFIX + ".json")
    if os.path.exists(dst) and not FORCE:
        return (marker, "skip", None)
    try:
        o  = rdata.conversion.convert(rdata.parser.parse_file(os.path.join(D, f)))
        tm = terms(o, "mu")
        ts = terms(o, "sigma")

        qc = None
        if tm["qc_col"]:
            qc = float(np.median(tm["X"][:, tm["cols"].index(tm["qc_col"])]))

        is_log = bool(tm["age_col"]) and "log" in tm["age_col"]
        grid_x = np.linspace(float(tm["x"].min()), float(tm["x"].max()), NAGE)
        ages   = np.exp(grid_x) if is_log else grid_x

        out = {
            "marker": marker,
            "n": int(np.asarray(o["N"]).ravel()[0]),
            "family": str(np.asarray(o["family"]).ravel()[0]),
            "converged": bool(np.asarray(o["converged"]).ravel()[0]),
            "age_term": tm["age_col"] or "",
            "qc_ref": (round(qc, 5) if qc is not None else None),
            "qc_coef_mu":    (round(cf(tm["coef"], tm["qc_col"]), 6) if tm["qc_col"] else None),
            "qc_coef_sigma": (round(cf(ts["coef"], ts["qc_col"]), 6) if ts["qc_col"] else None),
            "nu":  round(float(np.asarray(o["nu.coefficients"]).ravel()[0]), 6),
            "tau": round(float(np.exp(np.asarray(o["tau.coefficients"]).ravel()[0])), 6),
            "age": [round(float(x), 4) for x in ages],
        }
        for lab, fem in (("m", 0), ("f", 1)):
            out[lab] = {
                "mu":    [round(float(x), 5) for x in eta(tm, grid_x, fem, qc)],
                "sigma": [round(float(x), 5) for x in np.exp(eta(ts, grid_x, fem, qc))],
            }

        out["site_ref"] = "moyenne des sites"
        for par, t, lab in (("mu", tm, "site_offsets_mu"), ("sigma", ts, "site_offsets_sigma")):
            off = blup(o, par)
            if not off and t["site_cols"]:
                pfx = os.path.commonprefix(t["site_cols"])
                m0  = site_mean(t["coef"], t["site_cols"])
                off = {c[len(pfx):]: cf(t["coef"], c) - m0 for c in t["site_cols"]}
            out[lab] = {k: round(v, 5) for k, v in sorted(off.items())}
            if par == "mu":
                vals = list(off.values())
                # ecart-type entre sites : sert au retrecissement lors d'un
                # recalibrage sur des controles locaux
                out["site_sd"] = round(float(np.std(vals, ddof=1)), 5) if len(vals) > 1 else 0.0

        with open(dst, "w") as fh:
            json.dump(out, fh, separators=(",", ":"))
        return (marker, "ok", None)
    except Exception as e:
        return (marker, "err", f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 95.
    os.makedirs(OUT, exist_ok=True)
    files = sorted(f for f in os.listdir(D) if f.endswith(".rds") and keep(f[:-4]))
    todo  = files if FORCE else [f for f in files
                                 if not os.path.exists(os.path.join(OUT, f[:-4] + SUFFIX + ".json"))]
    print(f"{len(files)} marqueurs selectionnes, {len(todo)} a traiter")
    if not todo:
        sys.exit(0)
    t0, n, errs = time.time(), 0, []
    with Pool(3) as p:
        for m, st, e in p.imap_unordered(one, todo, chunksize=2):
            n += 1
            if st == "err":
                errs.append((m, e))
            if time.time() - t0 > budget:
                p.terminate(); break
    print(f"+{n} traites, {len(errs)} erreurs, restant ~{len(todo)-n}")
    for m, e in errs[:5]:
        print("  ERR", m, e)
