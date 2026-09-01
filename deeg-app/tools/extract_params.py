#!/usr/bin/env python3
"""
Extrait les parametres normatifs (mu, sigma, nu, tau) des objets gamlss SHASHo2
et les evalue sur une grille d'age reguliere.

Tourne en local, sans R. Ne sort AUCUNE donnee individuelle : les ages observes
des participants servent uniquement a interpoler le terme lisse, puis sont jetes.

Structure du modele (lue dans les objets) :
    mu    ~ ps(age) + as.factor(female_bin) + as.factor(Unique_Site_ID) + ratio_ch_good   [lien identite]
    sigma ~ idem                                                                          [lien log]
    nu    ~ 1  [identite]
    tau   ~ 1  [log]

Reconstruction :
    eta = b0 + b_age*age + s(age) + b_sex*I(female) + b_site + b_qc*qc
ou s(age) est le terme lisse P-spline, interpole depuis (age_i, s_i).

Usage :
    python3 extract_params.py <MODEL_DIR> <OUT_DIR> [--n-age 300] [--check]
"""
import argparse, os, sys, warnings
import numpy as np
import rdata
import re

warnings.filterwarnings("ignore")

SEX_COL  = "as.factor(female_bin)1"
SITE_PFX = "as.factor(Unique_Site_ID)"
QC_COL   = "ratio_ch_good"
AGE_COL  = "ps(age)"


# --- selection des marqueurs -------------------------------------------------
# cf_alpha, pw_alpha, bw_alpha et les variantes "knee" ne sont valides qu'au
# niveau global (_mean) : leurs versions regionales sont ecartees.
REGIONAL_EXCLUDE = re.compile(r"(cf_alpha|pw_alpha|bw_alpha|knee)")


def is_regional(marker):
    return "_lh_" in marker or "_rh_" in marker


def keep_marker(marker):
    if not is_regional(marker):
        return True                      # tous les marqueurs globaux *_mean
    return not REGIONAL_EXCLUDE.search(marker)


def load(path):
    return rdata.conversion.convert(rdata.parser.parse_file(path))


def cols_of(o, par):
    x = o[f"{par}.x"]
    return [str(c) for c in x.coords[list(x.coords)[-1]].values]


def parts(o, par):
    """Renvoie (age, s_age, coefs dict, design cols, X)."""
    X = np.asarray(o[f"{par}.x"], dtype=float)
    cols = cols_of(o, par)
    coef = np.asarray(o[f"{par}.coefficients"], dtype=float)
    s = o.get(f"{par}.s")
    s = np.asarray(s, dtype=float).ravel() if s is not None else np.zeros(X.shape[0])
    age = X[:, cols.index(AGE_COL)]
    return age, s, dict(zip(cols, coef)), cols, X


def smooth_interp(age, s, grid):
    """Interpole le terme lisse sur la grille. s ne depend que de l'age."""
    order = np.argsort(age)
    a, sv = age[order], s[order]
    # moyenne des doublons d'age
    ua, inv = np.unique(a, return_inverse=True)
    us = np.bincount(inv, weights=sv) / np.bincount(inv)
    return np.interp(grid, ua, us)  # extrapolation plate hors plage (grille bornee de toute facon)


def eta_on_grid(o, par, grid, sex, site, qc):
    age, s, coef, cols, X = parts(o, par)
    b = coef.get("(Intercept)", 0.0)
    e = b + coef.get(AGE_COL, 0.0) * grid + smooth_interp(age, s, grid)
    if sex == 1:
        e = e + coef.get(SEX_COL, 0.0)
    if site is not None:
        e = e + coef.get(SITE_PFX + site, 0.0)  # 0.0 = site de reference
    e = e + coef.get(QC_COL, 0.0) * qc
    return e


def model_params(o, grid, sex, site, qc):
    mu = eta_on_grid(o, "mu", grid, sex, site, qc)                 # lien identite
    sigma = np.exp(eta_on_grid(o, "sigma", grid, sex, site, qc))   # lien log
    nu = float(np.asarray(o["nu.coefficients"]).ravel()[0])        # identite, ~1
    tau = float(np.exp(np.asarray(o["tau.coefficients"]).ravel()[0]))  # log, ~1
    return mu, sigma, nu, tau


def check(o):
    """Compare la reconstruction aux valeurs ajustees stockees. Doit etre ~0."""
    out = {}
    for par, link in (("mu", "identity"), ("sigma", "log")):
        age, s, coef, cols, X = parts(o, par)
        c = np.array([coef[k] for k in cols])
        eta = X @ c + s
        fv = np.asarray(o[f"{par}.fv"], dtype=float).ravel()
        pred = eta if link == "identity" else np.exp(eta)
        out[par] = float(np.max(np.abs(pred - fv)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model_dir"); ap.add_argument("out_dir")
    ap.add_argument("--n-age", type=int, default=300)
    ap.add_argument("--check", action="store_true", help="valide la reconstruction et s'arrete")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    files = sorted(f for f in os.listdir(a.model_dir) if f.lower().endswith(".rds"))
    n_all = len(files)
    files = [f for f in files if keep_marker(f[:-4])]
    print(f"{n_all} modeles, {n_all - len(files)} ecartes (cf/pw/bw_alpha et knee regionaux), {len(files)} retenus")
    if a.limit: files = files[:a.limit]
    os.makedirs(a.out_dir, exist_ok=True)

    if a.check:
        worst = {"mu": 0.0, "sigma": 0.0}
        for f in files:
            e = check(load(os.path.join(a.model_dir, f)))
            for k in worst: worst[k] = max(worst[k], e[k])
            print(f"{f:55s} mu {e['mu']:.2e}  sigma {e['sigma']:.2e}")
        print(f"\nerreur max : mu {worst['mu']:.3e}  sigma {worst['sigma']:.3e}")
        return

    rows, manifest = [], []
    for i, f in enumerate(files, 1):
        marker = f[:-4]
        try:
            o = load(os.path.join(a.model_dir, f))
        except Exception as e:
            print(f"[SKIP] {marker}: {e}"); continue

        age, _, coef, cols, _ = parts(o, "mu")
        grid = np.linspace(age.min(), age.max(), a.n_age)
        sites = [c[len(SITE_PFX):] for c in cols if c.startswith(SITE_PFX)]
        qc_ref = float(np.median(np.asarray(o["mu.x"], dtype=float)[:, cols.index(QC_COL)]))

        for sex in (0, 1):
            mu, sigma, nu, tau = model_params(o, grid, sex, None, qc_ref)
            for g, m, s in zip(grid, mu, sigma):
                rows.append((marker, float(g), sex, float(m), float(s), nu, tau))

        manifest.append(dict(
            marker=marker,
            family=str(np.asarray(o["family"]).ravel()[0]),
            n_train=int(np.asarray(o["N"]).ravel()[0]),
            age_min=float(age.min()), age_max=float(age.max()),
            qc_ref=qc_ref,
            site_ref=sorted(set(sites))[0] if sites else "",
            n_sites=len(sites) + 1,
            aic=float(np.asarray(o["aic"]).ravel()[0]),
            site_offsets="|".join(f"{s}={coef.get(SITE_PFX+s, 0.0):.6g}" for s in sites),
        ))
        if i % 50 == 0: print(f"  {i}/{len(files)}")

    import csv, gzip
    with gzip.open(os.path.join(a.out_dir, "normative_params.csv.gz"), "wt", newline="") as fh:
        w = csv.writer(fh); w.writerow(["marker","age","female","mu","sigma","nu","tau"])
        w.writerows(rows)
    with open(os.path.join(a.out_dir, "model_manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys())); w.writeheader(); w.writerows(manifest)
    print(f"\n{len(manifest)} modeles -> {a.out_dir}")


if __name__ == "__main__":
    main()
