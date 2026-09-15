#!/usr/bin/env python3
"""Regroupe les JSON de parametres en un seul fichier, sur une grille d'age
logarithmique, pour le depot d'un profil complet (une requete au lieu de N).

Les JSON portent mu et sigma sur une grille lineaire propre a chaque marqueur ;
le bundle les reechantillonne sur 120 points en log entre les memes bornes, ce
qui est ce que la page attend (paramsAt interpole en log).

Variables d'environnement :
  OUT    dossier des JSON  (defaut ~/mnt/D-EEG/static/data/normative)
  NAGE   points de grille  (defaut 120)
"""
import os, json, csv
import numpy as np

OUT  = os.path.expanduser(os.environ.get("OUT", "~/mnt/D-EEG/static/data/normative"))
NAGE = int(os.environ.get("NAGE", "120"))


def main():
    man = json.load(open(os.path.join(OUT, "manifest.json")))
    markers = [m["marker"] for m in man["models"]]

    norm = {}
    npath = os.path.join(OUT, "norm_params.csv")
    if os.path.exists(npath):
        for r in csv.DictReader(open(npath)):
            try:
                norm[r["marker"]] = (float(r["mean"]), float(r["std"]))
            except (TypeError, ValueError):
                pass

    out, missing = {}, []
    for mk in markers:
        p = os.path.join(OUT, mk + ".json")
        if not os.path.exists(p):
            missing.append(mk); continue
        d = json.load(open(p))
        age = np.asarray(d["age"], float)
        a0, a1 = float(age[0]), float(age[-1])
        grid = np.exp(np.linspace(np.log(a0), np.log(a1), NAGE))
        rec = {"a0": round(a0, 4), "a1": round(a1, 4),
               "nu": d["nu"], "tau": d["tau"],
               "site_sd": d.get("site_sd", 0.0)}
        for sex in ("m", "f"):
            mu = np.interp(grid, age, np.asarray(d[sex]["mu"], float))
            sg = np.interp(grid, age, np.asarray(d[sex]["sigma"], float))
            rec[sex] = [[round(float(v), 5) for v in mu],
                        [round(float(v), 5) for v in sg]]
        if mk in norm:
            rec["mean"], rec["std"] = norm[mk]
        out[mk] = rec

    bundle = {"n_age": NAGE, "markers": out}
    dst = os.path.join(OUT, "profile_bundle.json")
    with open(dst, "w") as fh:
        json.dump(bundle, fh, separators=(",", ":"))
    mb = os.path.getsize(dst) / 1e6
    print(f"bundle : {len(out)} marqueurs, {NAGE} ages, {mb:.2f} Mo")
    if missing:
        print("JSON absents :", missing)


if __name__ == "__main__":
    main()
