#!/usr/bin/env python3
"""Carte cerebrale moyenne d'un groupe, a partir des z deja calcules.

Lit un fichier de z-scores par sujet, moyenne par region au sein du groupe
demande, et colorie l'atlas DK. Tourne en local : aucune donnee individuelle
ne sort, seules des moyennes de groupe sont ecrites.

Variables d'environnement :
  Z       csv de z-scores            (defaut ~/mnt/Downloads/z_scores_randomqc.csv)
  ATLAS   svg de l'atlas             (defaut ~/mnt/D-EEG/static/images/dk_atlas.svg)
  OUT     dossier de sortie          (defaut ~/mnt/Desktop/brainmaps_group)
  GROUP   colonne du groupe          (defaut "Autism Spectrum Disorder")
  REF     colonne de reference       (defaut "No Diagnosis Given")
  FAMS    familles, separees par ,   (defaut toutes celles presentes)
  LIM     echelle +/- , 0 = auto     (defaut 0)
"""
import os, re, csv, math
import numpy as np, pandas as pd

Z     = os.path.expanduser(os.environ.get("Z", "~/mnt/Downloads/z_scores_randomqc.csv"))
ATLAS = os.path.expanduser(os.environ.get("ATLAS", "~/mnt/D-EEG/static/images/dk_atlas.svg"))
OUT   = os.path.expanduser(os.environ.get("OUT", "~/mnt/Desktop/brainmaps_group"))
GROUP = os.environ.get("GROUP", "Autism Spectrum Disorder")
REF   = os.environ.get("REF", "No Diagnosis Given")
FAMS  = [f for f in os.environ.get("FAMS", "").split(",") if f]
LIM   = float(os.environ.get("LIM", "0"))

SUF = "_zscore"
RE_COL = re.compile(r"^(.+?_(?:lh|rh))_(.+)" + SUF + r"$")

# bleu <-> gris <-> rouge, point neutre gris (theme clair)
NEG, MID, POS = (42, 120, 214), (240, 239, 236), (227, 73, 72)
ABSENT = "#e7e6e1"


def color(v, lim):
    t = max(-1.0, min(1.0, v / lim)) if lim else 0.0
    end = NEG if t < 0 else POS
    a = abs(t)
    return "rgb(%d,%d,%d)" % tuple(round(m + (e - m) * a) for m, e in zip(MID, end))


def bh(p):
    """Benjamini-Hochberg : renvoie les q-valeurs."""
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def main():
    from scipy import stats
    head = pd.read_csv(Z, nrows=0).columns.tolist()
    fams = {}
    for c in head:
        m = RE_COL.match(c)
        if m:
            fams.setdefault(m.group(2), {})[m.group(1)] = c
    if FAMS:
        fams = {f: v for f, v in fams.items() if f in FAMS}
    fams = {f: v for f, v in fams.items() if len(v) >= 20 and "knee" not in f}

    meta = [GROUP, REF]
    cols = meta + [c for v in fams.values() for c in v.values()]
    df = pd.read_csv(Z, usecols=[c for c in cols if c in head], low_memory=False)
    g = df[GROUP] == True
    r = df[REF] == True
    print(f"{GROUP} : n={int(g.sum())}   {REF} : n={int(r.sum())}")
    if not g.sum():
        raise SystemExit("groupe vide")

    svg0 = open(ATLAS).read()
    os.makedirs(OUT, exist_ok=True)
    summary = []

    for fam, regions in sorted(fams.items()):
        rows = []
        for reg, col in regions.items():
            a = df.loc[g, col].dropna().values
            b = df.loc[r, col].dropna().values
            if len(a) < 20:
                continue
            t, p = stats.ttest_ind(a, b, equal_var=False)
            rows.append((reg, a.mean(), b.mean(), a.mean() - b.mean(), len(a), p))
        if not rows:
            continue
        q = bh([x[5] for x in rows])
        lim = LIM or max(0.05, np.ceil(max(abs(x[3]) for x in rows) * 20) / 20)

        svg = svg0
        for (reg, ma, mb, d, n, p), qq in zip(rows, q):
            pat = re.compile(r'(<path class="reg" id="%s" )' % re.escape(reg))
            style = 'style="fill:%s%s" ' % (color(d, lim),
                                            ';stroke:#111;stroke-width:4' if qq < .05 else '')
            svg = pat.sub(r'\1' + style, svg, count=1)
            summary.append({"family": fam, "region": reg, "n": n,
                            "mean_z_group": round(ma, 4), "mean_z_ref": round(mb, 4),
                            "diff": round(d, 4), "p": f"{p:.3g}", "q_fdr": f"{qq:.3g}"})
        svg = svg.replace('style="fill:', 'style="fill:')       # no-op lisible
        svg = svg.replace('fill="#3a3a3a"', 'fill="%s"' % ABSENT)
        stops = "".join('<stop offset="{}%" stop-color="{}"/>'.format(
            k * 10, color(-lim + 2 * lim * k / 10, lim)) for k in range(11))
        legend = (
            '<g transform="translate(60,1210)">'
            '<defs><linearGradient id="lg">' + stops + '</linearGradient></defs>'
            '<rect x="0" y="0" width="1480" height="26" rx="13" fill="url(#lg)"/>'
            '<text x="0" y="62" font-family="sans-serif" font-size="34" fill="#444">'
            + "{:+.2f}".format(-lim) + '</text>'
            '<text x="740" y="62" font-family="sans-serif" font-size="34" fill="#444" text-anchor="middle">0</text>'
            '<text x="1480" y="62" font-family="sans-serif" font-size="34" fill="#444" text-anchor="end">'
            + "{:+.2f}".format(lim) + '</text>'
            '<text x="740" y="-14" font-family="sans-serif" font-size="34" fill="#222" text-anchor="middle">'
            + "{} &#8212; {} vs {}, difference de z moyen (contour : q&lt;0.05)".format(fam, GROUP, REF)
            + '</text></g>')
        svg = svg.replace('viewBox="-6 -6 1618 1188"', 'viewBox="-6 -6 1618 1300"')
        svg = svg.replace('</svg>', legend + '</svg>')
        dst = os.path.join(OUT, fam + ".svg")
        open(dst, "w").write(svg)
        nsig = int((q < .05).sum())
        print(f"  {fam:22s} {len(rows):3d} regions, echelle +/-{lim:.2f}, "
              f"{nsig:3d} significatives apres FDR")

    with open(os.path.join(OUT, "summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader(); w.writerows(summary)
    print(f"\n{len(fams)} cartes dans {OUT}")


if __name__ == "__main__":
    main()
