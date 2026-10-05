#!/usr/bin/env python3
"""Templates and synthetic examples for the charts page.

Subjects are drawn from the published normative parameters, so the files hold
no participant data and score as expected on the page: a profile from the
reference distribution, a cohort from a new site with a known offset and a
known case shift, and a controls file from that same site for calibration.

Usage:
    python3 deeg-app/tools/make_examples.py [--root REPO] [--seed N]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path

REGIONAL_FAMILY = "periodic_alpha"
SITE_OFFSET = 0.2  # control SD added to mu for every subject of the synthetic site
CASE_SHIFT = {
    "periodic_delta_mean": 0.4,
    "periodic_gamma_mean": 0.4,
    "periodic_alpha_mean": -0.4,
    "periodic_beta_mean": -0.3,
    "cf_alpha_mean": -0.3,
}
REGIONAL_CASE_SHIFT = -0.3
N_PER_GROUP = 30
N_CONTROLS = 30
QC_RANGE = (0.45, 0.65)
PROFILE_AGE, PROFILE_SEX = 12.4, "f"


@dataclass(frozen=True)
class Model:
    marker: str
    a0: float
    a1: float
    nu: float
    tau: float
    mean: float
    std: float
    grid: dict[str, list[list[float]]]  # "m"/"f" -> [mu, sigma] on a log-age grid
    qc_ref: float | None
    qc_coef_mu: float
    qc_coef_sigma: float

    def params_at(self, age: float, sex: str, qc: float | None) -> tuple[float, float]:
        mu_g, sg_g = self.grid[sex]
        n = len(mu_g)
        t = (math.log(age) - math.log(self.a0)) / (math.log(self.a1) - math.log(self.a0)) * (n - 1)
        i = max(0, min(n - 2, math.floor(t)))
        w = t - i
        mu = mu_g[i] + w * (mu_g[i + 1] - mu_g[i])
        sigma = sg_g[i] + w * (sg_g[i + 1] - sg_g[i])
        if qc is not None and self.qc_ref is not None:
            d = qc - self.qc_ref
            mu += self.qc_coef_mu * d
            sigma *= math.exp(self.qc_coef_sigma * d)
        return mu, sigma

    def raw(self, z: float, age: float, sex: str, qc: float | None, shift: float = 0.0) -> float:
        """Raw value of a subject at standard-normal score z, mu shifted by `shift` control SD."""
        mu, sigma = self.params_at(age, sex, qc)
        y = mu + shift + sigma * self.tau * math.sinh((math.asinh(z) + self.nu) / self.tau)
        return y * self.std + self.mean


def model_from(marker: str, rec: dict) -> Model:
    return Model(
        marker=marker,
        a0=rec["a0"],
        a1=rec["a1"],
        nu=rec["nu"],
        tau=rec["tau"],
        mean=rec["mean"],
        std=rec["std"],
        grid={"m": rec["m"], "f": rec["f"]},
        qc_ref=rec.get("qc_ref"),
        qc_coef_mu=rec.get("qc_coef_mu") or 0.0,
        qc_coef_sigma=rec.get("qc_coef_sigma") or 0.0,
    )


def load_whole_brain(norm: Path, qc: bool) -> dict[str, Model]:
    recs = json.loads((norm / "profile_bundle.json").read_text())["markers"]
    out: dict[str, Model] = {}
    for name, rec in recs.items():
        if name.endswith("_qc") != qc:
            continue
        marker = name[:-3] if qc else name
        out[marker] = model_from(marker, rec)
    return out


def load_regional(norm: Path, family: str, qc: bool) -> dict[str, Model]:
    recs = json.loads((norm / "regional" / f"{family}.json").read_text())["markers"]
    kind, _, band = family.rpartition("_")
    out: dict[str, Model] = {}
    for key, rec in recs.items():
        if key.endswith("_qc") != qc:
            continue
        region = key[:-3] if qc else key
        marker = f"{region}_{band}_{kind}" if kind else f"{region}_{family}"
        out[marker] = model_from(marker, rec)
    return out


def is_whole_brain(marker: str) -> bool:
    return marker.endswith("_mean")


def case_shift(marker: str) -> float:
    if marker in CASE_SHIFT:
        return CASE_SHIFT[marker]
    return 0.0 if is_whole_brain(marker) else REGIONAL_CASE_SHIFT


def fmt(v: float) -> str:
    return f"{v:.6g}"


def write_csv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


def write_templates(out: Path, whole: list[str]) -> None:
    write_csv(out / "profile_template.csv", ["marker", "value"], [{"marker": m, "value": ""} for m in whole])
    write_csv(out / "cohort_template.csv", ["age", "sex", "ratio_ch_good", "TD", "ASD", *whole], [])
    write_csv(out / "controls_template.csv", ["marker", "value", "age", "sex"], [])


def write_profile(path: Path, rng: random.Random, plain: dict[str, Model]) -> None:
    rows = [
        {"marker": m, "value": fmt(model.raw(rng.gauss(0.0, 1.0), PROFILE_AGE, PROFILE_SEX, None))}
        for m, model in plain.items()
    ]
    write_csv(path, ["marker", "value"], rows)


def write_cohort(path: Path, rng: random.Random, plain: dict[str, Model], adj: dict[str, Model]) -> None:
    markers = list(plain)
    rows = []
    for group in ("TD", "ASD"):
        for _ in range(N_PER_GROUP):
            age = round(rng.uniform(6, 18), 1)
            sex = "F" if rng.random() < 0.5 else "M"
            qc = round(rng.uniform(*QC_RANGE), 2)
            row = {
                "age": str(age),
                "sex": sex,
                "ratio_ch_good": str(qc),
                "TD": str(group == "TD").lower(),
                "ASD": str(group == "ASD").lower(),
            }
            for m in markers:
                shift = SITE_OFFSET + (case_shift(m) if group == "ASD" else 0.0)
                row[m] = fmt(adj.get(m, plain[m]).raw(rng.gauss(0.0, 1.0), age, sex.lower(), qc, shift))
            rows.append(row)
    write_csv(path, ["age", "sex", "ratio_ch_good", "TD", "ASD", *markers], rows)


def write_controls(path: Path, rng: random.Random, whole: dict[str, Model]) -> None:
    rows = []
    for _ in range(N_CONTROLS):
        age = round(rng.uniform(5, 40), 1)
        sex = "F" if rng.random() < 0.5 else "M"
        for m, model in whole.items():
            rows.append(
                {
                    "marker": m,
                    "value": fmt(model.raw(rng.gauss(0.0, 1.0), age, sex.lower(), None, SITE_OFFSET)),
                    "age": str(age),
                    "sex": sex,
                }
            )
    write_csv(path, ["marker", "value", "age", "sex"], rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2], help="repository root")
    p.add_argument("--seed", type=int, default=20251005)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    norm = args.root / "static" / "data" / "normative"
    out = args.root / "static" / "data" / "templates"
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    plain = load_whole_brain(norm, qc=False) | load_regional(norm, REGIONAL_FAMILY, qc=False)
    adj = load_whole_brain(norm, qc=True) | load_regional(norm, REGIONAL_FAMILY, qc=True)
    plain = dict(sorted(plain.items(), key=lambda kv: (not is_whole_brain(kv[0]), kv[0])))
    whole = [m for m in plain if is_whole_brain(m)]
    write_templates(out, whole)
    write_profile(out / "example_profile.csv", rng, plain)
    write_cohort(out / "example_cohort.csv", rng, plain, adj)
    write_controls(out / "example_controls.csv", rng, {m: plain[m] for m in whole})
    print(f"wrote {len(list(out.glob('*.csv')))} files to {out} (seed {args.seed})")


if __name__ == "__main__":
    main()
