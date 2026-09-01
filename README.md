# D-EEG — Heterogeneity of brain dynamics in genetic and psychiatric conditions

Project page and interactive normative EEG charts.

- Paper page: <https://dub21.github.io/D-EEG/>
- Interactive charts: <https://dub21.github.io/D-EEG/charts.html>

## Interactive charts

`charts.html` serves normative trajectories for **1312 EEG measures** (34 cortical
regions × 2 hemispheres of the Desikan-Killiany atlas, plus whole-brain averages) and
places new subjects on them.

Everything runs in the browser. There is no server, no external library, and **no data
is ever uploaded** — a subject's values are read locally and discarded when the page
closes.

### Placing your own data

Two-column CSV, one row per measure, raw values:

```
marker,value
bankssts_lh_alpha_periodic,0.2919
precuneus_rh_gamma_connectivity,0.3746
```

Marker names are the model's own; the full list is in
[`static/data/normative/norm_params.csv`](static/data/normative/norm_params.csv).
Enter age and sex in the form, choose the file, and the page returns a z-score and a
centile per measure, sorted by deviation, with the full profile downloadable as CSV.

**Your features must come from the same pipeline as the reference cohort**
([PPSPrep](https://github.com/ppsp-team/PPSPrep) for preprocessing,
[Markers](https://github.com/Dub21/Markers) for feature extraction), with the settings
described in the paper: 2-second epochs, eLORETA source reconstruction, Desikan-Killiany
parcellation. A different epoch length, reference, or unit convention will still produce
plausible-looking z-scores that are wrong. The page does not currently check this.

## Normative models

Fitted with [PyNM](https://github.com/ppsp-team/PyNM) over R's GAMLSS, family
**SHASHo2** (four parameters: location, scale, skewness, kurtosis):

```
mu    ~ ps(age) + as.factor(female_bin) + as.factor(Unique_Site_ID) + ratio_ch_good
sigma ~ ps(age) + as.factor(female_bin) + as.factor(Unique_Site_ID) + ratio_ch_good
nu    ~ 1
tau   ~ 1
```

Trained on participants with no diagnosis (per-model n in `manifest.json`, typically
~1000, age 0.5 to 66 years). A subject's standardised score is

```
z = sinh( tau · asinh( (y - mu) / (sigma · tau) ) - nu )    with z ~ N(0,1)
centile = Phi(z)
```

where `y = (raw value - mean) / std` using the constants in `norm_params.csv`.

## Repository layout

```
index.html                       paper page
charts.html                      interactive charts, self-contained
static/data/normative/
  <marker>.json                  1312 files: mu and sigma over 200 ages × 2 sexes,
                                 nu, tau, per-site offsets
  manifest.json                  catalogue and per-model fit diagnostics
  norm_params.csv                1305 standardisation constants
  profile_bundle.json            all markers on a 120-point log age grid,
                                 fetched only when a profile CSV is submitted
static/data/stats/               Shapiro, SMSE, MSLL, skewness, kurtosis per model
deeg-app/tools/                  scripts that turn fitted .rds models into the JSON above
deeg-app/refit/PATCH.md          specification for refitting site as a random effect
```

Fitted `.rds` model objects are **not** in this repository and must not be added: they
embed the training data (per-participant values, ages, sites and residuals). `.gitignore`
blocks them.

## Known limitations

- **Site effect.** Site is a fixed effect, so it cannot be extrapolated to a site the
  model never saw. Its median range across models is about **1 standard deviation**, and
  exceeds 1 SD for 53% of them. The published curves are centred on the mean of the
  training sites, which is neutral but leaves that uncertainty on any new subject's
  centile. `deeg-app/refit/PATCH.md` specifies the random-effect refit that would fix it.
- **Seven whole-brain measures have no standardisation constants**
  (`exponent_mean`, `offset_mean`, the five `entropy_*_mean`) and accept
  already-standardised values only.
- **Convergence.** 54% of models reached the GAMLSS 20-iteration cap without meeting the
  convergence criterion. The measured consequence is an internal inconsistency in sigma
  of about 0.002% (median), and residual diagnostics are indistinguishable between
  converged and non-converged models. Affected models are flagged in the page.
- **Data quality.** Models include the proportion of good channels as a covariate. When a
  submitted profile does not provide it, the training median is used.

## Related repositories

- [PPSPrep](https://github.com/ppsp-team/PPSPrep) — preprocessing
- [Markers](https://github.com/Dub21/Markers) — feature extraction
- [PyNM](https://github.com/ppsp-team/PyNM) — normative modelling

## Acknowledgments
Parts of this project page were adopted from the [Nerfies](https://nerfies.github.io/) page.

## Website License
<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
