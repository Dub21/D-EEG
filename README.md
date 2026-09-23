# D-EEG — Heterogeneity of brain dynamics in genetic and psychiatric conditions

Project page and interactive normative EEG charts.

- Paper page: <https://dub21.github.io/D-EEG/>
- Interactive charts: <https://dub21.github.io/D-EEG/charts.html>

## Interactive charts

`charts.html` serves normative trajectories for **20 whole-brain EEG measures** and
places new subjects on them: periodic power, connectivity (wPLI) and permutation entropy
on five frequency bands, plus the aperiodic exponent and offset and the three alpha-peak
parameters.

Each measure comes in two fitted versions. The default carries no data-quality covariate
and is calibrated for a recording whose quality is unknown. The second includes the
proportion of retained channels as a covariate in mu and sigma; the page uses it only when
that value is supplied, and never substitutes a default for it — plugging in the training
median would narrow sigma without correcting the subject's own quality offset.

Regional models (34 Desikan-Killiany regions × 2 hemispheres) are being refitted to the
specification below and are not published at present.

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
[`static/data/normative/manifest.json`](static/data/normative/manifest.json).
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
mu    ~ ps(age) + as.factor(female_bin) + random(as.factor(Unique_Site_ID))
sigma ~ ps(age) + as.factor(female_bin) + random(as.factor(Unique_Site_ID))
nu    ~ 1
tau   ~ 1
```

Site is a **random** intercept, so the published curves are those of the average site and
a new site can be placed on them. Fitted with `n.cyc = 200`; 20 of the 25 models meet the
convergence criterion.

A second set of models adds the proportion of retained channels as a parametric covariate
in mu and sigma. Its effect is marker-dependent: negligible for beta and theta power,
but worth up to 0.47 SD at the extremes of the quality distribution for delta connectivity
and gamma power.

Trained on participants with no diagnosis (per-model n in `manifest.json`, typically
~1000, age 0.5 to 66 years). A subject's standardised score is

```
z = sinh( tau · asinh( (y - mu) / (sigma · tau) ) - nu )    with z ~ N(0,1)
centile = Phi(z)
```

where `y = (raw value - centre) / scale`, the centre being the **median** of the
control values and the scale their standard deviation, using the constants in
`norm_params.csv`. The column is named `mean` there for backward compatibility.

## Repository layout

```
index.html                       paper page
charts.html                      interactive charts, self-contained
static/data/normative/
  <marker>.json                  mu and sigma over 200 ages × 2 sexes, nu, tau,
                                 per-site offsets (BLUPs of the random intercept)
  manifest.json                  catalogue and per-model fit diagnostics
  norm_params.csv                standardisation constants
  profile_bundle.json            all markers on a 120-point log age grid,
                                 fetched only when a profile CSV is submitted
static/data/stats/               Shapiro, SMSE, MSLL, skewness, kurtosis per model
deeg-app/tools/extract_json.py   fitted .rds models -> per-marker JSON
deeg-app/tools/build_manifest.py JSON present -> manifest.json
```

Fitted `.rds` model objects are **not** in this repository and must not be added: they
embed the training data (per-participant values, ages, sites and residuals). `.gitignore`
blocks them.

## Known limitations

- **Site effect.** Site is a random intercept estimated on 14 training sites, with a
  between-site standard deviation of about **0.3 SD**. The published curves are those of
  the average site. A new site's own offset is not known, and that uncertainty remains on
  any new subject's centile until it is estimated from local controls.
- **Convergence.** 13 of the 20 published models reach the GAMLSS iteration cap without
  meeting the convergence criterion. Affected models are flagged under the chart.
- **Age coverage.** Below about 2 years the models rest on a handful of participants, so
  the confidence band widens sharply there and the amplitude is not well determined.
- **Regional models.** Only whole-brain measures are published while the regional models
  are refitted to this specification.

## Related repositories

- [PPSPrep](https://github.com/ppsp-team/PPSPrep) — preprocessing
- [Markers](https://github.com/Dub21/Markers) — feature extraction
- [PyNM](https://github.com/ppsp-team/PyNM) — normative modelling

## Acknowledgments
Parts of this project page were adopted from the [Nerfies](https://nerfies.github.io/) page.

## Website License
<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
