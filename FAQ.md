# D-EEG charts: frequently asked questions

Short answers to the questions the page cannot answer on its own. The model
specification and the scoring formula are in the [README](README.md).

## What does the page compute?

For each measure it places a value on the normative trajectory fitted on
typically developing participants and returns a z-score and a centile. The
trajectory is a GAMLSS model with a SHASHo2 distribution, so the z-score is
exact under skew and heavy tails rather than a simple (value − median) / SD.

## Which files can I upload?

Three shapes, all CSV, all read in the browser:

- **One subject.** Two columns, `marker,value`, one row per measure. Age and
  sex come from the form.
- **A cohort.** One row per subject with `age`, `sex` (M/F, or `female_bin`),
  optionally `ratio_ch_good`, any `true`/`false` column as a diagnosis group,
  and one column per measure.
- **Controls for calibration.** Four columns, `marker,value,age,sex`, one row
  per control and per measure.

Templates and synthetic examples for all three are under
[`static/data/templates/`](static/data/templates/) and in the page's
*Download files* panel. Measure names are the model's own, for instance
`periodic_alpha_mean` or `bankssts_lh_alpha_periodic`; the full list is the
`marker` column of `norm_params.csv`.

## Is my data sent anywhere?

No. The page is static, loads only its own parameter files, and uses no
third-party script or analytics. Values are read from the file you choose,
scored in your browser, and discarded when the tab closes. You can confirm
this in the browser's network panel or by running the page offline.

## Can I score a single subject from a site that was not in the training set?

Yes, with a caveat. The curves are those of the average training site. A new
site carries an offset the page cannot know from one subject, and between-site
differences range from negligible to about one control SD depending on the
measure. The z-score is therefore relative to the average site. Calibrating on
local controls removes most of that uncertainty.

## How many controls do I need to calibrate?

The page estimates one offset per measure from your controls and shrinks it
toward zero as the sample gets small, so few controls move the curves little
rather than wrongly. It needs at least 5 controls per measure to estimate
anything, and flags any estimate made from fewer than 100 controls, the number
Bethlehem et al. recommend for the MRI brain charts. Sexes are pooled, and the
offset shifts the location of the curve only.

## Do my features have to come from the same pipeline?

Yes. The reference values were extracted with PPSPrep and Markers using the
settings in the paper: 2-second epochs, average reference, eLORETA source
reconstruction, Desikan-Killiany parcellation, SpecParam for the aperiodic
fit. Another epoch length, reference, montage or unit convention gives
plausible numbers that are wrong by an unknown amount. The page does not
detect this.

## What is the quality field, and should I fill it?

It is the proportion of channels retained after preprocessing, between 0 and
1. When you give it, the page switches to the model fitted with that
covariate and evaluates it at your value. Leave it empty if you do not know
it: the default model is calibrated for a recording of unknown quality, and a
guessed value would narrow the curves without correcting your own offset.

## What do the shaded bands mean?

They are centiles of the reference population at each age: the darker band is
the 25th to 75th centile, the lighter one the 5th to 95th. They describe how
much typically developing participants vary, not how uncertain the curve is.
No confidence interval on the curve itself is published yet.

## Why does a chart say the model did not converge?

The fit stopped at the GAMLSS iteration cap before meeting the convergence
criterion. This affects 13 of the 20 whole-brain models. Their curves are
smooth and their residuals pass the diagnostics shown under the chart, but
treat them as provisional. The flag is `conv` in `manifest.json`.

## Are regional measures supported?

Yes. Each of the 68 Desikan-Killiany regions has its own model for periodic
power, connectivity, entropy and the aperiodic parameters, selectable under
*Region* on the charts. A cohort or profile file with regional columns is
mapped onto the atlas under *Brain maps › Your own map*; for a cohort, regions
surviving FDR are outlined. Per-region fit diagnostics are not shown.

## Which ages are covered?

From about 6 months to 66 years. Below 2 years the models rest on few
participants, so the bands widen and the amplitude is poorly determined.
Values outside the range are reported as such and not scored.

## Can I get the fitted model objects or the training data?

Not from this repository. The GAMLSS objects embed the training data and are
excluded deliberately; what is published is the set of parameters needed to
score a new subject. Healthy Brain Network data are public, Simons Searchlight
data are available through SFARI Base, and the other cohorts on request to the
corresponding author of the paper.

## How do I cite this, and can I contribute data?

Cite the paper listed on the [project page](https://dub21.github.io/D-EEG/).
To contribute recordings acquired with a 128-channel EGI system, or to ask
for a feature, contact the corresponding author.

## Where is the code?

This repository holds the website and the scripts that turn fitted models
into the published parameter files. Preprocessing is
[PPSPrep](https://github.com/ppsp-team/PPSPrep), feature extraction is
[Markers](https://github.com/Dub21/Markers), and normative modelling is
[PyNM](https://github.com/ppsp-team/PyNM).
