# Deconvolve a GC/MS chromatogram and export its compound spectra

Read an ANDI-MS netCDF chromatogram, estimate a smoothly time-varying
point spread function, recover log-concave compounds with
[`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html),
and write a NIST-compatible MSP spectral library. All original mass
channels are retained in the final spectra.

## Usage

``` r
# S3 method for class 'GCMSdeco'
residuals(
  object,
  type = c("pearson", "response"),
  time_range = NULL,
  mass = NULL,
  ...
)

GCMSdeco(
  cdf_file,
  msp_file = NULL,
  family = quasipoisson("identity"),
  lambda = "gic",
  nwindows = 10L,
  alkanes_file = NULL,
  time_range = NULL,
  calibration_range = NULL,
  baseline_df = NULL,
  precluster = 32L,
  dispersion = NULL,
  theta = NULL,
  correlation = c("auto", "none", "lag1"),
  ncores = NULL,
  control = list(),
  psf_control = list(),
  verbose = FALSE,
  xlsx_file = NULL,
  sample_name = NULL,
  library = NULL,
  search = list(),
  output_dir = NULL,
  coupled = TRUE,
  batch = list()
)

# S3 method for class 'GCMSdeco'
coef(object, ..., normalize = c("none", "l1"))

# S3 method for class 'GCMSdeco'
fitted(object, ...)

# S3 method for class 'GCMSdeco'
nobs(object, ...)

# S3 method for class 'GCMSdeco'
logLik(object, ...)

# S3 method for class 'GCMSdeco'
AIC(object, ..., k = 2)

# S3 method for class 'GCMSdeco'
BIC(object, ...)

# S3 method for class 'GCMSdeco'
GIC(object, ...)

# S3 method for class 'GCMSdeco'
EBIC(object, ...)

# S3 method for class 'GCMSdeco'
MBIC(object, ...)

# S3 method for class 'GCMSdeco'
print(x, ...)
```

## Arguments

- object, x:

  A fitted `GCMSdeco` object.

- type:

  Residual type: `"pearson"` uses the fitted marginal variance;
  `"response"` returns observed minus fitted latent intensity.

- time_range:

  Optional increasing pair of sample retention-time bounds in minutes.
  Infinite bounds are accepted. This can exclude the solvent interval or
  select a chromatographic window.

- mass:

  Optional nominal m/z values for residual extraction.

- ...:

  Further method arguments.

- cdf_file:

  A directory of sample subdirectories with a `standards` subdirectory,
  or a local ANDI-MS `.cdf` path or HTTP(S) URL. The optional ncdf4
  package is required. Masses are binned to the nearest integer m/z;
  traces between the minimum and maximum recorded mass are retained.

- msp_file:

  Output `.msp` path. Parent directories are created if necessary; an
  existing output file is replaced.

- family:

  Identity-link noise model accepted by
  [`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html).
  The default `quasipoisson("identity")` estimates one pooled dispersion
  across channels, suitable for overdispersed instrumental intensities.

- lambda:

  Information-criterion preset or a nonnegative numeric penalty. For an
  individual chromatogram the numeric value is the total charge per
  compound. The directory workflow uses this charge in each initial fit.
  Its joint refinement uses the numeric value per selection degree of
  freedom. Both charges precede any temporal adjustment.

- nwindows:

  Number of overlapping blind-PSF windows, default `10L`. `1L` uses a
  shared PSF. The returned sparse design interpolates the PSFs smoothly
  across the chromatogram.

- alkanes_file:

  Optional ANDI-MS alkane-standard chromatogram. Isolated peaks
  calibrate the initial PSF through
  [`logconNMF::estimate_psf()`](https://rdrr.io/pkg/logconNMF/man/estimate_psf.html),
  with automatic support, sample refinement and resampling to the sample
  scan interval.

- calibration_range:

  Optional retention-time bounds in minutes for the standard. The whole
  standard is used when omitted.

- baseline_df:

  Dimension of the nonnegative cubic B-spline background. Use `1L` for a
  constant background in a short window. A spline uses at least four
  basis functions. `NULL` chooses a low-dimensional basis from the scan
  count and typical PSF width; an integer of at least four fixes its
  dimension.

- precluster:

  A saved
  [`channel_groups()`](https://rdrr.io/pkg/logconNMF/man/channel_groups.html)
  map shares the same aggregation across analyses. `NULL` or `FALSE`
  uses every channel throughout. `TRUE` uses up to 32 groups; an integer
  supplies another maximum. Size-adjusted Poisson clustering,
  aggregation, selection and fusion operate in that reduced space,
  followed by refinement on all original channels.

- dispersion:

  Known common Gaussian variance, quasipoisson scale, Gamma dispersion,
  NB1 variance-to-mean ratio or common Gamma scale. `NULL` estimates a
  pooled value from temporal differences and freezes it. Parameters
  supplied in a
  [signal_families](https://tomwenseleers.github.io/L0glm/reference/signal_families.html)
  object are used when this argument is omitted. Poisson and NB2 use
  unit dispersion.

- theta:

  NB2 shape, shared by all original channels. `NULL` estimates and
  freezes a temporal moment estimate; a positive value fixes it; a
  negative value initializes at its absolute value and enables pooled
  likelihood updates. A fixed shape in the supplied NB family is
  respected.

- correlation:

  `"auto"` or `"lag1"` estimates an optional lag-one working adjustment
  to the selection penalty. `"none"` uses independent observations; a
  number between zero and 0.49 supplies the lag-one correlation.

- ncores:

  Total logical-CPU budget. `NULL` uses physical cores. Directory
  analyses run up to `min(samples,ncores)` sample workers, dividing the
  native thread budget between them. With at least 16 samples,
  `ncores=16L` uses 16 workers with one native thread each. BLAS threads
  are limited during batch fitting and caller settings are restored. For
  a single chromatogram, this sets its native OpenMP thread count.

- control:

  Named list: `broad_repeats=1L`, `full_sweeps=1L`, `final_polishes=1L`,
  `native_sweeps=6L`, `block_sweeps=6L`, `psf_width=NULL`,
  `screen_factor=1`, `censor=0.5`, `linkage="ward"`,
  `joint_trials=NULL`, `profile_cycles=4L`, `fusion_passes=0L`,
  `profile_fusion=TRUE`. `joint_trials` uses direct Poisson fusion for
  Poisson data and joint local refinement for the other families.
  `psf_width` supplies the typical FWHM in scans. When omitted, local
  updates use the middle design column and broad-profile fusion uses the
  median of up to 21 representative design columns. `screen_factor`
  scales the initial selection penalty. `censor` is the Gamma detection
  limit: recorded zero intensities are left-censored below this value;
  positive values must be at least the limit. Zero requests an
  uncensored Gamma model. A supplied Gamma family constructor provides
  its own detection limit unless overridden here. `linkage` may be
  `"ward"` or `"bregman"`. `joint_trials` allows overlapping compounds
  to adjust together during local fusion refinement. `profile_cycles`
  bounds the smooth shape/background initialization; positive values
  also enable fine native profile/spectrum refinement. Zero skips these
  shape-polishing cycles. `fusion_passes` optionally adds broad-anchor
  and compatible-tail searches during initialization. `profile_fusion`
  enables algebraic profile fusion, retained contained-fragment
  proposals and a final broad-fragment pass on every original channel.
  Broad fusions are compared with an equally refined unmerged local
  model. The entire support of every updated component is included in
  the local fit, with the remaining components and background held in an
  offset. `original_cleanup=TRUE` enables a final original-channel
  deletion, local fusion and residual-discovery pass. Set it to `FALSE`
  to skip that pass.

- psf_control:

  Named PSF settings. Entries `method`, `standard`, `precluster`,
  `weights` and `offset` are passed directly to
  [`logconNMF::estimate_psf()`](https://rdrr.io/pkg/logconNMF/man/estimate_psf.html).
  Other entries form its `control` list, for example
  `list(kernel_est=psf, maxiter=0L)` keeps a supplied PSF fixed and
  `list(L=61L, maxiter=8L, tol=1e-3)` estimates a PSF on the sample.
  `list(standard=attr(GCMS_data2, "alkanes_TIC"))` uses the bundled
  calibration TIC; its numeric names supply retention times.
  `precluster=1L` uses the total-ion signal, `FALSE` the original
  channels, and `32L` at most 32 groups. A fitted
  `estimate_psf(method="standard")` result can be passed as `standard`
  with `method="standard"` to reuse a solvent-trimmed ladder calibration
  across samples, rebuilding each sparse design on its own grid. The PSF
  uses the same noise family and fixed noise estimates as the final
  factorization, with its information adjustment included in the
  weights.

- verbose:

  Print stage timings, objective changes and compound counts.

- xlsx_file:

  Optional Excel output path. The workbook contains compound areas and
  one spectrum column on the first sheet, followed by run metadata
  including the sample name. Fragments use `m/z (intensity)`, sorted by
  decreasing intensity and scaled to base peak 9999; fragments rounding
  to zero are omitted from this display. Requires the optional openxlsx
  package. `NULL` omits Excel export.

- sample_name:

  Sample identifier, defaulting to the CDF filename without its
  extension. Included in the compound table, MSP entries and workbook.

- library:

  Optional GCMSlibrary for retention-index-constrained identification.

- search:

  Named options for library_search(), such as phase="DB5",
  ri_tolerance=30, n_matches=5 and method="stein_scott". The alkane file
  supplies RI calibration automatically; a precomputed calibration can
  be supplied here. Identification is stored in the returned object and
  included in MSP/Excel exports.

  Directory input returns a
  [GCMSbatch](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  object, with acquisition-ordered sample groups, per-sample exports and
  a consensus area table and spectral library.

- output_dir:

  Directory for batch results; defaults to a sibling directory named
  after the input folder with suffix `_GCMSprime`.

- coupled:

  Share spectra and refine sample-specific profiles with
  [`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html)
  for directory input.

- batch:

  Named directory options; see
  [GCMSbatch](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md).

- normalize:

  Spectral coefficient normalization: `"none"` returns integrated
  channel intensities; `"l1"` gives a unit-sum spectrum per compound,
  with its magnitude carried by the compound area.

- k:

  Complexity multiplier for
  [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html).

## Value

An object of class `GCMSdeco` containing `model` (the full
factorization), a `compounds` table with names, retention times,
original-file apex scan indices, integrated intensities and ion counts,
the unified PSF estimate and sparse design in `psf`, optional standard
calibration in `psf_standard`, `background_df`, `psf_width_scans`, an
average precision diagnostic, retention-time and mass axes, original
scan indices, sample identity, output file paths and stage timings.
Coefficient, fitted-value, log-likelihood and information-criterion
methods delegate to `model`.

## Details

**File formats.** ANDI-MS variables `scan_index`, `point_count`,
`scan_acquisition_time`, `mass_values` and `intensity_values` define the
scan spectra. Scan indices are zero-based. Padding between scans is
accepted. Intensities must be finite and nonnegative. Times are
converted to minutes; absent units use the ANDI-MS seconds convention.
The selected retention-time interval is read directly from the point
arrays.

MSP entries are ordered by profile apex time, with a unique compound
name, retention time and integrated ion intensity in the comments,
followed by mass/intensity pairs. Exported spectra have base peak 999
and retain every positive fitted ion.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns unnormalized
fitted spectra. With `normalize="l1"`,
[`coef()`](https://rdrr.io/r/stats/coef.html) returns spectra summing to
one. The MSP and R object retain every positive fitted ion. The
spreadsheet's compact spectrum text is a rounded display alongside the
integrated areas. Their compound profiles are \\P\_{ir}=H\_{ir}\sum_j
S\_{rj}\\, preserving \\HS=PT\\. The workbook reports `area`, equal to
\\\sum_j S\_{rj}\\, in intensity-scans, and `area_min`, the trapezoidal
integral of \\P\_{\cdot r}\\ over retention time, in intensity-minutes.

The returned object stores the factors and source CDF path. Observed
data are loaded from the CDF only when a plot needs them. Fitted signals
are reconstructed from the factors; the observed matrix is released
after fitting.

**PSF and background.**
[`logconNMF::estimate_psf()`](https://rdrr.io/pkg/logconNMF/man/estimate_psf.html)
estimates constant or smoothly time-varying PSFs from the sample,
optionally starting from an alkane standard. Original channels, outcome
groups and the total-ion signal are supported. Interpolated columns are
projected to log concavity on their local support and normalized to unit
sum. Their typical full width at half maximum, \\w\\, sets the
background dimension to \$\$q=\max\\4,\min\[40,3+\lceil
n/(200w)\rceil\]\\.\$\$ Background coefficients are nonnegative and
unpenalized. This keeps the baseline much smoother than an individual
elution peak. The calibrated median PSF width also sets the guard for
the final broad-profile fusion pass in
[`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html).
Set `control=list(profile_fusion=FALSE)` to disable that pass.

**Shared noise model.** Noise parameters are estimated once from the
original channels and reused in PSF estimation and factorization. The
PSF weights include the information adjustment \\1/\kappa\\, where
\\\kappa=1+2\widehat\rho_1\\. Aggregation adjusts the variance
parameters according to the selected family; NB2 and ordinary Gamma use
moment matching. See
[`logconNMF::estimate_psf()`](https://rdrr.io/pkg/logconNMF/man/estimate_psf.html)
for the distributional details. The returned `psf_working_weight`
summarises the average Poisson-scale precision for comparing noise
scales; the actual PSF likelihood uses the specified family.

`precluster=32L` performs selection and fusion in at most 32
size-adjusted outcome groups, followed by refinement on every original
channel. The final original-channel cleanup checks deletions, local
fusions and missed residual peaks using the complete compound-count
criterion. Its diagnostics are in `model$cleanup`;
`control=list(original_cleanup=FALSE)` skips it. See
[`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html)
for the selection criterion and shape constraints. The backend uses a
permissive component initialization, native log-concave profile/spectrum
refinement, and final rank selection with the requested
working-correlation adjustment. Small fragments contained in a broad
peak are assessed against equally refined unmerged local fits on every
original ion channel. Trial diagnostics are in `model$contained_fusion`;
the initialization criterion and history are in `model$initialization`.
The spreadsheet and returned coefficient rows use increasing apex
retention time. Rank selection, hard log concavity, censored Gamma
intensities and compound-count information criteria follow
[`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html).

## Experimental GC/MS examples

The [GC/MS deconvolution
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)
converts the bundled
[GCMS_data1](https://rdrr.io/pkg/logconNMF/man/GCMS_data1.html) and
[GCMS_data2](https://rdrr.io/pkg/logconNMF/man/GCMS_data2.html) matrices
to CDF, fits the full-channel models and exports their spectra. It
includes timings, information criteria and all plot types, alongside the
algorithm. **GCMS_data1: overlapping cuticular hydrocarbons.**
![Background-corrected observed mass traces and fitted
compounds.](figures/gcms-data1-traces.png)![Observed, corrected, fitted
and Pearson-residual
heatmaps.](figures/gcms-data1-heatmaps.png)![Recovered elution profiles
and mass spectra.](figures/gcms-data1-profiles-spectra.png)**GCMS_data2:
beer volatiles with an asymmetric acid peak.** ![Observed ion traces and
fitted beer-volatile
components.](figures/gcms-data2-traces.png)![Full-channel signal and
Pearson-residual
heatmaps.](figures/gcms-data2-heatmaps.png)![Representative compound
profiles and mass spectra.](figures/gcms-data2-profiles-spectra.png)

## References

NIST MS Search user manual, section on NIST text format:
<https://chemdata.nist.gov/mass-spc/ms-search/docs/Ver20Man.pdf>.

## See also

[`write_GCMS_cdf()`](https://tomwenseleers.github.io/GCMSprime/reference/write_GCMS_cdf.md),
[`plot.GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/plot.GCMSdeco.md),
[`logconNMF::logconNMF()`](https://rdrr.io/pkg/logconNMF/man/logconNMF-package.html),
[`logconNMF::estimate_psf()`](https://rdrr.io/pkg/logconNMF/man/estimate_psf.html),
[`L0glm_blind_deco()`](https://rdrr.io/pkg/logconNMF/man/L0glm_blind_deco.html)

## Examples

``` r
library(logconNMF)
if (FALSE) { # \dontrun{
library(L0glm)
data(GCMS_data1)
data(GCMS_data2)
out <- file.path(tempdir(), "L0glm-GCMS")
dir.create(out, showWarnings = FALSE)
write_GCMS_cdf(GCMS_data1, file.path(out, "GCMS_data1.cdf"))
write_GCMS_cdf(GCMS_data2, file.path(out, "GCMS_data2.cdf"))

fit1 <- GCMSdeco(file.path(out, "GCMS_data1.cdf"),
  file.path(out, "GCMS_data1.msp"), nwindows = 1L, precluster = 16L,
  baseline_df = 1L,
  correlation = "auto", ncores = 8L,
  psf_control = list(kernel_est = convkernel("gaussian", 65L,
    list(sigma = 8)), maxiter = 0L))
fit2 <- GCMSdeco(file.path(out, "GCMS_data2.cdf"),
  file.path(out, "GCMS_data2.msp"), nwindows = 1L, precluster = 32L,
  correlation = "auto", ncores = 8L,
  psf_control = list(kernel_est = convkernel("gaussian", 65L,
    list(sigma = 6)), maxiter = 0L))
fit1$compounds
GIC(fit2)
coef(fit2, normalize = "l1")

plot(fit1, type = "lines", nrows = 1L,
  minutes_per_panel = diff(range(fit1$time)) + .01, labels = TRUE)
plot(fit1, type = "heatmap")
plot(fit1, type = "profiles", nrows = 3L)
plot(fit2, type = "lines", file = file.path(out, "GCMS_data2_traces.pdf"),
  nrows = 1L, minutes_per_panel = diff(range(fit2$time)) + .01)
plot(fit2, type = "heatmap", file = file.path(out, "GCMS_data2_heatmaps.pdf"))
plot(fit2, type = "profiles", file = file.path(out, "GCMS_data2_profiles_and_spectra.pdf"))

# For a full chromatogram with a separate calibration standard:
beer <- GCMSdeco("beer.cdf", "beer.msp", alkanes_file = "alkanes.cdf",
  xlsx_file = "beer.xlsx", nwindows = 10L, precluster = 32L)

# Automatic PSF estimation from a bundled window and its named alkane TIC.
psf <- estimate_psf(GCMS_data2, nwindows = 3L)
dim(psf$X)

# A CDF contains the chromatogram; pass the R calibration attribute explicitly.
calibrated <- GCMSdeco(file.path(out, "GCMS_data2.cdf"),
  file.path(out, "GCMS_data2_calibrated.msp"), nwindows = 3L, precluster = 32L,
  psf_control = list(standard = attr(GCMS_data2, "alkanes_TIC")))
} # }

if (FALSE) { # \dontrun{
# Sixteen related samples, retaining all mass channels. Compounds occur in
# half the samples on average, with retention drift and concentration-dependent
# widths and tails. Here the simulated quasipoisson dispersion is known.
library(logconNMF)
sim <- simulate_GCMS_samples(samples=16L, n=1200L, k=48L,
  preset="time_area", presence=0.5, irregular=0, drift=2,
  family="quasipoisson", dispersion=2, seed=9118)
input <- file.path(tempdir(), "GCMSprime_samples16")
dir.create(input, showWarnings=FALSE)
for(s in sim$samples) {
  write_GCMS_cdf(s$Y, file.path(input, paste0(s$sample_id, ".cdf")),
    time=s$time, acquired_at=s$acquired_at)
}
rm(sim, s)                         # observations are read from the CDF files
experiment <- GCMSdeco(input, output_dir=paste0(input, "_results"),
  family=quasipoisson("identity"), dispersion=2, correlation="none",
  nwindows=2L, baseline_df=4L, precluster=32L, ncores=16L,
  psf_control=list(maxiter=3L),
  batch=list(coupling=list(cycles=2L, initial_cycles=1L, rounds=1L, q=32L)))
experiment$seconds
experiment$parallel
experiment$samples[c("sample_name", "compounds", "GIC")]
experiment$areas
plot(experiment, type="areas")
plot(experiment, sample=1L, type="lines")
plot(experiment, sample=1L, type="profiles", which=1:6)
saveRDS(experiment$channel_groups, file.path(input, "channel_groups.rds"))
} # }
```
