# GC/MS deconvolution for single samples

When two compounds elute together, their ion traces overlap and the
spectrum at any one scan contains contributions from both.
[`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
uses the changing ion patterns across scans to estimate a separate
elution profile, mass spectrum and integrated area for each component.
It takes an ANDI-MS `.cdf` file and writes a NIST-compatible `.msp`
library. Set `xlsx_file` to export the quantification table as well.

The fit separates compound profiles from a smooth background. A hard
log-concavity constraint allows the fitted profiles to be asymmetric or
broad. Final spectra retain every original mass channel.

The [single and multichannel deconvolution
tutorial](https://tomwenseleers.github.io/L0glm/articles/deconvolution.html)
introduces selection with a fixed peak shape. Here the peak shapes are
refined as part of the factorization, which is particularly useful for
overloaded or shark-fin-shaped peaks.

## Two experimental chromatograms

The package includes a cuticular-hydrocarbon window with six co-eluting
compounds (`GCMS_data1`, 119 scans × 671 mass channels) and a
beer-headspace window with a prominent broad acid peak (`GCMS_data2`,
1,764 scans × 518 mass channels). Their row names give retention times
in minutes; column names give nominal m/z. Both examples keep all mass
channels.

``` r

library(L0glm)
library(L0nmf)
library(GCMSprime)
data(GCMS_data1)
data(GCMS_data2)

out <- file.path(tempdir(), "GCMSprime-examples")
dir.create(out, showWarnings = FALSE)
write_GCMS_cdf(GCMS_data1, file.path(out, "GCMS_data1.cdf"))
write_GCMS_cdf(GCMS_data2, file.path(out, "GCMS_data2.cdf"))
```

[`write_GCMS_cdf()`](https://tomwenseleers.github.io/GCMSprime/reference/write_GCMS_cdf.md)
preserves the intensity matrix, including silent ions, and writes scan
times in seconds for ANDI-MS readers. The same
[`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
call then works with a matrix exported from R or a chromatogram exported
by instrument software:

``` r

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

fit1$compounds[c("retention_time", "area", "n_ions")]
GIC(fit2)
coef(fit2, normalize = "l1")       # relative intensities in all 518 channels
fit2$seconds                      # elapsed time by stage
```

These short windows use a supplied Gaussian PSF as the initial profile
dictionary. `maxiter=0L` keeps this PSF fixed;
[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.html)
then refines the individual compound profiles. For a full chromatogram,
the default `nwindows=10L` with positive `maxiter` estimates a smoothly
varying PSF. An alkane standard can supply its initial calibration
through `alkanes_file="alkanes.cdf"`. The hydrocarbon window uses a
constant background; the longer beer window uses the automatically
chosen low-df cubic spline.

### Estimate the PSF automatically

[`estimate_psf()`](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.html)
provides one interface for standard calibration and blind sample
estimation. Both bundled matrices include a full alkane TIC in
`attr(Y, "alkanes_TIC")`; its numeric names give retention times in
minutes. The estimator uses this calibration automatically, resamples to
the sample scan interval and refines the PSF on the sample:

``` r

psf <- estimate_psf(GCMS_data2, nwindows = 3L)
X <- psf$X                         # sparse scan-by-scan PSF design
psf$width_scans

# A constant PSF, calibration alone, or blind fitting without a standard:
constant <- estimate_psf(GCMS_data2, nwindows = 1L)
standard <- estimate_psf(GCMS_data2, method = "standard", nwindows = 3L)
blind <- estimate_psf(GCMS_data2, standard = NULL, method = "blind",
  nwindows = 3L, precluster = 32L)
```

`precluster=FALSE` uses all ion channels, `1L` uses their total, and
`32L` uses at most 32 size-adjusted groups. The same identity-link noise
models are available as in
[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.html),
including negative binomial and Gamma models. The [PSF help
page](https://tomwenseleers.github.io/L0nmf/reference/estimate_psf.html)
describes likelihood fitting, noise scaling and aggregation.

[`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
uses the same estimator. A CDF contains the sample chromatogram; its
accompanying R calibration attribute can be supplied explicitly:

``` r

calibrated <- GCMSdeco(file.path(out, "GCMS_data2.cdf"),
  file.path(out, "GCMS_data2_calibrated.msp"), nwindows = 3L, precluster = 32L,
  psf_control = list(standard = attr(GCMS_data2, "alkanes_TIC")))
```

Calibration peaks should reflect the sample’s chromatographic
conditions. The bundled beer calibration has broader peaks than the PSF
used in the plotted comparison, so the table and figures below use the
explicit Gaussian PSFs in `fit1` and `fit2` above. Automatic calibration
is available through the calls in this section.

The table reports median elapsed time across three fits on the
[benchmark
workstation](https://tomwenseleers.github.io/L0glm/articles/performance.html#reproduce-and-extend).
Each run uses eight OpenMP threads and one BLAS thread. Timings include
CDF reading, PSF initialization, factorization, original-channel
completion and MSP export. All three runs return the same profiles and
spectra.

| Dataset | Scans | Mass channels | Components | Seconds | Quasi-GIC | Adjusted weighted R² |
|:---|---:|---:|---:|---:|---:|---:|
| `GCMS_data1` | 119 | 671 | 6 | 5.14 | 246,473 | 0.98999 |
| `GCMS_data2` | 1,764 | 518 | 38 | 68.46 | 1,467,419 | 0.99498 |

Weighted R² uses weights `1/pmax(fitted_mean, 0.1)` and a separate
weighted mean for each mass channel. Its adjustment is conditional on
the fitted profiles, using the positive spectral and background
coefficients as the parameter count.

The default working-correlation adjustment is estimated from temporal
differences and held fixed during final rank comparisons. Both windows
have an estimated information factor of 1.98. The permissive
initialization retains candidates while their profiles are refined;
final broad-fragment decisions compare merged and equally refined
unmerged local models. The beer fit represents the extended acid
contribution with one log-concave profile, retaining all 518 mass
channels. Diagnostics are available in `fit2$model$contained_fusion` and
`fit2$model$initialization`.

### Cuticular hydrocarbons

![Background-corrected ion traces and fitted compounds in
GCMS_data1](../reference/figures/gcms-data1-traces.png)

Background-corrected ion traces and fitted compounds in GCMS_data1

![Observed, corrected, fitted and Pearson-residual heatmaps for
GCMS_data1](../reference/figures/gcms-data1-heatmaps.png)

Observed, corrected, fitted and Pearson-residual heatmaps for GCMS_data1

![Elution profiles and mass spectra for
GCMS_data1](../reference/figures/gcms-data1-profiles-spectra.png)

Elution profiles and mass spectra for GCMS_data1

The profile/spectrum panels show every recovered compound.

Download
[traces](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data1-traces.pdf),
[heatmaps](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data1-heatmaps.pdf),
[all profiles and
spectra](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data1-profiles-spectra.pdf),
or the [MSP
library](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data1.msp).

### Beer volatiles and a broad acid peak

![Background-corrected ion traces and fitted compounds in
GCMS_data2](../reference/figures/gcms-data2-traces.png)

Background-corrected ion traces and fitted compounds in GCMS_data2

![Observed, corrected, fitted and Pearson-residual heatmaps for
GCMS_data2](../reference/figures/gcms-data2-heatmaps.png)

Observed, corrected, fitted and Pearson-residual heatmaps for GCMS_data2

![Elution profiles and mass spectra for
GCMS_data2](../reference/figures/gcms-data2-profiles-spectra.png)

Elution profiles and mass spectra for GCMS_data2

The profile/spectrum preview includes the broadest component and five
other large components, in temporal order. The PDF contains every
recovered component. The broad acid contribution and neighbouring
components are shown on the original ion channels. The final residual
pass can retain additional small components when their likelihood gain
exceeds the compound penalty.

Download
[traces](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data2-traces.pdf),
[heatmaps](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data2-heatmaps.pdf),
[all profiles and
spectra](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data2-profiles-spectra.pdf),
or the [MSP
library](https://tomwenseleers.github.io/GCMSprime/gcms/gcms-data2.msp).

## Plot and export the results

``` r

plot(fit1, type = "lines", nrows = 1L,
  minutes_per_panel = diff(range(fit1$time)) + 0.01, labels = TRUE)
plot(fit1, type = "heatmap")
plot(fit1, type = "profiles", nrows = 3L)

plot(fit2, type = "lines", file = file.path(out, "GCMS_data2_traces.pdf"),
  nrows = 1L, minutes_per_panel = diff(range(fit2$time)) + 0.01)
plot(fit2, type = "heatmap", file = file.path(out, "GCMS_data2_heatmaps.pdf"))
plot(fit2, type = "profiles", file = file.path(out, "GCMS_data2_profiles_and_spectra.pdf"),
  nrows = 3L)
```

Grey traces show individual observed ion signals after subtracting the
fitted background. Coloured curves show the fitted compounds, scaled by
their strongest ion; the square-root display makes weak peaks easier to
inspect. The four raster heatmaps show the observations,
background-corrected observations, fitted chromatogram and Pearson
residuals. Red vertical lines mark compound apices. Profile-and-spectrum
panels follow increasing apex retention time; `order="area"` displays
the largest compounds first.

Supply `xlsx_file=file.path(out, "GCMS_data2.xlsx")` when fitting to
also export the compound table and integrated areas. Each spectrum
occupies one unwrapped cell, with fragments ordered by intensity. The
returned object stores the source CDF path rather than a second copy of
the observed chromatogram; plotting reloads the required scans from that
file.

## Fit a supplied PSF design directly

Use
[`splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.html)
when a suitable PSF design and background basis are already available.
This also gives direct control over the background for a short
chromatographic window:

``` r

Y <- GCMS_data1
X <- convmatrix(convkernel("gaussian", 65L, list(sigma = 8)), nrow(Y))
B <- matrix(1, nrow(Y), 1L)          # constant background in this short window
nmf <- splogconNMF(X, Y, background = B,
  family = quasipoisson("identity"), precluster = 16L,
  correlation = "auto", ncores = 8L)
nmf$rank
nmf$H                              # unit-sum, log-concave profiles
coef(nmf)                          # full-channel integrated intensities
GIC(nmf)
```

For `GCMS_data2`, a low-dimensional nonnegative spline background can be
constructed with
`splines::bs(as.numeric(rownames(Y)), df=4L, intercept=TRUE)`.
Quasipoisson is useful for instrumental intensities with extra
variation; Poisson, Gaussian, negative binomial and Gamma signal models
are also available. See [the likelihood and shape-constraint
details](https://tomwenseleers.github.io/L0nmf/articles/algorithm.html).

## Understand the method

[GC/MS deconvolution: mathematical
basis](https://tomwenseleers.github.io/GCMSprime/articles/algorithm.md)
explains known and estimated PSFs, compound-dependent profiles, noise,
background and best rank tuning. [L0nmf’s mathematical
article](https://tomwenseleers.github.io/L0nmf/articles/algorithm.html)
gives the objective functions, constraints and pseudocode.

The experimental data are described in [Van Oystaeyen et
al. (2014)](https://doi.org/10.1126/science.1244899) and [Schreurs et
al. (2024)](https://doi.org/10.1038/s41467-024-46346-0).

## Identify the recovered compounds

[Alkane calibration and library
search](https://tomwenseleers.github.io/GCMSprime/articles/library-search.md)
add retention indices, ranked spectral matches and annotations to the
Excel and MSP reports without refitting the chromatogram.
