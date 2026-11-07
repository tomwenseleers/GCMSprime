# Simulate related chromatograms with shared compound spectra

Generate exponentially modified Gaussian and optional triangular acid
peaks, smooth backgrounds, irregular scan grids, sample drift and
concentration- dependent peak widths and tails. Missing compounds can be
simulated explicitly.

## Usage

``` r
simulateGCMSSamples(
  samples = 3L,
  n = 600L,
  k = 26L,
  S = NULL,
  preset = c("time_area_acid", "same", "time", "time_area"),
  presence = 1,
  absent = NULL,
  irregular = 0.12,
  drift = 2,
  area = 2e+05,
  seed = 811,
  width_time = 0.5,
  width_area = 0.45,
  tail_time = 1,
  tail_area = 1.5,
  rt_area = 0.25,
  acid_fraction = 0.08,
  family = c("poisson", "quasipoisson", "nb", "gaussian", "Gamma", "nb1", "Gamma_scale"),
  dispersion = 1,
  theta = 20,
  retentionIndex = NULL,
  retentionMap = NULL,
  background = NULL
)
```

## Arguments

- samples:

  Number of samples.

- n:

  Number of scans, recycled across samples.

- k:

  Number of distinct compounds.

- S:

  Optional spectral library with one compound per row.

- preset:

  Peak-shape scenario: identical shapes, time-dependent shapes, time and
  concentration-dependent shapes, or the latter with acid peaks.

- presence:

  Probability that each compound occurs in each sample.

- absent:

  Optional logical compounds-by-samples absence matrix.

- irregular:

  Relative scan-width variation, from zero to 0.5.

- drift:

  Drift per sample, in seconds.

- area:

  Typical integrated compound intensity.

- seed:

  Random-number seed.

- width_time, width_area, tail_time, tail_area, rt_area:

  Effects of time and concentration on peak width, tail and retention
  time.

- acid_fraction:

  Fraction of broad triangular acid profiles.

- family:

  Signal distribution name supported by
  [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.html).

- dispersion, theta:

  Noise parameters.

- retentionIndex:

  Optional retention indices, one per row of `S`. When `retentionMap` is
  supplied, the default is the positional `RI` attribute of the spectral
  library.

- retentionMap:

  Optional calibration data frame with increasing `ri` and `time`
  columns. Times are reference retention times in minutes. Compound
  locations follow this map and sample drift warps the map for each
  injection. Each sample returns its known calibration.

- background:

  Optional nonnegative baseline intensity per unit exposure: a scalar,
  one value per mass channel, or a reference-time-by-channel matrix.
  Matrix rows are interpolated across the reference chromatogram for
  each sample's scan grid. The default uses a random low-dimensional
  spline.

## Value

Simulation truth and a list of samples with observations, original bin
edges and times (minutes), integrated-count means, true spectra,
profiles, areas, presence indicators and acquisition dates. A simulated
observation is integrated over its bin; `raw_exposure` is the original
bin duration in seconds. For
[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.html),
supply it as exposure. With `retentionMap`, each sample includes
`retentionIndex` and an increasing `retentionCalibration` table
incorporating its known drift. The `onset` field records
concentration-adjusted leading positions on the reference time grid; the
sample drift then maps them to observed time.

## Details

In concentration-dependent scenarios, peak width and tail increase with
concentration while the leading 1%-height position stays approximately
fixed. The `rt_area` coefficient optionally adds a small concentration-
dependent shift to that position. Triangular acid profiles broaden
towards later times from a fixed leading edge. With a single acid
profile it is placed in the middle of the chromatogram.

## See also

[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md),
[`logconNMF::simulateGCMSSamples()`](https://tomwenseleers.github.io/logconNMF/reference/simulateGCMSSamples.html)

## Examples

``` r
sim <- simulateGCMSSamples(samples=3L,n=c(180L,190L,200L),k=8L,
  presence=0.5,seed=3)
sim$present
```
