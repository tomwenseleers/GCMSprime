# Multisample GC/MS deconvolution, alignment and shared spectral recovery

Supply a directory to
[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
to process an experiment. Sample group subdirectories are retained, and
`standards/` contains alkane calibrations. The closest standard in
acquisition time calibrates each sample. Results share spectra through
[`logconNMF::logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.html)
while keeping original scan grids.

## Usage

``` r
# S3 method for class 'GCMSbatch'
fitted(object, sample = 1L, ...)

# S3 method for class 'GCMSbatch'
coef(object, sample = NULL, normalize = c("none", "l1"), ...)

# S3 method for class 'GCMSbatch'
plot(x, sample = NULL, ...)

# S3 method for class 'GCMSbatch'
logLik(object, ...)

# S3 method for class 'GCMSbatch'
GIC(object, ...)

# S3 method for class 'GCMSbatch'
AIC(object, ..., k = 2)

# S3 method for class 'GCMSbatch'
BIC(object, ...)

# S3 method for class 'GCMSbatch'
EBIC(object, ...)

# S3 method for class 'GCMSbatch'
MBIC(object, ...)

# S3 method for class 'GCMSbatch'
nobs(object, ...)

# S3 method for class 'GCMSbatch'
print(x, ...)
```

## Arguments

- sample:

  One sample index or identifier. `NULL` in
  [`coef()`](https://rdrr.io/r/stats/coef.html) returns the common
  spectral library. `NULL` in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) plots the
  coupled area table.

- ...:

  Arguments passed to the corresponding sample or coupled method.

- normalize:

  Spectrum normalization: unnormalized or unit L1.

- x, object:

  A multisample result returned by
  [`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md).

- k:

  Information-criterion multiplier.

## Value

A `GCMSbatch` object with `samples`, optional `coupled` fit, sample RDS
`paths`, consensus `areas`, reusable `channel_groups`, the `parallel`
thread allocation, output files and stage timings. Use
`readRDS(x$paths[i])` for a full per-sample
[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
object.

## Details

Directory layout can be `experiment/standards/alkanes.cdf`,
`experiment/beer1/sample001.cdf`, and `experiment/beer2/sample002.cdf`.
Analytical samples are sorted by acquisition timestamp, including its
time zone. The exported CDF timestamp is never used as an acquisition
date. If acquisition dates are absent, pass explicit relative sample
paths in `batch=list(sample_order=...)`.

`batch` supports `method="cosine"` or `"hellinger"`, `align=TRUE`,
`coupling=list()` for
[`logconNMF::logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.html)
controls, `sample_order=NULL`, `exposure=NULL`, `resume=FALSE`,
`plots=FALSE`, `workers=NULL` and `cache=TRUE`. `workers` limits
concurrent samples within the total `ncores` budget. `cache` controls
the aggregated-observation RAM cache during coupled fitting. `exposure`
may be a positive scalar, a list with one scan-wise vector per
acquisition-ordered sample, or a function of a sample's metadata and
time vector. The default preserves vendor intensities. Known integration
durations are appropriate when the data are integrated ion counts.
Sequential full-scan mass spectrometers allocate acquisition time over
m/z; exported peak counts and retained traces do not identify the
per-mass dwell time. No automatic division by the number of exported
channels is applied.

Samples run concurrently with their own sparse PSF and smooth
background. A common channel map is learned from summed, time-binned
chromatograms; pass a saved
[`logconNMF::channelGroups()`](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.html)
map as `precluster` to reuse reference groups. Matching uses the shared
aggregate space. Original-channel means and the likelihood are used for
coupled updates. Shared spectra are exported on the union of recorded
mass ranges; unmeasured channels receive zero likelihood weight. Missing
compounds retain `NA` areas. No gap filling is performed.
`coupled=FALSE` keeps the independently fitted results.

The output directory contains per-sample MSP, Excel and fitted RDS
files, an acquisition-ordered `samples.csv`, `areas.csv`, a consensus
MSP library, and `experiment.xlsx` with samples, compounds, areas and
packed spectra. `batch=list(plots=TRUE)` also writes each sample's trace
and heatmap PDFs. Observations and factors are stored per sample,
enabling streamed fitting; only aggregate observations are cached across
samples in RAM, and the returned object contains paths, tables and
compact common spectra. Original-channel spectral/background updates
read blocks of at most 32 channels. Each active independent fit has its
own temporary working buffers; these are released before the next
assigned sample. A 10,000-by-32 cache needs 2.56 MB per sample, or 2.56
GB for 1,000 samples, in addition to worker and factor memory. Keep the
output directory with the saved object. Explicit `resume=TRUE` reuses
independent fits only when the file metadata and fitting arguments
match; coupled refinement is recomputed.

## See also

[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md),
[`logconNMF::logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.html),
[`writeGCMSCDF()`](https://tomwenseleers.github.io/GCMSprime/reference/writeGCMSCDF.md)

## Examples

``` r
if (FALSE) { # \dontrun{
experiment <- gcmsDeco("chromatograms", output_dir="results",
  family=quasipoisson("identity"), precluster=32L,
  batch=list(method="hellinger",plots=TRUE))
experiment$areas
plot(experiment,sample=1,type="lines")
} # }
```
