# GC/MS deconvolution for multiple samples

[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
accepts an experiment directory and analyzes its chromatograms together.
The batch can share PSF calibration. Each sample has its own native scan
grid, smooth background and log-concave compound profiles. Simultaneous
factorisation shares mass spectra across samples and returns one
consensus compound library with a sample-by-compound area table.

## Organize and fit an experiment

``` text
chromatograms/
  standards/
    alkanes.cdf
  beer1/
    injection001.cdf
    injection002.cdf
  beer2/
    injection003.cdf
```

``` r

library(GCMSprime)

experiment <- gcmsDeco("chromatograms", output_dir="beer_results",
  family=quasipoisson("identity"), precluster=32L,
  batch=list(plots=TRUE))

experiment$samples
experiment$areas
coef(experiment, normalize="l1")
GIC(experiment)
plot(experiment, sample=1L, type="lines", file="sample1_traces.pdf")
```

The `standards` subdirectory is reserved for calibration chromatograms.
Analytical samples are sorted using their acquisition timestamps,
including time-zone information. Group labels come from their relative
directories. Each sample uses the nearest alkane standard in acquisition
time. A single standard can calibrate the whole experiment;
`alkanes_file=` explicitly selects a particular calibration.

For files without acquisition timestamps, supply their relative paths in
the intended run order:

``` r

experiment <- gcmsDeco("chromatograms", output_dir="beer_results",
  batch=list(sample_order=c("beer1/injection001.cdf",
                            "beer1/injection002.cdf",
                            "beer2/injection003.cdf")))
```

## Two distinct 16-sample benchmarks

The two experiments use different simulated chromatograms. Both retain
the native scan grids and report recovery against the generating
compounds.

| Setting | Small-scale software comparison | Large-scale resource validation |
|:---|---:|---:|
| Samples | 16 | 16 |
| Scans per sample | 1,200 | 10,000 |
| Recorded ion channels | 258 | 550 |
| Reference compounds | 100 | 1,000 |
| True occurrences | 769 | 7,917 |
| Samples containing the acid | 7 | 8 |
| PSF calibration | Blind, two windows | Standard-assisted, ten windows |
| Complete 16-sample analysis | **124.70 s (2.08 min)** | **1,360.13 s (22.67 min)** |

Times cover the full batch and exclude simulation. In the small
experiment, calibration and grouping took 15.13 s, grouped
initialization and joint refinement took 103.90 s, and result
construction and export took 5.41 s. The benchmark wrapper added 0.26 s.
The large experiment is documented below.

## Small-scale 16-sample software comparison

The small-scale fit (seed `9118`) uses two blind-PSF windows, three PSF
iterations, 32 common ion groups and up to 12 initial conditional sweeps
on the original channels. All profile figures for this simulation on
both websites use the same saved fit.

The **16 chromatograms** share a library of **100 reference compounds**.
Each spans a 20-minute retention window, with 1,200 scans and all 258
mass channels. Compounds are absent from about half the samples, giving
**769 true occurrences**. Library retention indices determine their
positions. One compound has a broad shark-fin-shaped acid profile; the
remaining profiles are exponentially modified Gaussians. Higher
concentrations extend peaks towards later times from approximately fixed
leading edges, and a small positive drift shifts successive injections
to the right.

The smooth B-spline background is derived from the experimental
`gcmsData2` window and contributes 3.03% of the integrated compound
signal. The simulation retains its native abundance scale and uses
quasi-Poisson noise with **dispersion 783.18**, estimated from that
window. GCMSprime estimates dispersion separately in each sample, with a
**median estimate of 829.55**. The generating retention-index maps are
used for scoring; the deconvolution receives neither the known areas nor
the true profiles.

The reproducible recipe includes the empirical background and complete
reference library:

``` r

library(GCMSbenchmark)
bench <- system.file("benchmarks", package="GCMSbenchmark", mustWork=TRUE)
source(file.path(bench, "simulated-16.R"), encoding="UTF-8")
small <- runBenchmark16("benchmark-small", method="GCMSprime", ncores=8L)
small$summary[c("analysis_s", "TP", "FP", "FN", "F1")]
# Optional figures from the saved results:
source(file.path(bench, "plot-simulated-16.R"), encoding="UTF-8")
plotBenchmark16("benchmark-small")
```

The analysis recovered **767 of 769 true occurrences**, with **42 false
positives** (peak F1 **0.972**), including all seven acid occurrences.
Robust adjusted log-area R² was **0.999995**, the median absolute
peak-area error was **0.10%**, and mean spectral ion-support F1 was
**0.995** at the shared 0.1% base-peak threshold.

The CDF-to-results analysis took **124.70 s** with four two-thread
sample workers on an HP ZBook Fury G8, Intel Core i9-11950H, eight
physical cores, 16 logical threads, Windows and R 4.5.0. CDF generation
is separate.

[![Total-ion area accuracy with robust log–log
fits](../reference/figures/simulated-16-area-predictions.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-area-predictions.png)

[![Acid profiles across all samples, truth dashed and estimates
solid](../reference/figures/simulated-16-acid-profiles.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-acid-profiles.png)

The [benchmarks against other
software](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.md)
give all seven workflows, matching rules, noise calibration, profile
plots and tables.

## Large-scale 16-sample resource benchmark

This independent simulation (seed `171117`) has **16 chromatograms,
10,000 scans, 550 recorded ion channels and 1,000 reference compounds**.
Their spectra occupy 258 channels. The remaining channels contain
background. The simulated peaks broaden with abundance and drift in
retention time. The data include irregular scan intervals, a background
estimated from `gcmsData2` and quasi-Poisson dispersion 783.18. A
separate narrow-peak standard supplies the time-varying PSF.

The conditional solver estimated **1,005 shared identities** from the
original ion measurements. All **7,917 true occurrences** were
recovered, including every occurrence of the shark-fin acid. The largest
absolute area error among the eight acid occurrences was **0.19%**.

| Quantity                          | Measured result |
|:----------------------------------|----------------:|
| True occurrences recovered        |   7,917 / 7,917 |
| False positives / false negatives |          16 / 0 |
| Occurrence F1                     |         0.99899 |
| Joint quasi-GIC                   |      34,070,769 |
| Median spectral cosine            |       0.9999989 |
| Mean within-spectrum ion F1       |         0.99516 |
| Mean profile cosine               |        0.999983 |
| Median absolute peak-area error   |          0.095% |
| Acid occurrences recovered        |           8 / 8 |

Scoring uses maximum-cardinality one-to-one matching within 0.25 minutes
and a minimum spectral cosine of 0.9. Ion-support F1 uses a shared
cutoff of 0.1% of base-peak intensity. Area and profile accuracy are
measured on matched occurrences. The small-scale software comparison
uses RI-constrained reciprocal hits.

| CDF workflow stage | Seconds |
|:---|---:|
| Reading data, preparing sample PSFs and initialization | 226.10 |
| Grouped PSF regression, alignment and joint original-channel refinement | 1092.54 |
| Constructing and exporting results | 41.49 |
| **Complete fit and exports** | **1,360.13** |

Of the 1,092.54 s in the second stage, 116.03 s were spent on grouped
PSF regressions and 976.14 s on the joint conditional fit. Simulation,
separate standard calibration and scoring against truth are outside
these timings. With 16 logical threads requested, the memory-aware
scheduler chose **4 workers with 4 threads each** for independent fits.
The joint updates use one BLAS thread. On the Windows HP ZBook Fury G8
(Core i9-11950H, eight physical cores, 64 GB RAM), the combined process
tree reached **7.01 GiB private memory** and **7.43 GiB working set**.
At least 28.00 GiB of physical memory remained available. Mean CPU use
was 2.17 logical cores.

![Measured process-tree memory and CPU use for sixteen full-sized
chromatograms](../reference/figures/fullsize-16-resources.png)

Measured process-tree memory and CPU use for sixteen full-sized
chromatograms

![Occurrence recovery and total-ion peak-area accuracy across sixteen
full-sized samples](../reference/figures/fullsize-16-recovery.png)

Occurrence recovery and total-ion peak-area accuracy across sixteen
full-sized samples

![Acid profiles: inferred solid and true dashed, blue for smaller areas
and red for larger areas](../reference/figures/fullsize-16-acid.png)

Acid profiles: inferred solid and true dashed, blue for smaller areas
and red for larger areas

The acid panel displays the **eight correctly matched acid-positive
samples from this separate 16-sample large-scale simulation**. The acid
is absent from the other eight samples, which are omitted from this
positive-case shape comparison; this panel does not assess false acid
assignments in absent samples. Colours rank the eight displayed true
areas from blue to red. The x-axis is seconds relative to each true
1%-height onset. In contrast, the small-scale figures include all 16
samples (seven acid-positive, nine with zero true acid signal), rank
colours over all 16 areas, and use retention index. Both displays retain
original count amplitudes without unit-height or unit-area
normalization. Fitting uses each sample’s original scan grid.

The [complete reproducible
script](https://github.com/tomwenseleers/GCMSprime/blob/master/benchmarks/benchmark_fullsize_16.R)
and [compact measured
results](https://github.com/tomwenseleers/GCMSprime/tree/master/benchmarks/results/fullsize16)
are supplied with the source. From the GCMSprime repository root:

``` r

# Run from a GCMSprime repository checkout; use a fresh output directory.
source("benchmarks/benchmark_fullsize_16.R")
benchmarkFullsize16("benchmark-large", stage="all", ncores=16L)
plotFullsize16Acid("benchmark-large")
```

This writes inputs, fitted results (including timings), scoring tables
and the acid PNG to `benchmark-large`. Separate `stage="prepare"`,
`"fit"` and `"score"` calls are available when restarting a stage. The
default directory, when omitted, is `codex_sandbox/fullsize16`.

## Outputs and memory use

The output directory contains:

- `samples.csv`: acquisition order, group, calibration file and sample
  result paths.
- `areas.csv`: integrated compound areas; undetected occurrences have
  `NA`.
- `consensus.msp`: the common spectral library, using the shared 0.1%
  base-peak cutoff.
- `experiment.xlsx`: sample information, consensus compounds, packed
  spectra and areas.
- `experiment.rds`: the experiment object.
- `samples/<sample_id>/`: fitted sample RDS, NIST-compatible MSP and
  Excel results, plus trace and heatmap PDFs when `batch$plots=TRUE`.

Spectral text cells are unwrapped and list fragments by decreasing
intensity. Sample objects retain factor matrices and CDF paths. The
grouped PSF regressions run in parallel within the `ncores` budget. The
count solver limits automatic worker concurrency using the scan count
and available memory, including Windows commit headroom. Each worker
releases its design before joint refinement. The conditional updates
hold original-ion observations and fitted means in memory and use the
native time grid of each sample. The classic solver can stream
original-channel blocks from disk. Use `batch=list(workers=4L)` to set
concurrency explicitly.

The same ion groups are used across samples for matching. They are
learned from a pooled, time-binned chromatogram or supplied through
`precluster` as a saved [channel
map](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.html),
for example one trained on reference spectra. Each sample’s PSF and scan
grid remain separate.

``` r

sample1 <- readRDS(experiment$paths[1L])
sample1$compounds
plot(sample1, type="profiles", file="sample1_profiles_spectra.pdf")
fitted(experiment, sample=1L)
```

The [simultaneous NMF
tutorial](https://tomwenseleers.github.io/logconNMF/articles/simultaneous-nmf.html)
gives the mathematical model, alignment, spectral grouping, local
updates and simulation results. Its validation includes three samples
with all 258 mass channels, differing scan grids,
concentration-dependent tails and compounds absent from half the
samples. All seven signal noise families are supported.

## Scan spacing and exposure

Observed matrices stay on their native scan grids. Monotone alignment
changes profile coordinates and bin widths; its Jacobian preserves
integrated peak areas. The model is evaluated directly at the original
scan times.

By default, GCMSprime models the vendor’s exported intensities with unit
exposure. Quadrupole acquisition distributes scan time across the mass
range. The effective time per mass also depends on scan step size, scan
averaging and acquisition settings; see the [Agilent method
documentation](https://openlab.help.agilent.com/en/mergedProjects/GCMSHelp/MEAcquiUI.htm).
ANDI files can contain arbitrary abundance units, so scan duration and
the number of exported mass peaks alone do not identify physical count
exposure.

If integrated counts and their durations are known, supply exposure
explicitly:

``` r

experiment <- gcmsDeco("chromatograms", output_dir="count_results",
  family=poisson("identity"),
  batch=list(exposure=function(metadata, time) metadata$scan_duration))
```

Here `scan_duration` is appropriate only when it represents the
integration duration of each reported count. A calibrated per-mass
duration can instead be returned by the function. Metadata vectors are
restricted to the fitted scans. With an identity link, exposure
multiplies the latent mean. It is incorporated into profile and
background design values during fitting.

## Calibrated beer example

An alkane standard supplies narrow peak shapes for both full
chromatograms. The conditional count workflow estimates a pooled
time-varying PSF and transfers it to each sample’s observed scan times.
In this example, `calibration_range=c(6, Inf)` excludes the solvent and
early contaminated region from PSF calibration; the analytical samples
retain their full time range and all 518 measured mass channels.

``` r

standard <- "chromatograms/standards/alkanes.cdf"
experiment <- gcmsDeco("chromatograms", output_dir="beer_results",
  alkanes_file=standard, calibration_range=c(6, Inf), precluster=32L,
  algorithm="conditional", ncores=8L, batch=list(workers=2L))

# Use the calibrated narrow-peak width for onset-based RI coordinates.
reference <- readRDS(experiment$paths[
  grepl("54040014_a", experiment$paths)][[1L]])
cal <- calibrateRI(standard, carbon=7:23, psf=reference$psf_standard)
experiment <- annotateRI(experiment, cal, phase="DB5",
  msp_file="beer_results/consensus.msp",
  xlsx_file="beer_results/experiment.xlsx")
```

A 26-column spline basis with ion-specific nonnegative coefficients
allows the background spectrum to change through each chromatogram. The
fit used four threads per worker and retained **643 shared identities
and 1,079 occurrences**.

| Sample         | Occurrences | Dispersion | Sample quasi-GIC | Adjusted weighted R² |
|:---------------|------------:|-----------:|-----------------:|---------------------:|
| 502131900304_a |         550 |     631.78 |       15,217,743 |              0.98301 |
| 54040014_a     |         529 |     617.21 |       15,419,546 |              0.98316 |

The joint quasi-GIC was **29,721,269**, with each shared spectrum
charged once. Deconvolution took **643.56 s** from CDF input through
export. Calibration and grouping took 138.31 s, grouped regressions and
joint refinement 495.80 s, and result construction and export 8.72 s.
Subsequent RI calibration took 4.05 s and RI annotation and export took
4.85 s. These timings exclude plots.

With identical Agilent NIST records, a ±30 RI window and cosine ≥0.85,
the two beers gave **213 and 203 distinct reciprocal identifications**.
The classic fit gave 163 and 147. Additional forward assignments to an
identity already called in the same sample numbered 87 and 73, compared
with 65 and 54 for the classic fit. These library hits provide external
validation for real samples of unknown composition.

Diagnostic ions identify all 17 alkanes from C7 to C23. The retention
coordinate is the 1%-height onset plus half the local narrow-peak width
at that height. The natural cubic RI curve passes through all accepted
alkane anchors exactly. MSP files contain phase-labelled DB5 indices;
the consensus library uses their medians. Values outside the accepted
calibration range remain unassigned.

The full beer fits still contain residual fragmentation and occasional
overly broad profiles in difficult overlapping regions. Component counts
are estimates under the fitted model, rather than a verified inventory
of chemical compounds.

The following panels zoom into the 4.5–6 minute region of these
full-chromatogram fits. They are not separately fitted cropped windows.
Grey lines show observed ions after subtracting the jointly modelled
background; coloured lines are the fitted compounds.

![Local full-fit profiles for beer sample
502131900304_a](../reference/figures/beer-502131900304_a-local.png)

Local full-fit profiles for beer sample 502131900304_a

![Local full-fit profiles for beer sample
54040014_a](../reference/figures/beer-54040014_a-local.png)

Local full-fit profiles for beer sample 54040014_a

The fitted components can be inspected sample by sample, and common
spectra can be searched against a library with each sample’s
alkane-derived retention indices. Pass `library=` and `search=` to
[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
to include those identifications in per-sample MSP and Excel exports.

Set `coupled=FALSE` to export independent fits. Set
`batch=list(resume=TRUE)` to reuse independently fitted samples whose
source-file metadata and fitting settings match the cache. Simultaneous
refinement is then computed from those fits.

See
[GCMSbatch](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
for options and [single-sample
deconvolution](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)
for diagnostic plots, profile fitting and spectral recovery.
