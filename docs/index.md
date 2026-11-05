# GCMSprime

GC/MS DECONVOLUTION · QUANTIFICATION · IDENTIFICATION

More compounds.\
Cleaner spectra.
================

Fast and scalable GC/MS deconvolution and quantification for single or
multiple samples, with inbuilt spectral library search.

[Get started
↗](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)[Search
spectral
libraries](https://tomwenseleers.github.io/GCMSprime/articles/library-search.md)

R · Full mass spectra · CDF input · MSP & Excel output

![GCMSprime hex logo showing resolved chromatographic peaks and their
mass spectra](logo-dark.svg)

[**Deconvolve**Resolve overlapping compounds and recover their
spectra](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)
[**Quantify**Integrated areas, fitted profiles and complete
reports](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
[**Identify**Fast spectral matching with column-specific retention
indices](https://tomwenseleers.github.io/GCMSprime/articles/library-search.md)

A GC–MS run records mixtures of ions whenever compounds overlap in time.
GCMSprime separates these signals into elution profiles and full mass
spectra, then integrates each resolved compound. It accepts ANDI-MS CDF
files and writes an MSP library of the recovered spectra and an Excel
quantification table when requested. For a sample series, it aligns
retention times and shares spectra across related occurrences.

An alkane ladder can calibrate peak shapes and retention indices.
Supplying a spectral library adds ranked candidate matches, with the
retention-index agreement and the reference source recorded alongside
each score.

## Installation

``` r

install.packages("remotes")
remotes::install_github("tomwenseleers/GCMSprime/GCMSprime")
```

## Automatic deconvolution for single samples

One call estimates the response shape, smooth background, compound
profiles, full mass spectra and integrated areas, then writes
NIST-compatible MSP and Excel reports. An alkane ladder can calibrate
the PSF and retention indices.

``` r

library(GCMSprime)
fit <- GCMSdeco("sample.cdf", "sample.msp", xlsx_file="sample.xlsx",
  alkanes_file="alkanes.cdf", precluster=32L)
plot(fit, type="lines", file="sample_traces.pdf")
plot(fit, type="heatmap", file="sample_heatmaps.pdf")
plot(fit, type="profiles", file="sample_profiles_spectra.pdf")
```

The compound table gives apex retention times, integrated areas and the
number of fitted ions. Use the trace and heatmap plots to inspect
overlapping peaks and residual structure, and use the profile and
spectrum plots to inspect individual components.

``` r

head(fit$compounds[c("compound", "retention_time", "area", "n_ions")])
fit$seconds
```

The fit allows asymmetric and broad profiles under a hard log-concavity
constraint. The rank criterion selects the component set, and final
refinement uses every original mass channel.

### Six co-eluting hydrocarbons

The bundled cuticular-hydrocarbon window illustrates recovered elution
profiles and spectra. This short-window example uses a supplied PSF
dictionary and a constant background; individual compound shapes are
then refined.

``` r

data(GCMS_data1, package="L0nmf")
write_GCMS_cdf(GCMS_data1, "GCMS_data1.cdf")
example_fit <- GCMSdeco("GCMS_data1.cdf", "GCMS_data1.msp",
  nwindows=1L, baseline_df=1L, precluster=16L, correlation="auto", ncores=8L,
  psf_control=list(kernel_est=L0nmf::convkernel("gaussian", 65L,
    list(sigma=8)), maxiter=0L))
```

[![Recovered elution profiles and spectra for six co-eluting
hydrocarbons](reference/figures/gcms-data1-profiles-spectra.png)](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html#cuticular-hydrocarbons)

### Beer volatiles and a broad acid peak

A beer-headspace window contains a strongly asymmetric acid peak mixed
with narrower neighbours. Flexible log-concave profiles accommodate this
shape while the spectra use all 518 channels. The same workflow
represents the broad acid contribution as one log-concave profile among
the narrower peaks. The single-sample tutorial reports the measured
CDF-to-MSP times and fit metrics.

[![Beer chromatogram with fitted compounds and the broad acid
contribution](reference/figures/gcms-data2-traces.png)](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html#beer-volatiles-and-a-broad-acid-peak)

The [single-sample
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)
includes both datasets, fitted-signal heatmaps, profile/spectrum panels,
timings and downloadable MSP/PDF outputs. [GC/MS deconvolution:
mathematical
basis](https://tomwenseleers.github.io/GCMSprime/articles/algorithm.md)
explains the progression from known PSFs to estimated and
component-specific shapes.

## Automatic deconvolution for multiple samples

Point to a directory containing sample-group folders and a `standards/`
subdirectory with alkane calibrations. Acquisition timestamps determine
run order. The coupled model aligns compounds and shares their spectra
while retaining sample-specific profiles, backgrounds and areas.

``` r

experiment <- GCMSdeco("chromatograms", output_dir="experiment_results",
  precluster=32L, batch=list(plots=TRUE))
experiment$samples
experiment$areas
```

Samples can be fitted in parallel, with each active worker loading its
sample into RAM for the independent fit. During coupled refinement,
original observations and factor checkpoints are stored per sample on
disk. The fit reads blocks of original mass channels to evaluate the
likelihood and update spectra and backgrounds. With preclustering,
grouped traces are cached in RAM by default. Shared spectra and matching
summaries also remain in memory. Original scan grids are retained, and
compounds absent from a sample stay unassigned in its area table.

In the [16-sample
benchmark](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html#sixteen-simulated-chromatograms),
GCMSprime recovers all 48 shared compounds and all 380 sample
occurrences, with no false calls in **63.91 s** with 16 workers. It uses
32 groups for initialization and matching, then fits the spectra on all
258 original ions. The tutorial shows timings by stage, the recovered
areas and an application to two full beer chromatograms.

## Spectral library search

**Open and commercial libraries can be imported.** MoNA and Golm MSP
collections, GCMS-ID CSV predictions, and OpenChrom-exported Agilent or
Shimadzu libraries are supported. Large catalogues can remain in indexed
SQLite storage.

``` r

lib <- read_msp("library.msp", backend="sqlite", db_file="library.sqlite")
identified <- GCMSdeco("sample.cdf", "identified.msp",
  xlsx_file="identified.xlsx", alkanes_file="alkanes.cdf", library=lib,
  search=list(phase="DB5", ri_tolerance=30, n_matches=5,
              method=c("stein_scott", "hellinger")))
```

The compounds recovered from single- or multisample deconvolution can be
automatically annotated using **spectral match quality and retention
index**. Alkane-derived indices first restrict candidates to the
selected column class and tolerance; fast native search then ranks their
spectra. Multisample consensus IDs support a unique-compound list across
the experiment, with sample-specific RI evidence and full match
provenance retained.

Stein–Scott, cosine/angular, Hellinger and symmetric probability/count
comparisons are available, with optional ion-rarity weighting. Excel and
MSP reports include the ranked matches and annotations. The
[library-import
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/spectral-libraries.md)
covers formats and storage; [retention indices and library
search](https://tomwenseleers.github.io/GCMSprime/articles/library-search.md)
shows single- and multisample identification. The [identification
theory](https://tomwenseleers.github.io/GCMSprime/articles/identification-theory.md)
gives the formulas, calibration and comparison benchmarks.

## Citation

An article describing the methods implemented in GCMSprime is in
preparation. If you use the package before the article is formally
published, we would appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026). *GCMSprime: Scalable and Automated
GC–MS Deconvolution, Quantification, and Spectral Library Search in R*.
R package version 0.1.0, <https://tomwenseleers.github.io/GCMSprime/>.

The [Authors and Citation
page](https://tomwenseleers.github.io/GCMSprime/authors.html#citation)
also provides the citation in BibTeX format. Please also cite the
methodological sources linked in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize
and refactor code and help draft the documentation. The key ideas and
methodological decisions are the authors’ own, and the authors take full
responsibility for the work presented here.
