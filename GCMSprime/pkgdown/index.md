# GCMSprime

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">GC/MS DECONVOLUTION · QUANTIFICATION · IDENTIFICATION</p>
<h1>Accurate integration.<br><span>Cleaner spectra.</span></h1>
<p class="l0-lead">Fast and scalable GC/MS deconvolution and quantification for single samples or large multisample experiments to separate overlapping peaks, with inbuilt spectral library search.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/gcms-deconvolution.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="articles/library-search.html">Search spectral libraries</a></div>
<p class="l0-hero-note">R · Full mass spectra · CDF input · MSP &amp; Excel output</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="GCMSprime hex logo showing resolved chromatographic peaks and their mass spectra" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="articles/gcms-deconvolution.html"><strong>Deconvolve</strong><span>Resolve overlapping compounds and recover their spectra</span></a>
<a href="reference/gcmsDeco.html"><strong>Quantify</strong><span>Integrated areas, fitted profiles and complete reports</span></a>
<a href="articles/library-search.html"><strong>Identify</strong><span>Fast spectral matching with column-specific retention indices</span></a>
</div>

A GC–MS run records mixtures of ions whenever compounds overlap in time.
GCMSprime separates these signals into elution profiles and full mass spectra,
then integrates each resolved compound. It accepts ANDI-MS CDF files and
writes an MSP library of the recovered spectra and an Excel quantification
table when requested. For a sample series, it aligns retention times and
shares spectra across related occurrences.

An alkane ladder can calibrate peak shapes and retention indices. Supplying a
spectral library adds ranked candidate matches, with the retention-index
agreement and the reference source recorded alongside each score.

See the [16-sample benchmark](#benchmarks-against-other-software) for peak-area
and spectral recovery across seven workflows.

## Installation

```r
install.packages("remotes")
remotes::install_github(c(
  "tomwenseleers/L0glm/L0glm",
  "tomwenseleers/logconNMF/logconNMF",
  "tomwenseleers/GCMSprime/GCMSprime"
))
```

For the conditional Poisson and quasi-Poisson solver, install NumPy, SciPy,
Numba and threadpoolctl in one Python environment:

```sh
python -m pip install numpy scipy numba threadpoolctl
```

Set `options(logconNMF.python = "/path/to/python")` if needed.
`logconNMF::conditionalSolverAvailable()` checks the interpreter. The first use
compiles and caches numerical kernels. `algorithm="auto"` selects a compatible
solver. Use `algorithm="conditional"` to require this backend explicitly.

## Automatic deconvolution for single samples

One call estimates the response shape, smooth background, compound profiles,
full mass spectra and integrated areas, then writes NIST-compatible MSP and
Excel reports. An alkane ladder can calibrate the PSF and retention indices.

```r
library(GCMSprime)
fit <- gcmsDeco("sample.cdf", "sample.msp", xlsx_file="sample.xlsx",
  alkanes_file="alkanes.cdf", precluster=32L)
plot(fit, type="lines", file="sample_traces.pdf")
plot(fit, type="heatmap", file="sample_heatmaps.pdf")
plot(fit, type="profiles", file="sample_profiles_spectra.pdf")
```

The compound table gives apex retention times, integrated areas and the
number of fitted ions. Use the trace and heatmap plots to inspect overlapping
peaks and residual structure, and use the profile and spectrum plots to
inspect individual components.

```r
head(fit$compounds[c("compound", "retention_time", "area", "n_ions")])
fit$seconds
spectra <- coef(fit, normalize="l1")  # unit-L1 spectra, shared 0.1% ion cutoff
```

The fit allows asymmetric and broad profiles under a hard log-concavity
constraint. The rank criterion selects the component set, and final
refinement uses every original mass channel.

### Six co-eluting hydrocarbons

The bundled cuticular-hydrocarbon window illustrates recovered elution profiles
and spectra. This short-window example uses a supplied PSF dictionary and a
constant background; individual compound shapes are then refined.

```r
data(gcmsData1, package="logconNMF")
writeGCMSCDF(gcmsData1, "GCMS_data1.cdf")
example_fit <- gcmsDeco("GCMS_data1.cdf", "GCMS_data1.msp",
  nwindows=1L, baseline_df=1L, precluster=7L, correlation="auto", ncores=8L,
  control=list(initial_relaxed=FALSE, fusion_passes=0L),
  psf_control=list(kernel_est=logconNMF::convKernel("gaussian", 65L,
    list(sigma=9)), maxiter=0L))
```

Among six-factor solutions with chemically plausible retention times, a
sensitivity grid over PSF widths and precluster counts selected this fit for
its combination of low GIC, accurate apex times, compact profiles and strong
enrichment of the compounds' diagnostic ions.

[![Grey ion traces, six resolved profiles and matching spectra for six co-eluting hydrocarbons](reference/figures/gcms-data1-decomposition.png)](articles/gcms-deconvolution.html#cuticular-hydrocarbons)

### Beer volatiles and a broad acid peak

A beer-headspace window contains a strongly asymmetric acid peak mixed with
narrower neighbours. Flexible log-concave profiles accommodate this shape while
the spectra use all 518 channels. The same workflow represents the broad acid
contribution as one log-concave profile among the narrower peaks. The
single-sample tutorial reports the measured CDF-to-MSP times and fit metrics.

[![Beer chromatogram with fitted compounds and the broad acid contribution](reference/figures/gcms-data2-traces.png)](articles/gcms-deconvolution.html#beer-volatiles-and-a-broad-acid-peak)

The [single-sample tutorial](articles/gcms-deconvolution.html) includes both
datasets, fitted-signal heatmaps, profile/spectrum panels, timings and downloadable
MSP/PDF outputs. [GC/MS deconvolution: mathematical basis](articles/algorithm.html)
explains the progression from known PSFs to estimated and component-specific shapes.

## Automatic deconvolution for multiple samples

Point to a directory containing sample-group folders and a `standards/`
subdirectory with alkane calibrations. Acquisition timestamps determine run
order. The simultaneous model aligns compounds and shares their spectra while
retaining sample-specific profiles, backgrounds and areas.

```r
experiment <- gcmsDeco("chromatograms", output_dir="experiment_results",
  precluster=32L, batch=list(plots=TRUE))
experiment$samples
experiment$areas
```

Grouped PSF regressions divide the core budget between independent samples and native threads. The scheduler estimates how many regressions fit in available memory, including Windows commit headroom. After those workers exit, a single process shares the original ion counts and fitted means. It updates profiles and backgrounds across samples and fits disjoint ion columns of the shared spectra concurrently. Peak merges and exchanges retain their original order. The classic solver can stream original-channel blocks from disk.
The full-size fit recovered **7,917 of 7,917 occurrences** with **16 false calls**
in **11.24 min** across 16 chromatograms of 10,000 scans and 550 channels.
The [two 16-sample experiments](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html#benchmark-scales) give the measured resource use, timings and recovery at both acquisition sizes.

In the [small-scale 16-sample software comparison](articles/benchmarks.html), GCMSprime recovers
**767 of 769 true occurrences in 70.81 s**, including the acid compound in
all seven samples where it occurs. Each chromatogram has 100 candidate
compounds, 1,200 scans and 258 mass channels, with drift, concentration-dependent
broadening and a smooth background derived from experimental GC/MS data.
Initialization uses 32 ion groups and the final fit uses every original ion.
The grouped regressions used four two-thread workers. The joint fit reused the eight-thread budget on the original ion channels.

## Benchmarks against other software

On **16 chromatograms with 100 reference compounds**, GCMSprime recovers
**767 of 769 true peak occurrences**, with **42 false positives**
(peak F1 **0.972**). Its total-ion peak areas show **robust adjusted R² 0.999995**
and **median absolute area error 0.10%**.

[![True and estimated total-ion peak areas with robust fits](reference/figures/simulated-16-area-predictions.png)](reference/figures/simulated-16-area-predictions.png)

All methods use total-ion count equivalents, accounting for native model-ion
fractions and integration-time units. The robust log–log fit and prediction
band show the main error distribution. Points have opacity 0.5, so darker areas
show denser overlap. The dashed 1:1 reference is visible in every panel.

[![Reciprocal spectral cosine similarities for seven workflows](reference/figures/simulated-16-spectral-violin.png)](reference/figures/simulated-16-spectral-violin.png)

[![Within-spectrum F1 using the shared 0.1% base-peak cutoff](reference/figures/simulated-16-spectral-support.png)](reference/figures/simulated-16-spectral-support.png)

The shared **0.1% base-peak cutoff** gives GCMSprime **mean ion-support
F1 0.995**, compared with **0.860** for the strongest comparator. This is
also the default cutoff for GCMSprime's reported spectra and spectral exports.
Matching requires reciprocal best spectral hits within **30 RI units**.
The reference includes compounds absent from individual samples.

See [benchmarks against other software](articles/benchmarks.html) for timings,
TP/FP counts, square-root signal heatmaps, recovered compound areas and
regular/acid elution profiles. The
[GCMSbenchmark package](https://github.com/tomwenseleers/GCMSprime/tree/master/GCMSbenchmark)
provides the seeded simulation and reproducible comparison.

## Spectral library search

**Open and commercial libraries can be imported.** MoNA and Golm MSP collections,
GCMS-ID CSV predictions, and OpenChrom-exported Agilent or Shimadzu libraries
are supported. Large catalogues can remain in indexed SQLite storage.

```r
lib <- readMSP("library.msp", backend="sqlite", db_file="library.sqlite")
identified <- gcmsDeco("sample.cdf", "identified.msp",
  xlsx_file="identified.xlsx", alkanes_file="alkanes.cdf", library=lib,
  search=list(phase="DB5", ri_tolerance=30, n_matches=5,
              method=c("stein_scott", "hellinger")))
```

The compounds recovered from single- or multisample deconvolution can be
automatically annotated using **spectral match quality and retention index**.
Alkane-derived indices first restrict candidates to the selected column class
and tolerance; fast native search then ranks their spectra. Multisample
consensus IDs support a unique-compound list across the experiment, with
sample-specific RI evidence and full match provenance retained.

Stein–Scott, cosine/angular, Hellinger and symmetric probability/count
comparisons are available, with optional ion-rarity weighting. Excel and MSP
reports include the ranked matches and annotations. The
[library-import tutorial](articles/spectral-libraries.html) covers formats and
storage; [retention indices and library search](articles/library-search.html)
shows single- and multisample identification. The
[identification theory](articles/identification-theory.html) gives the formulas,
calibration and comparison benchmarks.


## Citation

An article describing the methods implemented in GCMSprime is in preparation.
If you use the package before the article is formally published, we would
appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026).
*GCMSprime: Scalable and Automated GC–MS Deconvolution, Quantification, and Spectral Library Search in R*.
R package version 0.1.0.9000,
<https://tomwenseleers.github.io/GCMSprime/>.

The [Authors and Citation page](authors.html#citation) also provides the
citation in BibTeX format. Please also cite the methodological sources linked
in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize and
refactor code and help draft the documentation. The key ideas and methodological
decisions are the authors' own, and the authors take full responsibility for the
work presented here.
