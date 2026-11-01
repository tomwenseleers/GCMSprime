# GCMSprime

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">GC/MS DECONVOLUTION · QUANTIFICATION · IDENTIFICATION</p>
<h1>More compounds.<br><span>Cleaner spectra.</span></h1>
<p class="l0-lead">State-of-the-art GC/MS deconvolution, quantification and library search.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/gcms-deconvolution.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="articles/library-search.html">Search spectral libraries</a></div>
<p class="l0-hero-note">R · Full mass spectra · CDF input · MSP &amp; Excel output</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="GCMSprime hex logo showing resolved chromatographic peaks and their mass spectra" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="articles/gcms-deconvolution.html"><strong>Deconvolve</strong><span>Resolve overlapping compounds and recover their spectra</span></a>
<a href="reference/GCMSdeco.html"><strong>Quantify</strong><span>Integrated areas, fitted profiles and complete reports</span></a>
<a href="articles/library-search.html"><strong>Identify</strong><span>Fast spectral matching with column-specific retention indices</span></a>
</div>

GCMSprime is part of a focused three-package toolkit: [L0glm](https://tomwenseleers.github.io/L0glm/) supplies fast sparse GLM selection, [L0nmf](https://tomwenseleers.github.io/L0nmf/) supplies sparse factorization and PSF estimation, and [GCMSprime](https://tomwenseleers.github.io/GCMSprime/) turns those methods into a chromatographic workflow.

## Installation

```r
install.packages("remotes")
remotes::install_github("tomwenseleers/GCMSprime/GCMSprime")
```

## Fully automatic deconvolution

```r
library(GCMSprime)
fit <- GCMSdeco("sample.cdf", "sample.msp", xlsx_file="sample.xlsx",
                alkanes_file="alkanes.cdf", precluster=32L)
fit
plot(fit, type="lines", file="sample_traces.pdf")
plot(fit, type="heatmap", file="sample_heatmaps.pdf")
plot(fit, type="profiles", file="sample_profiles_spectra.pdf")
```

The workflow calibrates a time-varying PSF, models a smooth baseline, resolves overlapping compounds and refines their spectra on every original mass channel. Results are ordered by apex retention time. MSP files are readable by NIST MS Search; Excel reports contain integrated areas and compact spectra.

## Experimental examples

The [GC/MS tutorial](articles/gcms-deconvolution.html) includes two experimental chromatogram windows, complete fitted signals, spectra, timings and exports.

![Co-eluting hydrocarbon profiles](reference/figures/gcms-data1-traces.png)

![Beer volatile profiles](reference/figures/gcms-data2-traces.png)

The underlying [factorization and PSF methods](https://tomwenseleers.github.io/L0nmf/) are available independently in L0nmf.

## Spectral libraries

Import MSP libraries into sparse memory or SQLite storage, preserving raw peaks, chemical identifiers, synonyms and retention-index provenance:

```r
lib <- read_msp("library.msp", backend="sqlite", db_file="library.sqlite")
spectra(lib, which=1:1000, nominal=TRUE, normalize="l2")
library_metadata(lib, which=1:10)
retention_indices(lib, which=1:1000, phase="DB5")
```

[Library formats and storage](articles/spectral-libraries.html) covers OpenChrom MSP exports, MoNA/Golm MSP and GCMS-ID CSV predictions.

## Identify recovered compounds

Calibrate retention indices from a diagnostic-ion alkane ladder and combine column-specific RI constraints with fast spectral matching:

```r
cal <- calibrate_ri("alkanes.cdf", carbon=7:23)
identified <- identify_compounds(fit, lib, calibration=cal,
  phase="DB5", ri_tolerance=30, n_matches=5,
  method=c("stein_scott", "hellinger"),
  msp_file="sample_identified.msp", xlsx_file="sample_identified.xlsx")
```

For metric spectral distances, use `method="angular"` or `"stein_scott_angular"`.
Reusable ion-rarity weights come from `library_ion_weights(lib)` and are passed through `ion_weights=`.
Multinomial matching balances unknown reference counts to the sample total automatically.

[RI calibration and library search](https://tomwenseleers.github.io/GCMSprime/articles/library-search.html) compares five similarity measures and explains the peak-onset timing convention, retained ion totals and scalable SQLite search.

## Citation and acknowledgements

An article describing these methods is in preparation. Please cite the package and the methodological sources linked in the help pages when using this work. OpenAI Codex/GPT-6 Astra and 6.1 Sol assisted code optimization, refactoring and documentation; the authors developed the key ideas and take responsibility for the work.
