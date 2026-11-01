# GCMSprime

GC/MS DECONVOLUTION · QUANTIFICATION · IDENTIFICATION

More compounds.\
Cleaner spectra.
================

State-of-the-art GC/MS deconvolution, quantification and library search.

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

GCMSprime is part of a focused three-package toolkit:
[L0glm](https://tomwenseleers.github.io/L0glm/) supplies fast sparse GLM
selection, [L0nmf](https://tomwenseleers.github.io/L0nmf/) supplies
sparse factorization and PSF estimation, and
[GCMSprime](https://tomwenseleers.github.io/GCMSprime/) turns those
methods into a chromatographic workflow.

## Installation

``` r

install.packages("remotes")
remotes::install_github("tomwenseleers/GCMSprime/GCMSprime")
```

## Fully automatic deconvolution

``` r

library(GCMSprime)
fit <- GCMSdeco("sample.cdf", "sample.msp", xlsx_file="sample.xlsx",
                alkanes_file="alkanes.cdf", precluster=32L)
fit
plot(fit, type="lines", file="sample_traces.pdf")
plot(fit, type="heatmap", file="sample_heatmaps.pdf")
plot(fit, type="profiles", file="sample_profiles_spectra.pdf")
```

The workflow calibrates a time-varying PSF, models a smooth baseline,
resolves overlapping compounds and refines their spectra on every
original mass channel. Results are ordered by apex retention time. MSP
files are readable by NIST MS Search; Excel reports contain integrated
areas and compact spectra.

## Experimental examples

The [GC/MS
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)
includes two experimental chromatogram windows, complete fitted signals,
spectra, timings and exports.

![Co-eluting hydrocarbon
profiles](reference/figures/gcms-data1-traces.png)

Co-eluting hydrocarbon profiles

![Beer volatile profiles](reference/figures/gcms-data2-traces.png)

Beer volatile profiles

The underlying [factorization and PSF
methods](https://tomwenseleers.github.io/L0nmf/) are available
independently in L0nmf.

## Spectral libraries

Import MSP libraries into sparse memory or SQLite storage, preserving
raw peaks, chemical identifiers, synonyms and retention-index
provenance:

``` r

lib <- read_msp("library.msp", backend="sqlite", db_file="library.sqlite")
spectra(lib, which=1:1000, nominal=TRUE, normalize="l2")
library_metadata(lib, which=1:10)
retention_indices(lib, which=1:1000, phase="DB5")
```

[Library formats and
storage](https://tomwenseleers.github.io/GCMSprime/articles/spectral-libraries.md)
covers OpenChrom MSP exports, MoNA/Golm MSP and GCMS-ID CSV predictions.

## Identify recovered compounds

Calibrate retention indices from a diagnostic-ion alkane ladder and
combine column-specific RI constraints with fast spectral matching:

``` r

cal <- calibrate_ri("alkanes.cdf", carbon=7:23)
identified <- identify_compounds(fit, lib, calibration=cal,
  phase="DB5", ri_tolerance=30, n_matches=5,
  method=c("stein_scott", "hellinger"),
  msp_file="sample_identified.msp", xlsx_file="sample_identified.xlsx")
```

For metric spectral distances, use `method="angular"` or
`"stein_scott_angular"`. Reusable ion-rarity weights come from
`library_ion_weights(lib)` and are passed through `ion_weights=`.
Multinomial matching balances unknown reference counts to the sample
total automatically.

[RI calibration and library
search](https://tomwenseleers.github.io/GCMSprime/articles/library-search.html)
compares five similarity measures and explains the peak-onset timing
convention, retained ion totals and scalable SQLite search.

## Citation

An article describing the methods implemented in GCMSprime is in
preparation. If you use the package before the article is formally
published, we would appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026). *GCMSprime: Sparse GC/MS Deconvolution
and Spectral Libraries*. R package version 0.1.0,
<https://tomwenseleers.github.io/GCMSprime/>.

The [Authors and Citation
page](https://tomwenseleers.github.io/GCMSprime/authors.html#citation)
also provides the citation in BibTeX format. Please also cite the
methodological sources linked in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize
and refactor code and help draft the documentation. The key ideas and
methodological decisions are the authors’ own, and the authors take full
responsibility for the work presented here.
