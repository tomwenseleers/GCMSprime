# Package index

## Deconvolution, quantification and reports

- [`residuals(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`coef(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`fitted(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`nobs(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`logLik(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`AIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`BIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`GIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`EBIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`MBIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  [`print(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  : Deconvolve a GC/MS chromatogram and export its compound spectra
- [`fitted(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`coef(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`plot(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`logLik(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`GIC(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`AIC(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`BIC(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`EBIC(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`MBIC(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`nobs(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  [`print(`*`<GCMSbatch>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  : Multisample GC/MS deconvolution, alignment and shared spectral
  recovery
- [`plot(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/plot.GCMSdeco.md)
  : Plot GC/MS signals, compound profiles and recovered spectra
- [`write_GCMS_cdf()`](https://tomwenseleers.github.io/GCMSprime/reference/write_GCMS_cdf.md)
  : Write a chromatogram matrix as an ANDI-MS CDF file

## Import and inspect spectral libraries

- [`read_msp()`](https://tomwenseleers.github.io/GCMSprime/reference/read_msp.md)
  [`open_library()`](https://tomwenseleers.github.io/GCMSprime/reference/read_msp.md)
  [`print(`*`<GCMSlibrary>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/read_msp.md)
  : Import a mass spectral library
- [`read_gcms_id()`](https://tomwenseleers.github.io/GCMSprime/reference/read_gcms_id.md)
  : Import a GCMS-ID predicted spectral library
- [`library_metadata()`](https://tomwenseleers.github.io/GCMSprime/reference/library_metadata.md)
  : Inspect spectral-library metadata
- [`spectra()`](https://tomwenseleers.github.io/GCMSprime/reference/spectra.md)
  : Extract spectral-library peaks or a sparse mass matrix
- [`retention_indices()`](https://tomwenseleers.github.io/GCMSprime/reference/retention_indices.md)
  : Extract retention indices and their provenance

## Calibrate retention indices and identify compounds

- [`calibrate_ri()`](https://tomwenseleers.github.io/GCMSprime/reference/calibrate_ri.md)
  [`predict(`*`<RIcalibration>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/calibrate_ri.md)
  [`print(`*`<RIcalibration>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/calibrate_ri.md)
  : Calibrate retention indices from an n-alkane standard
- [`library_search()`](https://tomwenseleers.github.io/GCMSprime/reference/library_search.md)
  [`as.data.frame(`*`<GCMSmatches>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/library_search.md)
  [`print(`*`<GCMSmatches>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/library_search.md)
  : Identify deconvoluted spectra using retention indices and spectral
  libraries
- [`library_ion_weights()`](https://tomwenseleers.github.io/GCMSprime/reference/library_ion_weights.md)
  : Estimate ion-rarity weights from a spectral library
- [`identify_compounds()`](https://tomwenseleers.github.io/GCMSprime/reference/identify_compounds.md)
  : Annotate recovered GC/MS compounds and export library matches

## Shared tools and package overview

- [`reexports`](https://tomwenseleers.github.io/GCMSprime/reference/reexports.md)
  [`GIC`](https://tomwenseleers.github.io/GCMSprime/reference/reexports.md)
  [`EBIC`](https://tomwenseleers.github.io/GCMSprime/reference/reexports.md)
  [`MBIC`](https://tomwenseleers.github.io/GCMSprime/reference/reexports.md)
  : Objects exported from other packages
- [`GCMSprime`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSprime-package.md)
  [`GCMSprime-package`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSprime-package.md)
  : Scalable and Automated GC–MS Deconvolution, Quantification, and
  Spectral Library Search in R
