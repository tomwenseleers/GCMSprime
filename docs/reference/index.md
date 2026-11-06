# Package index

## Deconvolution, quantification and reports

- [`residuals(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`coef(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`fitted(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`nobs(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`logLik(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`AIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`BIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`GIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`EBIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`MBIC(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  [`print(`*`<GCMSdeco>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
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
- [`writeGCMSCDF()`](https://tomwenseleers.github.io/GCMSprime/reference/writeGCMSCDF.md)
  : Write a chromatogram matrix as an ANDI-MS CDF file

## Import and inspect spectral libraries

- [`readMSP()`](https://tomwenseleers.github.io/GCMSprime/reference/readMSP.md)
  [`openLibrary()`](https://tomwenseleers.github.io/GCMSprime/reference/readMSP.md)
  [`print(`*`<GCMSlibrary>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/readMSP.md)
  : Import a mass spectral library
- [`readGCMSID()`](https://tomwenseleers.github.io/GCMSprime/reference/readGCMSID.md)
  : Import a GCMS-ID predicted spectral library
- [`libraryMetadata()`](https://tomwenseleers.github.io/GCMSprime/reference/libraryMetadata.md)
  : Inspect spectral-library metadata
- [`spectra()`](https://tomwenseleers.github.io/GCMSprime/reference/spectra.md)
  : Extract spectral-library peaks or a sparse mass matrix
- [`retentionIndices()`](https://tomwenseleers.github.io/GCMSprime/reference/retentionIndices.md)
  : Extract retention indices and their provenance

## Calibrate retention indices and identify compounds

- [`calibrateRI()`](https://tomwenseleers.github.io/GCMSprime/reference/calibrateRI.md)
  [`predict(`*`<RIcalibration>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/calibrateRI.md)
  [`print(`*`<RIcalibration>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/calibrateRI.md)
  : Calibrate retention indices from an n-alkane standard
- [`annotateRI()`](https://tomwenseleers.github.io/GCMSprime/reference/annotateRI.md)
  : Add alkane-calibrated retention indices without a library search
- [`librarySearch()`](https://tomwenseleers.github.io/GCMSprime/reference/librarySearch.md)
  [`as.data.frame(`*`<GCMSmatches>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/librarySearch.md)
  [`print(`*`<GCMSmatches>`*`)`](https://tomwenseleers.github.io/GCMSprime/reference/librarySearch.md)
  : Identify deconvoluted spectra using retention indices and spectral
  libraries
- [`libraryIonWeights()`](https://tomwenseleers.github.io/GCMSprime/reference/libraryIonWeights.md)
  : Estimate ion-rarity weights from a spectral library
- [`identifyCompounds()`](https://tomwenseleers.github.io/GCMSprime/reference/identifyCompounds.md)
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
