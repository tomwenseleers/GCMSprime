# GCMSprime 0.1.0

## Chromatograms and calibration

* ANDI-MS CDF import and export preserve scan timing, acquired mass channels
  and sample metadata.
* Standard-assisted or blind constant/time-varying PSF estimation, pooled
  noise calibration and smooth background initialization.

## Deconvolution and quantification

* Single-sample sparse log-concave NMF resolves overlapping compounds,
  flexible profiles and full mass spectra with automatic best rank tuning.
* Optional aggregation to a small number of ion groups accelerates
  initialization; original-channel likelihood refinement completes the fit.
* Directory-based multisample analysis preserves acquisition order and group
  labels, selects alkane standards, aligns samples and estimates shared spectra
  with sample-specific areas and profiles.

## Spectral libraries and identification

* Streaming MSP and GCMS-ID CSV import with sparse in-memory or SQLite storage,
  chemical metadata, synonyms and retention-index provenance.
* Diagnostic-ion alkane assignment, smooth onset-based retention-index
  calibration and column-specific constraints for DB1, DB5 and DBWAX classes.
* Native spectral search with Stein–Scott, cosine/angular, Hellinger,
  Jensen–Shannon and multinomial comparisons, plus optional ion-rarity weights.

## Output and diagnostics

* NIST-compatible MSP, Excel compound/area tables and consensus multisample
  libraries, with retained ion totals and ranked library matches.
* Observed/fitted/Pearson-residual heatmaps, chromatogram traces, profiles and
  spectra in temporal order, and compact objects that reload observations.
* Reproducible caching, experiment manifests and an integrated workflow from
  sample files to quantified and annotated compounds.
