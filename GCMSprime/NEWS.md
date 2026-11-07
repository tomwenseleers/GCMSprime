# GCMSprime 0.1.0.9000

## Chromatograms and calibration

* ANDI-MS CDF import and export preserve scan timing, acquired mass channels
  and sample metadata.
* Standard-assisted or blind constant/time-varying PSF estimation, pooled
  noise calibration and smooth background initialization.
* Reusable TIC-based standard calibration preserves physical scan times and
  builds each sample's sparse PSF without repeatedly fitting the standard.

* Contained broad-peak fragments use native local refinement and matched
  unmerged controls, with permissive initialization and correlation-aware
  rank selection on every original ion channel.

## Deconvolution and quantification

* Single-sample sparse log-concave NMF resolves overlapping compounds,
  flexible profiles and full mass spectra with automatic best rank tuning.
* Optional aggregation to a small number of ion groups accelerates
  initialization; original-channel likelihood refinement completes the fit.
* Directory-based multisample analysis preserves acquisition order and group
  labels, selects alkane standards, aligns samples and estimates shared spectra
  with sample-specific areas and profiles.

* `ncores` is a total logical-CPU budget for concurrent samples and native
  fitting. A shared reusable ion-group map and bounded coupled updates keep
  memory use controlled while retaining full spectra and sample-specific PSFs.
* Stage timings, thread-allocation metadata and a sixteen-sample simulated
  example accompany serial/parallel result checks.
* `simulateGCMSSamples()` provides RI-aware simulation with shared spectra,
  concentration-dependent broadening, optional acid peaks and sample drift.

## Spectral libraries and identification

* Streaming MSP and GCMS-ID CSV import with sparse in-memory or SQLite storage,
  chemical metadata, synonyms and retention-index provenance.
* Diagnostic-ion alkane assignment, natural-cubic onset-based retention-index
  calibration and column-specific constraints for DB1, DB5 and DBWAX classes.
* Retention-index annotation can be exported directly to phase-labelled NIST
  MSP and Excel files, with extrapolated values explicitly flagged.
* Native spectral search with Stein–Scott, cosine/angular, Hellinger,
  Jensen–Shannon and multinomial comparisons, plus optional ion-rarity weights.

## Output and diagnostics

* A common 0.1% base-peak ion cutoff is used for reported spectra, library
  queries and MSP/Excel output. Set ion_threshold=0 to retain all positive
  ions; fitted peak areas, profiles and information criteria are preserved.

* NIST-compatible MSP, Excel compound/area tables and consensus multisample
  libraries, with retained ion totals and ranked library matches.
* Observed/fitted/Pearson-residual heatmaps, chromatogram traces, profiles and
  spectra in temporal order, and compact objects that reload observations.
* Reproducible caching, experiment manifests and an integrated workflow from
  sample files to quantified and annotated compounds.
