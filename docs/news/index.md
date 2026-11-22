# Changelog

## GCMSprime 0.1.0.9007

- Single-sample and batch Poisson and quasi-Poisson analyses use the
  conditional solver included in logconNMF.

## GCMSprime 0.1.0.9006

- Multisample deconvolution uses the original-ion overlap refit in
  logconNMF. Peak-area benchmarks report ordinary least-squares adjusted
  R² and pointwise 95% prediction intervals on log-transformed areas.

## GCMSprime 0.1.0.9005

- Compatible count-model batches read CDFs and estimate noise in
  parallel. After the grouped PSF regressions finish, the joint
  original-ion fit uses the same total core budget for sample profiles
  and shared spectra. `batch$workers` controls concurrent file imports
  and grouped regressions. `ncores=1L` runs the complete pipeline
  serially.

## GCMSprime 0.1.0.9004

- [`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  uses logconNMF’s conditional count solver for compatible analyses.
  Repeated injections can share a count-conserving PSF calibration,
  grouped initialization and a sparse original-ion spectrum per
  compound. The single-sample and 16-sample tutorials report timings and
  recovery from the fitted examples.
- `algorithm` selects the solver and `objective.trace` records its
  conditional objective.

## GCMSprime 0.1.0.9003

- Full-chromatogram component initialization retries Cholesky
  workspace-budget failures with diagonally preconditioned Eigen CG,
  reporting the switch and retaining the memory limit and selection
  settings. Requires logconNMF 0.1.0.9003.
- `gcmsDeco(control=list(initial_solver=...))` exposes L0glm solver
  controls, including `method`, `precond` and `max.memory` in GiB.
  `initial_fallback=FALSE` disables the automatic retry (GCMSprime issue
  [\#1](https://github.com/tomwenseleers/GCMSprime/issues1)).
- With `batch=list(cache=FALSE)`, coupled profile refinement reads local
  windows from disk instead of creating temporary whole-sample
  observation caches.

## GCMSprime 0.1.0.9002

### Chromatograms and calibration

- ANDI-MS CDF import and export preserve scan timing, acquired mass
  channels and sample metadata.

- Standard-assisted or blind constant/time-varying PSF estimation,
  pooled noise calibration and smooth background initialization.
  Automatic background knots follow the narrow-peak PSF scale and allow
  slowly changing spectra. Algorithm-versioned cache keys prevent reuse
  of obsolete initial fits.

- Reusable TIC-based standard calibration preserves physical scan times
  and builds each sample’s sparse PSF without repeatedly fitting the
  standard.

- Contained broad-peak fragments use native local refinement and matched
  unmerged controls, with permissive initialization and
  correlation-aware rank selection on every original ion channel.

### Deconvolution and quantification

- Single-sample sparse log-concave NMF resolves overlapping compounds,
  flexible profiles and full mass spectra with automatic best rank
  tuning.

- Optional aggregation to a small number of ion groups accelerates
  initialization and compatible Poisson/quasi-Poisson restart screens.
  Original-channel candidate/control refinement and likelihood checks
  complete the fit, retaining the full spectra and fixed nuisance model.

- Directory-based multisample analysis preserves acquisition order and
  group labels, selects alkane standards, aligns samples and estimates
  shared spectra with sample-specific areas and profiles.

- Coupled refinement reuses likelihoods and bounded
  count-data/normalizer caches, using the optimized logconNMF backend.
  Caches are released on return and can be disabled with
  `batch=list(cache=FALSE)`. The original-channel criterion, iteration
  budgets, grouping rules and shape constraints are retained.

- `ncores` is a total logical-CPU budget for concurrent samples and
  native fitting. Automatic directory budgets use logical cores; single
  chromatograms use physical cores. A shared reusable ion-group map and
  bounded coupled updates keep memory use controlled while retaining
  full spectra and sample-specific PSFs.

- Stage timings, thread-allocation metadata and a sixteen-sample
  simulated example accompany serial/parallel result checks.

- [`simulateGCMSSamples()`](https://tomwenseleers.github.io/GCMSprime/reference/simulateGCMSSamples.md)
  provides RI-aware simulation with shared spectra,
  concentration-dependent broadening, optional acid peaks and sample
  drift.

### Spectral libraries and identification

- Streaming MSP and GCMS-ID CSV import with sparse in-memory or SQLite
  storage, chemical metadata, synonyms and retention-index provenance.
- Diagnostic-ion alkane assignment, natural-cubic onset-based
  retention-index calibration and column-specific constraints for DB1,
  DB5 and DBWAX classes.
- Retention-index annotation can be exported directly to phase-labelled
  NIST MSP and Excel files, with extrapolated values explicitly flagged.
- Native spectral search with Stein–Scott, cosine/angular, Hellinger,
  Jensen–Shannon and multinomial comparisons, plus optional ion-rarity
  weights.

### Output and diagnostics

- A common 0.1% base-peak ion cutoff is used for reported spectra,
  library queries and MSP/Excel output. Set ion_threshold=0 to retain
  all positive ions; fitted peak areas, profiles and information
  criteria are preserved.

- NIST-compatible MSP, Excel compound/area tables and consensus
  multisample libraries, with retained ion totals and ranked library
  matches.

- Observed/fitted/Pearson-residual heatmaps, chromatogram traces,
  profiles and spectra in temporal order, and compact objects that
  reload observations.

- Reproducible caching, experiment manifests and an integrated workflow
  from sample files to quantified and annotated compounds.

- Nested stage timing avoids forced garbage collection between local
  fits; ordinary R memory management remains active.
