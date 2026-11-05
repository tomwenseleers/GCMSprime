# Reproduce the examples

- [The two bundled CDF windows](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
  include the fitting calls, timing protocol and downloadable spectra.
- `Rscript benchmarks/16_samples_parallel.R` simulates the 16 chromatograms,
  writes them to CDF, fits each worker budget and checks support and spectra.
  The [multisample tutorial](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html)
  reports that simulation and the calibrated beer workflow.
- Run `Rscript benchmarks/library_search_weights.R path/to/MoNA-export-GC-MS_Spectra.msp output`
  to compare spectral scores. The
  [run record](results/library_search_weighted_run.csv) identifies the MSP
  snapshot and thread budget. The
  [identification guide](https://tomwenseleers.github.io/GCMSprime/articles/identification-theory.html#accuracy-and-speed)
  reports the held-out accuracy and median search times.

The November 2026 benchmark audit used GCMSprime
[`e10b8bf`](https://github.com/tomwenseleers/GCMSprime/commit/e10b8bf),
L0nmf [`1c85840`](https://github.com/tomwenseleers/L0nmf/commit/1c85840)
and L0glm [`d7bdece`](https://github.com/tomwenseleers/L0glm/commit/d7bdece).
The workstation and timing boundaries are specified in each tutorial.
