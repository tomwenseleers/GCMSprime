# Reproduce the examples

- [The two bundled CDF windows](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
  include the fitting calls, timing protocol and downloadable spectra.
- [The full-chromatogram comparison](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.html) uses 16 shared inputs, 10,000 scans, 550 channels and 1,000 true compounds. [Simulation, software settings and figures](../GCMSbenchmark/README.md).

- Run `Rscript benchmarks/library_search_weights.R path/to/MoNA-export-GC-MS_Spectra.msp output`
  to compare spectral scores. The
  [run record](results/library_search_weighted_run.csv) identifies the MSP
  snapshot and thread budget. The
  [identification guide](https://tomwenseleers.github.io/GCMSprime/articles/identification-theory.html#accuracy-and-speed)
  reports the held-out accuracy and median search times.
