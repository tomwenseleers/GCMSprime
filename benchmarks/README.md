# Reproduce the examples

- [The two bundled CDF windows](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
  include the fitting calls, timing protocol and downloadable spectra.
- `Rscript benchmarks/16_samples_parallel.R` checks worker consistency
  on 16 samples with 48 reference compounds and 1,200 scans. The conditional
  solver recovered all 380 occurrences with zero false positives. The
  shared spectra, areas and GIC agreed across 1, 8 and 16 requested threads.
  Wall times were 71.70, 37.59 and 37.49 seconds, respectively.
  [Worker-count results](results/16_samples_parallel.csv)
- `benchmark_fullsize_16.R` reproduces the 16-sample resource validation with
  10,000 scans, 550 recorded mass channels and about 500 compounds per sample.
  From the repository root, run:

  ```r
  source("benchmarks/benchmark_fullsize_16.R")
  benchmarkFullsize16(stage="prepare")
  benchmarkFullsize16(stage="fit")
  benchmarkFullsize16(stage="score")
  ```

  Large CDF inputs and fitted stores go in the ignored `codex_sandbox/fullsize16`
  directory by default; `directory=` selects another persistent location.
  The [multisample tutorial](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html)
  reports current resource use and quality. Compact published results are in
  [`results/fullsize16`](results/fullsize16).
- The separate RI-aware seven-workflow comparison uses 100 reference compounds,
  1,200 scans, 258 channels and one acid-like compound. Its complete recipe and
  results are in [`GCMSbenchmark`](../GCMSbenchmark).
- Run `Rscript benchmarks/library_search_weights.R path/to/MoNA-export-GC-MS_Spectra.msp output`
  to compare spectral scores. The
  [run record](results/library_search_weighted_run.csv) identifies the MSP
  snapshot and thread budget. The
  [identification guide](https://tomwenseleers.github.io/GCMSprime/articles/identification-theory.html#accuracy-and-speed)
  reports the held-out accuracy and median search times.

The workstation, worker budgets, matching rules and timing boundaries are
specified alongside each result. The full-size rerun saves `sessionInfo.txt`.

On Windows, the optional standard-library Python monitor samples the whole
R worker process tree once per second. After preparing the inputs:

```text
python benchmarks/monitor_process.py codex_sandbox/fullsize16/resources -- C:/R/R-4.5.0/bin/x64/Rscript.exe --vanilla benchmarks/benchmark_fullsize_16.R fit
```

Use the Rscript path for the installed R version. The CSV records sampled
working/private memory, available physical RAM and cumulative CPU time.
