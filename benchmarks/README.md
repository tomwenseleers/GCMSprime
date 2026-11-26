# Reproduce the examples

- [The two bundled CDF windows](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
  include the fitting calls, timing protocol and downloadable spectra.
- `Rscript benchmarks/16_samples_parallel.R` checks worker consistency
  on 16 samples with 48 reference compounds and 1,200 scans. The conditional
  solver recovered **true rank 48 exactly (estimated rank 48)** and all 380 occurrences with zero false positives. The
  shared spectra, areas and GIC agreed across 1, 8 and 16 requested threads.
  Wall times were 65.67, 31.59 and 30.47 seconds, respectively.
  [Worker-count results](results/16_samples_parallel.csv)
- `benchmark_fullsize_16.R` reproduces the 16-sample resource validation with
  10,000 scans, 550 recorded mass channels and about 500 compounds per sample.
  The fast local search recovered the true rank of 1,000 exactly.
  From the repository root, run:

  ```r
  source("benchmarks/benchmark_fullsize_16.R")
  benchmarkFullsize16(stage="prepare")
  benchmarkFullsize16(stage="fit", ncores=16L, workers=6L, coupling=list(penalty_search=TRUE))
  benchmarkFullsize16(stage="score")
  ```

  Large CDF inputs and fitted stores go in the ignored `codex_sandbox/fullsize16`
  directory by default; `directory=` selects another persistent location.
  The [multisample tutorial](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html)
  reports current resource use and quality. Compact published results are in
  [`results/fullsize16`](results/fullsize16).
- In the separate RI-aware seven-workflow comparison, the optional support
  search recovered **true rank 100 exactly (estimated rank 100)**. The simulation uses 100 reference compounds,
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

Recreate the published CPU and memory figure directly from the recorded stage timings and process-tree samples.

```r
source("benchmarks/plot_fullsize_resources.R")
plotFullsize16Resources(file="fullsize-16-resources.png", physical_cores=8L, logical_cores=16L)
```

Recreate the occurrence counts and matched-area plot (violet diamonds mark acid peaks) from the same published tables.

```r
source("benchmarks/plot_fullsize_recovery.R")
plotFullsize16Recovery()
```
