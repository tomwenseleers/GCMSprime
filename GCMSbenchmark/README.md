# GCMSbenchmark

Compare GC/MS deconvolution workflows on common ANDI-MS input files.
The package collects timings, component calls, spectra, total-ion area
equivalents and exported elution profiles. Native exports are retained.

## Installation

~~~r
remotes::install_github("tomwenseleers/GCMSprime", subdir="GCMSbenchmark")
library(GCMSbenchmark)
~~~

The [installation guide](INSTALL.md) covers all R packages, external engines
and platform requirements. AMDIS execution requires Windows.

## Benchmark on 16 full chromatograms

The same 16 simulated chromatograms were supplied to every program. Each contains **10,000 scans × 550 nominal mass channels** (m/z 30–579), covering approximately 50 minutes. Together they contain **88 million simulated ion intensities**. Each of the **1,000 distinct compounds** occurs in eight randomly chosen samples. The design therefore contains **8,000 positive sample–compound pairs**. Rank here counts the distinct generating compounds across the batch, each with a time-local elution profile in each positive sample. Of the compounds, 990 have exponentially modified Gaussian (EMG) peaks and ten have a rising shark-fin profile with an abrupt trailing edge. Higher concentrations broaden both types. The acquisition clock allows a smooth retention shift of up to 2.93 seconds between the first and last injections and approximately 2% variation in scan spacing.

GCMSprime, using simultaneous logconNMF, **recovered the true rank of 1,000 exactly** in **35.95 seconds per sample**. It recovered **8,000 of 8,000** compound occurrences and made **no false-positive calls**. The recovered peaks included **80 of the 80 acid occurrences**. MSHub returned **1,414 deconvolved components**, 41.4% above the true rank. The other five workflows produce peak lists or aligned feature groups. Their outputs do not define a common factorisation rank.

| Workflow | Seconds/sample | Estimated rank (true 1,000) | TP | FP | FN | Acid TP / 80 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GCMSprime | 35.95 | **1,000 (exact)** | 8,000 | 0 | 0 | 80 / 80 |
| MS-DIAL | 16.02 | — | 7,728 | 11,685 | 272 | 73 / 80 |
| mzmine GC-EI | 6.92 | — | 4,202 | 185 | 3,798 | 29 / 80 |
| eRah | 25.94 | — | 7,172 | 5,615 | 828 | 70 / 80 |
| AMDIS | 28.99 | — | 7,878 | 9,122 | 122 | 80 / 80 |
| xcms + CAMERA | 4.22 | — | 7,261 | 890 | 739 | 31 / 80 |
| MSHub | 15.05 | 1,414 | 7,179 | 12,184 | 821 | 50 / 80 |

Each recovered spectrum competes against all 1,000 references, including compounds absent from that sample. A true positive requires a reciprocal best match with spectral cosine ≥0.70, within ±30 simulated RI units, to a compound present in that sample. The score allows one credit per sample and reference. Extra calls, including duplicates, count as false positives. Sensitivity uses all 8,000 occurrences as its denominator. Acid recovery uses all 80.

Each program ran once, in sequence, on the same 64 GB Windows laptop with eight physical cores and 16 logical processors. Seconds per sample means wall time divided by 16 for the whole batch. The comparison includes any required input conversion or staging and the program’s native analysis. GCMSprime’s time includes reading and fitting an independent PSF standard, the coupled fit and native results export. Simulation, common-format result extraction and accuracy scoring are outside the timed comparison. One run per program gives a descriptive timing for these files and this machine.

Install the programs and set their local executable paths in the named list `control`. The included R recipe generates the common CDFs and their scoring key, fits the programs sequentially, then scores the saved fits. [Engine setup and control list](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/INSTALL.md).

```r
bench <- system.file("benchmarks", package = "GCMSbenchmark", mustWork = TRUE)
source(file.path(bench, "simulated-16.R"))
directory <- "full-chromatograms"
prepareFullBenchmark(directory)
methods <- c("GCMSprime", "MS-DIAL", "mzmine", "MSHub", "AMDIS", "xcms+CAMERA", "eRah")
for (m in methods) fitFullBenchmark(directory, m, control[[m]], ncores = 16L)
for (m in methods) scoreFullBenchmark(directory, m)
collectFullBenchmark(directory, methods, file.path(directory, "results"))
source(file.path(bench, "plot-simulated-16.R"))
plotBenchmark16(file.path(directory, "results"), file.path(directory, "figures"))
```

[CDF checksums](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/inputs.csv) · [Reference spectra](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/reference-library.rds) · [Occurrence and area matrices](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/settings.rds) · [Summary and timing components](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/benchmark.csv) · [Matched peaks](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/matches.csv) · [Per-acid recovery](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/acid-by-compound.csv) · [Recovery by abundance quartile](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/sensitivity-by-area.csv) · [Software versions](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/versions.csv) · [Serial run log](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/run-order.csv). The numerical input hashes cover scan times, masses and intensities. The whole-file hashes also include the CDF export timestamp.

[All figures](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.html).

## Engine setup

| Method | Requirements and controls |
| --- | --- |
| `GCMSprime` | Install GCMSprime and logconNMF. Pass `gcmsDeco()` options such as `precluster`, `baseline_df` and `psf_control` through `...`. |
| `MS-DIAL` | Set `control[["MS-DIAL"]]` to a list with `exe` and optionally `powershell` (PowerShell 7 worked with the Windows console distribution), `method_file` or named `parameters`. The default input is mzML. For direct CDF use `format = "cdf"` and add the official NetCDF runtime through `netcdf` if needed. |
| `mzmine` | Set `control$mzmine$exe` to `mzmine_console.exe` and `user_file` to the existing machine-bound `.mzuser4` file when the CLI requires it. The included mzmine 4.11 batch uses the ADAP-inspired GC-EI workflow and exports its quantification CSV and spectral MGF locally. Supply a different `.mzbatch` with `batch`. |
| `eRah` | Install eRah and optionally future for multiple samples. Use `control$eRah$deconvolution` and `control$eRah$alignment` to set the arguments to `setDecPar()` and `setAlPar()`. |
| `AMDIS` | On Windows, set `control$AMDIS$exe` to `AMDIS_32.exe`. The adapter stages mzXML and runs serially. To rescore existing files, supply one `.ELU` in `control$AMDIS$elu` for each CDF in the same order. |
| `xcms+CAMERA` | Install both Bioconductor packages. xcms detects ion peaks with a matched filter and CAMERA groups them within each sample by coelution and EIC correlation. It returns summed ion integrals, pseudospectra, and representative-ion traces. Use `control[["xcms+CAMERA"]]` to set `detection` (for example `fwhm` in seconds or `snthresh`), `fwhm_grouping` (arguments to `CAMERA::groupFWHM()`), `grouping` (arguments to `CAMERA::groupCorr()`, or `FALSE` for FWHM grouping alone), or `min_ions` (default 3). |
| `MSHub` | For a batch of at least 10 samples, set `control$MSHub$exe` to Python and `source` to the standalone [MSHub repository](https://bitbucket.org/iAnalytica/mshub_process/src/master/). The runtime needs NumPy, SciPy, h5py, scikit-learn and netCDF4. An optional `pythonpath` selects an isolated dependency directory. |

## Supply simulation truth

Each truth record contains `id`, `apex` in minutes, `area` in integrated
counts and a logical `acid` indicator. Supply `H` with scans as rows,
`P` with compounds as rows, and the original `time` and `mass` grids.
Key the records by input-file basename without its extension.

The default `matching="assignment"` maximizes the number of one-to-one
matches within `rt_tolerance` and above `min_cosine`, then total cosine.
For RI-constrained reciprocal scoring, use `matching="reciprocal"` and
`ri_tolerance=30`. Supply a full `reference` library with `id`, `P`, `mass`
and `ri`, plus an increasing `time`/`ri` calibration data frame as
`ri_calibration` in each truth record. Aligned outputs also need
`reference_ri_calibration`. The benchmark retains every call's best
eligible library hit and reciprocal flag in `identifications`.

Profile cosine compares shapes on the original time grid. Ordinary adjusted R² comes from a least-squares regression on the positive matched log₁₀ areas. The field `area_pairs` gives the number of pairs. Duplicate calls are not summed into a matched area.
The summary reports their count alongside total false positives.
