# GCMSbenchmark

Run several GC/MS deconvolution workflows on the same ANDI-MS files. The
benchmark collects component calls, spectra, elution profiles and native
quantification where the backend exports them. Simulated truth is optional.

Install the package from the directory that contains this README.

```r
install.packages(".", repos = NULL, type = "source")
library(GCMSbenchmark)

b <- GCMSbenchmark(
  cdf_file = "input/sample.cdf",
  output_dir = "benchmark-run",
  method = c("GCMSprime", "eRah"),
  ncores = 8L
)
b
b$methods$GCMSprime$peaks
```

Use a fresh output directory for each run. The benchmark keeps native exports
alongside the standardized data, summary CSV and any per-peak truth matches.
It records a separate status and message for each backend.

## Benchmark on 16 simulated chromatograms

Each chromatogram has 1,200 scans and 258 mass channels. The 48 simulated
compounds occur in different subsets of the samples, giving 380 true peak
occurrences under retention-time drift, concentration-dependent peak
shapes and overdispersed count noise.
**GCMSprime recovered all 380 occurrences with zero false calls** and
an adjusted R² of **0.99997** for integrated peak areas on the
log₁₀ scale.

[![Estimated versus true log10 peak areas, with 95% prediction intervals](../GCMSprime/man/figures/simulated-16-area-predictions.png)](../GCMSprime/man/figures/simulated-16-area-predictions.png)

Each point is a positive, one-to-one matched area. The estimated areas
use the units reported by each engine. The dashed 1:1 line marks direct
agreement for integrated counts.

[![Individual best spectral cosine hits and their distributions on a logit axis](../GCMSprime/man/figures/simulated-16-spectral-violin.png)](../GCMSprime/man/figures/simulated-16-spectral-violin.png)

Each point is a true occurrence’s best spectral match in the sample.
The right panel requires a candidate within ±7.2 seconds. The axis is
logit-spaced and labelled in cosine units.

The [GCMSprime software comparison](https://tomwenseleers.github.io/GCMSprime/articles/software-comparison.html)
describes the workflows and the criteria used to score them.

### Reproduce the figures or rerun the benchmark

Install from this directory or run
`remotes::install_github("tomwenseleers/GCMSprime", subdir = "GCMSbenchmark")`.
The package includes the simulation recipe, engine settings, scoring tables
and plotting code. To recreate the figures from the saved tables, install
`ggplot2` and run the following.

```r
bench <- system.file("benchmarks", package = "GCMSbenchmark", mustWork = TRUE)
source(file.path(bench, "plot-simulated-16.R"), encoding = "UTF-8")
plotBenchmark16(file.path(bench, "results"), "benchmark-figures")
```

To regenerate the CDFs and fit the simulation, install GCMSprime, logconNMF
and the selected engines. The following runs the two R workflows.

```r
source(file.path(bench, "simulated-16.R"), encoding = "UTF-8")
fit <- runBenchmark16("benchmark-16", method = c("GCMSprime", "xcms+CAMERA"),
  ncores = 4L)
fit
plotBenchmark16("benchmark-16")
```

For all seven workflows, use `method = "all"` and point the runner to your
local installations. The paths below are examples. To run a subset, pass
its method names and any controls it needs.

```r
control <- list(
  "MS-DIAL" = list(exe = "C:/tools/MS-DIAL/MSDIALCUI.exe", powershell = "pwsh"),
  mzmine = list(exe = "C:/tools/mzmine/mzmine_console.exe",
    user_file = "C:/path/to/your.mzuser4"),
  AMDIS = list(exe = "C:/tools/AMDIS32/AMDIS_32.exe"),
  MSHub = list(exe = "C:/tools/python/python.exe", source = "C:/tools/mshub_process")
)
fit <- runBenchmark16("benchmark-16-all", method = "all", control = control,
  ncores = 16L)
plotBenchmark16("benchmark-16-all")
```

The benchmark used Windows and R 4.5.0, MS-DIAL 5.5.261105, mzmine 4.11.0,
eRah 2.2.0, AMDIS 2.71, xcms 4.6.4, CAMERA 1.64.0, and standalone MSHub
revision `59fe5d5d7ad5`. The MS-DIAL method file and mzmine batch are
included. AMDIS uses its local deconvolution settings. The reported run
used component width 16 scans and the `ONSITE.INI` settings
`SENSIT=60`, `RESOLUTION=0` and `PEAKSHAPE=0`.
The xcms + CAMERA filter and grouping settings were selected on the odd
samples and evaluated separately on the even samples.

The runner regenerates the data with seed 9118, saves the true profiles and
spectra, and keeps each engine’s native exports. The shared scoring rule
requires a retention-time difference of at most 7.2 seconds and spectral
cosine of at least 0.70. For MSHub, the true apices are evaluated on its mean
alignment clock. The best-hit figure allows a recovered spectrum to be the
best match for more than one true occurrence.

## Engine setup

| Method | Requirements and controls |
| --- | --- |
| `GCMSprime` | Install GCMSprime and logconNMF. Pass `gcmsDeco()` options such as `precluster`, `baseline_df` and `psf_control` through `...`. |
| `MS-DIAL` | Set `control[["MS-DIAL"]]` to a list with `exe` and optionally `powershell` (PowerShell 7 worked with the Windows console distribution), `method_file` or named `parameters`. The default input is mzML. For direct CDF use `format = "cdf"` and add the official NetCDF runtime through `netcdf` if needed. |
| `mzmine` | Set `control$mzmine$exe` to `mzmine_console.exe` and `user_file` to the existing machine-bound `.mzuser4` file when the CLI requires it. The included mzmine 4.11 batch uses the ADAP-inspired GC-EI workflow and exports its quantification CSV and spectral MGF locally. Supply a different `.mzbatch` with `batch`. |
| `eRah` | Install eRah and optionally future for multiple samples. Use `control$eRah$deconvolution` and `control$eRah$alignment` to set the arguments to `setDecPar()` and `setAlPar()`. |
| `AMDIS` | Set `control$AMDIS$exe` to `AMDIS_32.exe`. The adapter stages mzXML and runs serially. To rescore existing files, supply one `.ELU` in `control$AMDIS$elu` for each CDF in the same order. |
| `xcms+CAMERA` | Install both Bioconductor packages. xcms detects ion peaks with a matched filter and CAMERA groups them within each sample by coelution and EIC correlation. It returns summed ion integrals, pseudospectra, and representative-ion traces. Use `control[["xcms+CAMERA"]]` to set `detection` (for example `fwhm` in seconds or `snthresh`), `fwhm_grouping` (arguments to `CAMERA::groupFWHM()`), `grouping` (arguments to `CAMERA::groupCorr()`, or `FALSE` for FWHM grouping alone), or `min_ions` (default 3). |
| `MSHub` | For a batch of at least 10 samples, set `control$MSHub$exe` to Python and `source` to the standalone [MSHub repository](https://bitbucket.org/iAnalytica/mshub_process/src/master/). The runtime needs NumPy, SciPy, h5py, scikit-learn and netCDF4. An optional `pythonpath` selects an isolated dependency directory. |

Each truth record supplies an ID, an apex in minutes, an area in integrated
counts and an acid indicator for each compound. Supply `H` with scans as
rows, `P` with compounds as rows, and the original `time` and `mass` grids.
Use a named list of these records, keyed by CDF basename without its
extension. Optionally include `reference_apex` for methods whose outputs
use aligned retention times.

The scorer maximises the number of one-to-one matches within the retention
tolerance and above the spectral cosine cutoff, then maximises total
spectral cosine. Profile cosine compares shapes on the original time
grid. Areas in integrated counts can be compared directly. Other backend
signals are compared after centring the log area ratios and, with repeated
samples, after centring within each compound.

The adjusted R² for peak areas uses a log–log regression on positive matched
areas. The summary reports how many matches contribute to this regression.
