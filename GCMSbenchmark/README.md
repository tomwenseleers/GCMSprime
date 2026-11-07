# GCMSbenchmark

Compare GC/MS deconvolution workflows on the same ANDI-MS files. The package
collects timings, component calls, spectra, native peak areas and exported
elution profiles. Simulated truth adds reproducible recovery metrics.

## Installation

```r
remotes::install_github("tomwenseleers/GCMSprime", subdir = "GCMSbenchmark")
library(GCMSbenchmark)
b <- GCMSbenchmark("sample.cdf", "benchmark-run",
  method = c("GCMSprime", "eRah"), ncores = 8L)
b
```

The [installation guide](INSTALL.md) covers R packages, MS-DIAL, mzmine,
AMDIS, Python, MSHub and NetCDF. AMDIS execution requires Windows; ELU
exports can be imported on other platforms. Use a fresh output directory.
Each engine has a separate status and diagnostic message.

## Benchmark on 16 simulated chromatograms

Each sample has 1,200 scans and 258 mass channels. The 48 reference spectra
retain their library retention indices. Compounds are absent from about
half the samples, giving **380 true occurrences**. One compound has a broad
acid-like profile; the others are exponentially modified Gaussian peaks.
Higher concentrations broaden peaks towards later times from approximately
fixed starts. Successive samples drift slightly to the right. A smooth
background and quasi-Poisson noise with dispersion 2 complete the simulation.

GCMSprime recovered **372 of 380 occurrences with just 2 false positives**.
Its peak areas achieved adjusted R² **0.993** on the log₁₀ scale.

| Method | Seconds | TP | FP | F1 |
|:--|--:|--:|--:|--:|
| GCMSprime | 70.48 | 372 | 2 | 0.987 |
| MS-DIAL | 21.08 | 364 | 56 | 0.910 |
| mzmine | 41.83 | 348 | 91 | 0.850 |
| eRah | 10.73 | 347 | 258 | 0.705 |
| AMDIS | 36.04 | 351 | 492 | 0.574 |
| xcms + CAMERA | 26.97 | 349 | 34 | 0.915 |
| MSHub | 14.56 | 331 | 422 | 0.584 |

Times are wall-clock analysis times on an HP ZBook Fury G8 with an Intel
Core i9-11950H (8 physical cores, 16 logical threads), Windows and R 4.5.0.
Input conversion times are recorded separately. AMDIS runs serially and
xcms + CAMERA uses four workers. GCMSprime receives the generating dispersion
of 2; the known RI maps are used only for scoring.

[![Estimated versus true log10 peak areas with 95% prediction intervals](../GCMSprime/man/figures/simulated-16-area-predictions.png)](../GCMSprime/man/figures/simulated-16-area-predictions.png)

Points are RI-constrained reciprocal true matches with positive areas.
Each engine retains its native area units. The dashed 1:1 line applies
to integrated counts. Prediction intervals come from the log–log regression.
Model-ion areas depend on spectral response as well as integrated abundance.
The result tables also report area errors across samples after removing
each compound's native response factor.

[![Reciprocal spectral cosine hits on a logit-spaced axis](../GCMSprime/man/figures/simulated-16-spectral-violin.png)](../GCMSprime/man/figures/simulated-16-spectral-violin.png)

The library includes all 48 compounds, including those absent from a
sample. A pair must be a reciprocal best spectral hit within 30 RI units
and have cosine similarity of at least 0.70. Crosses identify hits to absent
compounds. Every other call is a false positive. Equal cosine scores are
ordered by RI difference and then ID.

AMDIS returned **843 components**, of which **351** were reciprocal
hits to true occurrences. **489 additional calls** shared their best eligible
reference hit with another call in the same sample. Its area regression
has adjusted R² **0.840**, conditional on the 351 selected true matches;
these are the individual returned areas, without adding duplicate calls.

[![Observed and true mean signal for one simulated sample](../GCMSprime/man/figures/simulated-16-sample-heatmap.png)](../GCMSprime/man/figures/simulated-16-sample-heatmap.png)

[![True and returned acid profiles across samples](../GCMSprime/man/figures/simulated-16-acid-profiles.png)](../GCMSprime/man/figures/simulated-16-acid-profiles.png)

The [benchmark article](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.html)
explains the scoring, profile exports and workflow differences.

### Reproduce the figures or rerun the comparison

The package includes the seeded simulation, engine settings, scoring tables
and plotting code. Recreate all four figures from the saved results with:

```r
bench <- system.file("benchmarks", package = "GCMSbenchmark", mustWork = TRUE)
source(file.path(bench, "plot-simulated-16.R"), encoding = "UTF-8")
plotBenchmark16(file.path(bench, "results"), "benchmark-figures")
```

To regenerate the CDFs and run the R workflows:

```r
source(file.path(bench, "simulated-16.R"), encoding = "UTF-8")
fit <- runBenchmark16("benchmark-16", method = c("GCMSprime", "xcms+CAMERA"),
  ncores = 4L)
plotBenchmark16("benchmark-16")
```

For all seven engines, supply their installed locations:

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

The run uses seed 9118 and saves truth, the reference library, native engine
exports, per-call identifications and per-pair area/profile metrics. The
reference includes compounds absent from each sample. Their presence is
checked after reciprocal matching. MSHub's aligned RTs use the generating
mean reference clock for RI conversion.

The tested versions are MS-DIAL 5.5.261105, mzmine 4.11.0, eRah 2.2.0,
AMDIS 2.71, xcms 4.6.4, CAMERA 1.64.0 and standalone MSHub revision
`59fe5d5d7ad5`. The MS-DIAL method file and mzmine batch are included.
AMDIS uses component width 16 scans and its local `ONSITE.INI` settings
`SENSIT=60`, `RESOLUTION=0`, `PEAKSHAPE=0`. The xcms + CAMERA width and
grouping settings were selected on odd samples of a regular-peak pilot
simulation and retained for this acid-containing comparison.

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

Profile cosine compares shapes on the original time grid. Area adjusted
R² comes from a log–log regression on positive true matches; `area_pairs`
gives its sample size. Duplicate calls are not summed into a matched area.
The summary reports their count alongside total false positives.
