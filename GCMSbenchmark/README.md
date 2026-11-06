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

| Method | Requirements and controls |
| --- | --- |
| `GCMSprime` | Install GCMSprime and logconNMF. Pass `gcmsDeco()` options such as `precluster`, `baseline_df` and `psf_control` through `...`. |
| `MS-DIAL` | Set `control[["MS-DIAL"]]` to a list with `exe` and optionally `powershell` (PowerShell 7 worked with the Windows console distribution), `method_file` or named `parameters`. The default input is mzML. For direct CDF use `format = "cdf"` and add the official NetCDF runtime through `netcdf` if needed. |
| `mzmine` | Set `control$mzmine$exe` to `mzmine_console.exe` and `user_file` to the existing machine-bound `.mzuser4` file when the CLI requires it. The included mzmine 4.11 batch uses the ADAP-inspired GC-EI workflow and exports its quantification CSV and spectral MGF locally. Supply a different `.mzbatch` with `batch`. |
| `eRah` | Install eRah and optionally future for multiple samples. Use `control$eRah$deconvolution` and `control$eRah$alignment` to set the arguments to `setDecPar()` and `setAlPar()`. |
| `AMDIS` | Set `control$AMDIS$exe` to `AMDIS_32.exe`. The adapter stages mzXML and runs serially. To rescore existing files, supply one `.ELU` in `control$AMDIS$elu` for each CDF in the same order. |
| `xcms+CAMERA` | Install both Bioconductor packages. xcms detects ion peaks with a matched filter and CAMERA groups them within each sample by coelution and EIC correlation. It returns summed ion integrals, pseudospectra, and representative-ion traces. Set `control[["xcms+CAMERA"]]` to override `detection`, `grouping`, or `min_ions` (default 3). |
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
