# Installing the benchmark engines

GCMSbenchmark runs the engines installed on your computer. Choose the methods
you want to compare, then use the steps below to install their requirements.
Application installers and licence registration are completed in the relevant
application.

## R packages

Use a current R version. Compiling L0glm, logconNMF and GCMSprime on Windows
requires the [Rtools version matching R](https://cran.r-project.org/bin/windows/Rtools/).
For R 4.5, use Rtools45. macOS requires the command-line developer tools and
a compatible Fortran compiler. Linux requires a C++17 toolchain and the
development libraries requested by the R package installers.

```r
install.packages(c("remotes", "BiocManager", "clue", "processx", "ggplot2",
                   "erah", "future", "ncdf4"))
remotes::install_github("tomwenseleers/L0glm", subdir = "L0glm")
remotes::install_github("tomwenseleers/logconNMF", subdir = "logconNMF")
remotes::install_github("tomwenseleers/GCMSprime", subdir = "GCMSprime")
remotes::install_github("tomwenseleers/GCMSprime", subdir = "GCMSbenchmark")
BiocManager::install(c("mzR", "xcms", "CAMERA", "BiocParallel"))
```

The R workflows use GCMSprime, eRah, and xcms + CAMERA. mzR supplies the
format conversion used for several external engines.
[Bioconductor installation](https://bioconductor.org/install/) ·
[eRah manual](https://cran.r-project.org/web/packages/erah/vignettes/erah-manual.html)

## MS-DIAL

Download the GC–MS-capable console distribution from the developers' site.
The benchmark was verified with the Windows .NET Framework console build
5.5.261105. Supply the executable and a PowerShell executable able to load
its export assemblies. PowerShell 7 was used for the benchmark.

```r
control[["MS-DIAL"]] <- list(
  exe = "C:/tools/MS-DIAL/MSDIALCUI.exe",
  powershell = "pwsh")
```

The adapter converts the CDFs to mzML. To pass CDFs directly, set
`format="cdf"` and, if required by your console build, install the official
NetCDF runtime of matching bitness. Set `netcdf` to the directory containing
its DLLs. This is used for the child process.

```r
control[["MS-DIAL"]]$format <- "cdf"
control[["MS-DIAL"]]$netcdf <- "C:/tools/netCDF/bin"
```

[MS-DIAL downloads](https://systemsomicslab.github.io/compms/msdial/main.html) ·
[Console instructions](https://systemsomicslab.github.io/compms/msdial/consoleapp.html) ·
[Unidata NetCDF-C downloads](https://downloads.unidata.ucar.edu/netcdf-c/)

## mzmine

Install the mzmine distribution for your operating system. It includes its
Java runtime. Open the application once and complete registration and sign-in
under Users. The console needs an available licence for the account.
The published batch was verified with mzmine 4.11.

```r
control$mzmine <- list(
  exe = "C:/tools/mzmine/mzmine_console.exe",
  user_file = "C:/path/to/your.mzuser4")
```

Use the local account file created by mzmine. Keep that file private.
Closing idle mzmine instances releases their licence claims. If the JVM
reserves too much memory for a small benchmark, a session-specific
`JAVA_TOOL_OPTIONS="-Xms256m -Xmx4g"` caps its heap.

[Installation](https://mzmine.github.io/mzmine_documentation/latest/getting_started.html) ·
[Licence plans](https://mzio.io/plans/)

## AMDIS

AMDIS is a Windows application. Install the standalone distribution from
NIST, then provide its executable. The adapter processes the files serially
and reads the resulting ELU exports.

```r
control$AMDIS <- list(exe = "C:/tools/AMDIS32/AMDIS_32.exe")
```

Set the documented deconvolution options in AMDIS before the benchmark.
The recipe records the settings used for the published run.
On other operating systems, `method="all"` omits AMDIS execution. ELU
files generated on Windows can be imported with
`control$AMDIS$elu`, one path per input CDF.

[NIST distribution](https://chemdata.nist.gov/dokuwiki/doku.php?id=chemdata%3Aamdis) ·
[AMDIS manual](https://www.nist.gov/system/files/documents/srd/AMDISMan.pdf)

## Python and MSHub

Install Python and create an isolated environment. Download the standalone
MSHub source and select revision `59fe5d5d7ad5` to match this benchmark.
Install NumPy, SciPy, h5py, scikit-learn, netCDF4 and cftime in that environment.
The source also uses matplotlib, pandas and statsmodels.

```text
python -m venv mshub-env
mshub-env/Scripts/python.exe -m pip install numpy scipy h5py scikit-learn netCDF4 cftime matplotlib pandas statsmodels
```

On Linux or macOS, use `mshub-env/bin/python`. Point the adapter to that
interpreter and the directory containing MSHub's `proc/` tree.

```r
control$MSHub <- list(
  exe = "C:/tools/mshub-env/Scripts/python.exe",
  source = "C:/tools/mshub_process")
```

An existing Python environment can use an isolated dependency directory
through `control$MSHub$pythonpath`. MSHub requires a batch of at least ten
samples. The adapter includes compatibility handling for NumPy aliases
used by the standalone source.

The published run used Python 3.13.9, NumPy 2.3.5, SciPy 1.16.3,
h5py 3.15.1, scikit-learn 1.7.2 and netCDF4 1.7.4.

[Standalone source](https://bitbucket.org/iAnalytica/mshub_process/src/master/) ·
[GNPS workflow guide](https://ccms-ucsd.github.io/GNPSDocumentation/gc-ms-deconvolution/)

## Check the selected setup

```r
packages <- c("L0glm", "logconNMF", "GCMSprime", "GCMSbenchmark", "erah",
              "xcms", "CAMERA", "mzR", "ncdf4", "ggplot2")
setNames(vapply(packages, requireNamespace, logical(1), quietly = TRUE), packages)
file.exists(control[["MS-DIAL"]]$exe)
file.exists(control$mzmine$exe)
file.exists(control$AMDIS$exe)
file.exists(file.path(control$MSHub$source, "proc", "io", "importmsdata.py"))
```

The benchmark records each engine's status and error message in its summary.
Use a fresh output directory after correcting an installation or control.
