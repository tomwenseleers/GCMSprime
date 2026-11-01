# Write a chromatogram matrix as an ANDI-MS CDF file

Export scan-by-mass intensities for
[`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
or external GC/MS software. The NetCDF classic file includes scan
indexing, total ion counts, mass ranges and ANDI-MS administrative and
instrument records.

## Usage

``` r
write_GCMS_cdf(
  y,
  cdf_file,
  time = as.numeric(rownames(y)),
  mass = as.numeric(colnames(y)),
  sample_name = tools::file_path_sans_ext(basename(cdf_file))
)
```

## Arguments

- y:

  Finite nonnegative numeric matrix, with scans in rows and nominal mass
  channels in columns. Zero channels are retained.

- cdf_file:

  Output `.cdf` path.

- time:

  Increasing scan times in minutes. Defaults to numeric row names.

- mass:

  Increasing nonnegative integer m/z values. Defaults to numeric column
  names.

- sample_name:

  Sample identity recorded in the header; defaults to the output
  filename without its extension.

## Value

Invisibly, the normalized output path.

## Details

Intensities are stored without rescaling or rounding. Integer data
fitting in a signed 32-bit integer use ANDI `Long` storage; other data
use `Double`. Scan times are converted to seconds in the file. Every
channel is written in every scan, including exact zeros. The file's
timestamp records export time; acquisition metadata are left
unspecified.

## See also

[`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md),
[`L0nmf::splogconNMF()`](https://rdrr.io/pkg/L0nmf/man/splogconNMF.html)

## Examples

``` r
# \donttest{
if (requireNamespace("ncdf4", quietly = TRUE)) {
  data(GCMS_data1)
  cdf <- tempfile(fileext = ".cdf")
  write_GCMS_cdf(GCMS_data1, cdf)
  unlink(cdf)
}
# }
```
