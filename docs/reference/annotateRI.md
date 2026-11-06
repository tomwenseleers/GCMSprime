# Add alkane-calibrated retention indices without a library search

Annotate a single-sample or multisample deconvolution with retention
indices and optionally export annotated MSP and Excel files. Profiles
and spectra are retained unchanged.

## Usage

``` r
annotateRI(
  x,
  calibration = NULL,
  phase = c("DB5", "DB1", "DBWAX"),
  extrapolate = FALSE,
  msp_file = NULL,
  xlsx_file = NULL
)
```

## Arguments

- x:

  A
  [`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
  or
  [GCMSbatch](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
  result.

- calibration:

  A
  [`calibrateRI()`](https://tomwenseleers.github.io/GCMSprime/reference/calibrateRI.md)
  result. `NULL` uses each sample's `alkanes_file` and calibrated PSF.

- phase:

  Column class, default `"DB5"`.

- extrapolate:

  Linearly extend the RI curve beyond its alkane anchors. The default
  returns `NA` outside the anchor range. Extended values are flagged in
  the compound table and MSP comments.

- msp_file, xlsx_file:

  Optional output paths. For multisample results these name the
  consensus outputs; the corresponding per-sample files are also
  refreshed. Updated sample RDS objects remain in the batch output
  directory.

## Value

The annotated object, including the calibration, retention coordinates,
indices, column classes and extrapolation flags. Multisample consensus
RIs are medians over detected occurrences.

## Details

Peak coordinates use the 1%-height onset plus half the calibrated
narrow-peak width, measured at 1% height in the local window. A natural
cubic spline interpolates the corresponding alkane coordinates and
indices. MSP exports use NIST's phase-labelled `Retention_index` syntax:
`SemiStdNP` for DB5, `StdNP` for DB1 and `StdPolar` for DBWAX.

## References

NIST MS Search 2.3 manual, pp. 53–54, RI import/export:
<https://chemdata.nist.gov/dokuwiki/lib/exe/fetch.php?media=chemdata%3Anist17%3Anistms_ver23man.pdf>

## See also

[`calibrateRI()`](https://tomwenseleers.github.io/GCMSprime/reference/calibrateRI.md),
[`identifyCompounds()`](https://tomwenseleers.github.io/GCMSprime/reference/identifyCompounds.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cal <- calibrateRI("alkanes.cdf", psf=standard_psf)
fit <- annotateRI(fit, cal, phase="DB5",
  msp_file="sample_RI.msp", xlsx_file="sample_RI.xlsx")
fit$compounds[c("retention_time", "retention_index", "ri_extrapolated")]
} # }
```
