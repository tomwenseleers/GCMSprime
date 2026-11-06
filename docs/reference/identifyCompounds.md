# Annotate recovered GC/MS compounds and export library matches

Add retention indices and ranked library matches to an existing
deconvolution, without refitting its profiles or spectra.

## Usage

``` r
identifyCompounds(x, library, ..., msp_file = NULL, xlsx_file = NULL)
```

## Arguments

- x:

  A gcmsDeco object.

- library:

  A GCMSlibrary.

- ...:

  Search options passed to librarySearch(), including calibration,
  phase, ri_tolerance, method and n_matches.

- msp_file:

  Optional new annotated MSP output path.

- xlsx_file:

  Optional Excel output path with compound, match and calibration
  sheets.

## Value

The gcmsDeco object with identification (a GCMSmatches object),
retention indices and best-match columns in compounds. The first
requested method supplies the best-match columns; all requested methods
remain in identification\$matches.

## Details

The original compound names, factorization and fitted ion totals are
preserved. Library names are tentative annotations. MSP files include
the RI, column class, original fitted ion total and best-match
description. Excel includes every returned match and the accepted alkane
calibration anchors. Explicit assignments of unlabelled reference RI
phases are recorded in the compound table, match table, Run sheet and
MSP annotation. Supply output paths to export companion files from a
saved fit.

## See also

\[calibrateRI()\], \[librarySearch()\], \[gcmsDeco()\]

## Examples

``` r
# Given a saved deconvolution and an imported spectral library:
if (FALSE) { # \dontrun{
cal <- calibrateRI("alkanes.cdf", carbon=7:23)
identified <- identifyCompounds(fit, lib, calibration=cal, phase="DB5",
  ri_tolerance=30, n_matches=5, method=c("stein_scott","hellinger"),
  msp_file="sample_identified.msp", xlsx_file="sample_identified.xlsx")
identified$identification$matches
} # }
```
