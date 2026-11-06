# Import a GCMS-ID predicted spectral library

Import a GCMS-ID predicted spectral library

## Usage

``` r
readGCMSID(
  file,
  backend = c("auto", "memory", "sqlite"),
  db_file = NULL,
  batch_size = 1000L,
  max_spectra = Inf
)
```

## Arguments

- file:

  Extracted GCMS-ID CSV file from its public downloads page.

- backend, db_file, batch_size:

  Storage options as in \[readMSP()\].

- max_spectra:

  Maximum number of CSV records to import.

## Value

A GCMSlibrary object with predicted spectra and predicted RI provenance.

## Details

The GCMS-ID columns compound_name, inchikey, compound_formula,
retention_index_value, exact_mass and peak_list are read in batches. The
StdNP, SemiStdNP and StdPolar name suffixes identify DB1, DB5 and DBWAX
phase classes. Records with different RI predictions are retained
separately. Unknown suffixes retain an unknown phase. All source columns
are preserved.

## See also

\<https://gcms-id.ca/downloads\>
