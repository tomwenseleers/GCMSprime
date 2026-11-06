# Inspect spectral-library metadata

Inspect spectral-library metadata

## Usage

``` r
libraryMetadata(x, which = NULL, raw = FALSE)
```

## Arguments

- x:

  A GCMSlibrary object.

- which:

  Spectrum indices; NULL selects every record.

- raw:

  Return a list of original tag/value tables, preserving repeated tags.

## Value

A metadata data frame, or a list of raw field tables.
