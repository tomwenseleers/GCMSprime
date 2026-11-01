# Extract spectral-library peaks or a sparse mass matrix

Extract spectral-library peaks or a sparse mass matrix

## Usage

``` r
spectra(
  x,
  which = NULL,
  format = c("matrix", "list"),
  nominal = FALSE,
  normalize = c("none", "l1", "l2")
)
```

## Arguments

- x:

  A GCMSlibrary object.

- which:

  Spectrum indices; select batches for large databases.

- format:

  Sparse matrix (mass rows, spectrum columns) or raw peak lists.

- nominal:

  Bin masses to the nearest integer; FALSE preserves exact masses.

- normalize:

  Optional unit-sum (l1) or unit-Euclidean-norm (l2) spectra.

## Value

A dgCMatrix or a list of two-column mass/intensity matrices.
