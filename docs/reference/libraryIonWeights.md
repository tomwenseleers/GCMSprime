# Estimate ion-rarity weights from a spectral library

Compute reusable smoothed inverse-document-frequency weights from
reference spectra. Rare ions receive larger weights than ubiquitous
fragments.

## Usage

``` r
libraryIonWeights(
  library,
  mass = 33:550,
  min_relative = 0.01,
  power = 1,
  ncores = NULL,
  batch_size = 4096L
)
```

## Arguments

- library:

  A GCMSlibrary.

- mass:

  Acquired nominal mass bins; default 33:550.

- min_relative:

  An ion is present when its intensity is positive and at least this
  fraction of the acquired-bin base peak, default 0.01.

- power:

  Exponent for inverse-frequency weighting, default 1. Zero gives
  uniform weights; larger values put more emphasis on rare ions.

- ncores:

  Native worker count; NULL uses physical cores.

- batch_size:

  Number of library spectra read per batch.

## Value

A named numeric vector of weights, normalized to mean one. Attributes
frequency, spectra, min_relative, power and seconds record its
derivation.

## Details

If L nonempty reference spectra are available and d_j contain ion j,
weights are proportional to \[1+log((L+1)/(d_j+1))\]^power. Frequencies
are counted across reference spectra, using the requested acquisition
bins. Zero-intensity spectra are omitted. Nominal duplicates are summed
before applying the relative threshold. Source intensities and totals
are preserved.

Build weights once and reuse them in librarySearch(ion_weights=...).
SQLite records are read in bounded batches; native frequency counts use
OpenMP. For validation, estimate weights from the reference catalogue
after excluding query records and near-duplicates. Uniform weights
remain the search default.

## See also

\[librarySearch()\]

## Examples

``` r
f <- tempfile(fileext=".msp")
writeLines(c("Name: A", "Num Peaks: 2", "43 100;57 20;",
             "Name: B", "Num Peaks: 2", "43 100;91 50;"), f)
lib <- readMSP(f)
w <- libraryIonWeights(lib, mass=c(43,57,91), ncores=1)
w
```
