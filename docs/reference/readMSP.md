# Import a mass spectral library

Stream NIST-compatible MSP files into sparse in-memory or disk-backed
storage. OpenChrom Agilent and Shimadzu exports and MoNA/Golm MSP fields
are accepted.

## Usage

``` r
readMSP(
  file,
  backend = c("auto", "memory", "sqlite"),
  db_file = NULL,
  batch_size = 1000L,
  source = basename(file),
  encoding = "UTF-8",
  origin = c("unknown", "experimental", "predicted"),
  strict = TRUE,
  max_spectra = Inf
)

openLibrary(path)

# S3 method for class 'GCMSlibrary'
print(x, ...)
```

## Arguments

- file:

  Local MSP file. Compressed files should first be extracted.

- backend:

  Storage: auto uses SQLite for files larger than 100 MB.

- db_file:

  SQLite output path; required for the sqlite backend. Auto storage uses
  a new file beside the MSP when this is omitted.

- batch_size:

  Number of records buffered per transaction.

- source:

  Library provenance, defaulting to the input filename.

- encoding:

  Text encoding; UTF-8 input with Windows-1252 fallback is supported.

- origin:

  Whether spectra are experimental, predicted or of unknown origin.

- strict:

  Check the declared number of peaks against the parsed pairs.

- max_spectra:

  Optional import limit, useful for previews.

- path:

  Existing spectral-library SQLite file.

- x:

  A GCMSlibrary object.

- ...:

  Additional printing arguments.

## Value

A GCMSlibrary object. Memory objects contain sparse spectra; disk
objects retain a path and open read-only connections when queried.

## Details

MSP tags are case-insensitive. Repeated tags, synonyms and unknown
fields remain available through libraryMetadata(raw=TRUE). Raw masses
and intensities are retained as double precision arrays. Sparse matrix
extraction sums duplicate masses; nominal binning and normalization are
explicit options. Explicit Pseudocount/Ion_total tags retain a supplied
full-spectrum total; librarySearch uses it for count-sensitive
comparisons without changing raw peaks. In-memory libraries store a
dgCMatrix plus compact metadata. Exceptional unsorted, duplicate or
zero-intensity peak records retain exact raw arrays.

SQLite stores one pair of little-endian float64 blobs per spectrum,
indexed metadata and a long retention-index table. Import is
transactional and uses bounded batches. Library files can be moved and
reopened with openLibrary(). Existing databases are never overwritten.
No live connection is serialized.

Retention indices retain their raw field, phase class and provenance.
NIST/OpenChrom comment summaries SemiStdNP, StdNP and StdPolar map to
DB5, DB1 and DBWAX; slash-delimited uncertainty fields remain in the raw
value. DB1, DB5 and DBWAX denote nonpolar, semi-standard nonpolar and
polar classes; a separate column/stationary-phase field also supplies
the RI phase class. An unspecified column remains unknown. Zero RI
placeholders are omitted from the index table. Supplied SMILES are
preserved without chemical canonicalization.

## Examples

``` r
f <- tempfile(fileext=".msp")
writeLines(c("Name: example", "Formula: C2H6O", "RI_DB5: 560",
             "Num Peaks: 3", "31 100; 45 60; 46 20;"), f)
lib <- readMSP(f)
lib
libraryMetadata(lib)
retentionIndices(lib)
spectra(lib)
```
