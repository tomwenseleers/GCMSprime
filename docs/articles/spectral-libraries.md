# Spectral library import

GCMSprime keeps raw peaks, chemical identifiers and retention-index
provenance together. Small libraries can be used in memory. Large
libraries are imported in batches into SQLite and queried without
loading every spectrum.

``` r

library(GCMSprime)
lib <- read_msp("library.msp", backend="sqlite", db_file="library.sqlite")
lib <- open_library("library.sqlite")
library_metadata(lib, which=1:10)
spectra(lib, which=1:1000, nominal=TRUE, normalize="l2")
retention_indices(lib, which=1:1000, phase="DB5")
```

Matrix rows identify mass values and columns identify library spectra.
Raw extraction with `format="list"` preserves the input mass/intensity
arrays, including duplicate masses. Sparse extraction sums duplicate
masses. Nominal binning and normalization are explicit choices.

## Source formats

- OpenChrom Agilent and Shimadzu MSP exports use the streaming native
  reader. Repeated `Synon` and `Comments` tags and unknown fields remain
  available with `library_metadata(lib, raw=TRUE)`.
- [MoNA](https://mona.fiehnlab.ucdavis.edu/downloads) provides
  NIST-compatible MSP exports, including its GC-MS collection. Those
  exports use
  [`read_msp()`](https://tomwenseleers.github.io/GCMSprime/reference/read_msp.md).
- [Golm](http://gmd.mpimp-golm.mpg.de/download/) MSP exports use the
  same interface. Source-specific annotations remain in the raw
  metadata; column and index conventions should be verified from the
  library’s accompanying metadata.
- [GCMS-ID](https://gcms-id.ca/downloads) provides predicted spectra and
  retention indices as CSV archives. After extraction, use
  `read_gcms_id("NP-Atlas.csv", backend="sqlite", db_file="NP-Atlas.sqlite")`.
  Their predicted origin is recorded explicitly.

## Retention indices and chemical information

Retention indices form a long table. Multiple index values can belong to
the same spectrum, with a phase class, raw column annotation, index
convention and provenance. DB1, DB5 and DBWAX identify nonpolar,
semi-standard nonpolar and polar classes. Unspecified columns retain an
unknown class. A zero RI placeholder is omitted from the index table
while remaining in the raw source metadata.

Names, CAS numbers, formulae, InChI, InChIKey and supplied SMILES have
normalized columns. Repeated synonyms and additional fields are
retained. Each source record remains distinct, including multiple
spectra or RI predictions for a compound. Supplied SMILES are preserved
in their original form.

## Storage and scale

SQLite stores indexed metadata, a separate retention-index table and one
pair of little-endian float64 peak blobs per spectrum. File-backed
objects retain a path and open read-only connections on demand. Imports
use transactions and publish the completed file only after parsing
succeeds.

A nominal-mass `dgCMatrix` requires approximately 12q+4(N+1) bytes for q
nonzero peaks across N spectra, excluding metadata. For one million
spectra with 80 nonzero peaks each this is about 964 MB; a dense matrix
with 1,200 mass rows would require about 9.6 GB. SQLite allows raw
high-resolution peaks and metadata to remain on disk while matrix blocks
are assembled for subsequent analysis.

Library matching and compound identification can build on this storage
layer. The import functions expose both raw peaks and normalized sparse
blocks for those analyses.

## Search recovered spectra

Use [RI calibration and spectral
search](https://tomwenseleers.github.io/GCMSprime/articles/library-search.md)
to identify recovered compounds with column-specific RI constraints,
Stein–Scott similarity and symmetric multinomial distances.
