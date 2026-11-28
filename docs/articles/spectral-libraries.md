# Spectral library import

GCMSprime can use **Agilent ChemStation (`.L`) and Shimadzu MS Database
(`.spc`)** libraries, including vendor-formatted editions of NIST. First
export the library from OpenChrom as a NIST text file (`.msp`), then
import that file with
[`readMSP()`](https://tomwenseleers.github.io/GCMSprime/reference/readMSP.md).
If you already have an `.msp` file, go straight to [Import the MSP file
in R](#import-the-msp-file-in-r).

## Convert an Agilent or Shimadzu library with OpenChrom

Install OpenChrom from the [OpenChrom download
page](https://www.openchrom.net/download) or use the command for your
operating system.

| Operating system    | Install command                                   |
|---------------------|---------------------------------------------------|
| Windows (WinGet)    | `winget install --id Lablicate.OpenChrom --exact` |
| Linux with Snap     | `sudo snap install openchrom`                     |
| macOS with Homebrew | `brew install --cask openchrom`                   |

1.  Start OpenChrom. When the **Vendor plugins missing** prompt asks
    whether to install converter plug-ins, choose **Yes**.
2.  Select **Agilent ChemStation GC-MS** for an Agilent library or
    **Shimadzu MS Database** for a Shimadzu library. You can select
    both. Install the selected converters and restart OpenChrom.
3.  Open **Data Analysis**. In its file browser, navigate to the library
    and double-click the Agilent `.L` library (which may appear as a
    folder) or the Shimadzu `.spc` database file.
4.  With the library open, choose **File → Save As…**. Select the NIST
    text (`*.msp`) output format and save the file, for example as
    `library.msp`.

## Import the MSP file in R

For a full NIST library, run
[`readMSP()`](https://tomwenseleers.github.io/GCMSprime/reference/readMSP.md)
once to create a SQLite database. GCMSprime imports spectra in batches
and keeps their raw peaks, chemical identifiers and retention-index
annotations together.

``` r

library(GCMSprime)
lib <- readMSP("library.msp", backend = "sqlite", db_file = "library.sqlite")
```

In later sessions, reopen the database and inspect a spectrum and its
available retention indices.

``` r

lib <- openLibrary("library.sqlite")
libraryMetadata(lib, which = 1L)
spectra(lib, which = 1L, nominal = TRUE, normalize = "l2")
retentionIndices(lib, which = 1L)
```

For a small MSP file, `readMSP("library.msp", backend = "memory")` keeps
the library in memory. Matrix rows identify mass values and columns
identify library spectra. Raw extraction with `format="list"` preserves
the input mass/intensity arrays, including duplicate masses. Sparse
extraction sums duplicate masses. Nominal binning and normalization are
explicit choices.

## Source formats

- For MSP files exported by OpenChrom, repeated `Synon` and `Comments`
  tags and unknown fields remain available with
  `libraryMetadata(lib, raw=TRUE)`.
- [MoNA](https://mona.fiehnlab.ucdavis.edu/downloads) provides
  NIST-compatible MSP exports, including its GC-MS collection. Those
  exports use
  [`readMSP()`](https://tomwenseleers.github.io/GCMSprime/reference/readMSP.md).
- [Golm](http://gmd.mpimp-golm.mpg.de/download/) MSP exports use the
  same interface. Source-specific annotations remain in the raw
  metadata; column and index conventions should be verified from the
  library’s accompanying metadata.
- [GCMS-ID](https://gcms-id.ca/downloads) provides predicted spectra and
  retention indices as CSV archives. After extraction, use
  `readGCMSID("NP-Atlas.csv", backend="sqlite", db_file="NP-Atlas.sqlite")`.
  Their predicted origin is recorded explicitly.

## Retention indices and chemical information

Agilent and Shimadzu editions of the same NIST release can carry
different retention indices for the same reference spectrum. Inspect the
exported phase annotations before an RI-constrained search. For example,
`retentionIndices(lib, which = 1L, phase = "DBWAX")` returns the
polar-column indices recorded for the first spectrum.

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
