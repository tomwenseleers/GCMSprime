# Scalable and Automated GC–MS Deconvolution, Quantification, and Spectral Library Search in R

Chromatogram deconvolution, compound spectra, integrated areas and
reproducible reports, built on sparse log-concave matrix factorization.
Spectral library tools preserve chemical metadata and retention-index
provenance.

## Details

[`gcmsDeco`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
resolves overlapping GC-MS signals from a single CDF file or a directory
of related samples. It returns fitted peak profiles, full spectra and
integrated areas, with MSP and optional Excel output. Directory fits
align samples and report a common compound list.

Use
[`calibrateRI`](https://tomwenseleers.github.io/GCMSprime/reference/calibrateRI.md)
to estimate retention indices from alkane standards and
[`identifyCompounds`](https://tomwenseleers.github.io/GCMSprime/reference/identifyCompounds.md)
to rank candidates from a reference library. A ranked match records
spectral and retention-index agreement for inspection.

## See also

Useful links:

- <https://tomwenseleers.github.io/GCMSprime/>

- <https://github.com/tomwenseleers/GCMSprime>

- Report bugs at <https://github.com/tomwenseleers/GCMSprime/issues>

## Author

**Maintainer**: Tom Wenseleers <tom.wenseleers@kuleuven.be>

Authors:

- Tom Wenseleers <tom.wenseleers@kuleuven.be>

- Christophe Vanderaa <christophe.vanderaa@kuleuven.be>
