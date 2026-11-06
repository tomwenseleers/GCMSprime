#' Scalable and Automated GC–MS Deconvolution, Quantification, and Spectral Library Search in R
#'
#' @details
#' \code{\link{gcmsDeco}} resolves overlapping GC-MS signals from a single
#' CDF file or a directory of related samples. It returns fitted peak profiles,
#' full spectra and integrated areas, with MSP and optional Excel output.
#' Directory fits align samples and report a common compound list.
#'
#' Use \code{\link{calibrateRI}} to estimate retention indices from alkane
#' standards and \code{\link{identifyCompounds}} to rank candidates from
#' a reference library. A ranked match records spectral and retention-index
#' agreement for inspection.
#'
#' @keywords internal
#' @import L0glm
#' @import logconNMF
#' @importFrom Matrix sparseMatrix Diagonal colSums rowSums head tail
#' @import stats graphics grDevices methods
#' @importFrom utils data flush.console
"_PACKAGE"

#' @importFrom L0glm GIC
#' @export
L0glm::GIC
#' @importFrom L0glm EBIC
#' @export
L0glm::EBIC
#' @importFrom L0glm MBIC
#' @export
L0glm::MBIC

#' @useDynLib GCMSprime, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL
