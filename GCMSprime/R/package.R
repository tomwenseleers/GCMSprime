#' Scalable and Automated GC–MS Deconvolution, Quantification, and Spectral Library Search in R
#'
#' @keywords internal
#' @import L0glm
#' @import L0nmf
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
