GCMS_threads <- function(nthreads) {
  if (is.null(nthreads)) nthreads <- parallel::detectCores(logical=FALSE)
  if (length(nthreads)==1L && is.na(nthreads)) nthreads <- 1L
  if (!is.numeric(nthreads)||length(nthreads)!=1L||!is.finite(nthreads)||
      nthreads<1||nthreads!=floor(nthreads)||nthreads>.Machine$integer.max)
    stop("nthreads must be a positive integer.",call.=FALSE)
  as.integer(nthreads)
}
