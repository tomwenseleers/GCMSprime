search_weights <- function(weights,mass) {
  if(is.null(weights)) return(rep(1,length(mass)))
  if(!is.numeric(weights)||!length(weights)||any(!is.finite(weights)|weights<0))
    stop("ion_weights must be finite nonnegative numbers.")
  if(!is.null(names(weights))) {
    ions<-suppressWarnings(as.numeric(names(weights)))
    if(any(!is.finite(ions)|ions<=0)||anyDuplicated(round(ions))) stop("Use distinct numeric mass names for ion_weights.")
    at<-match(mass,round(ions));w<-rep(1,length(mass));w[!is.na(at)]<-weights[at[!is.na(at)]]
  } else {
    if(!length(weights)%in%c(1L,length(mass))) stop("Unnamed ion_weights must have one value or one per acquired mass bin.")
    w<-rep(weights,length.out=length(mass))
  }
  if(!any(w>0)) stop("At least one acquired ion must have positive weight.")
  as.numeric(w)
}

#' Estimate ion-rarity weights from a spectral library
#'
#' Compute reusable smoothed inverse-document-frequency weights from reference
#' spectra. Rare ions receive larger weights than ubiquitous fragments.
#' @param library A GCMSlibrary.
#' @param mass Acquired nominal mass bins; default 33:550.
#' @param min_relative An ion is present when its intensity is positive and at
#'   least this fraction of the acquired-bin base peak, default 0.01.
#' @param power Exponent for inverse-frequency weighting, default 1. Zero gives
#'   uniform weights; larger values put more emphasis on rare ions.
#' @param ncores Native worker count; NULL uses physical cores.
#' @param batch_size Number of library spectra read per batch.
#' @return A named numeric vector of weights, normalized to mean one. Attributes
#'   frequency, spectra, min_relative, power and seconds record its derivation.
#' @details If L nonempty reference spectra are available and d_j contain ion j,
#'   weights are proportional to [1+log((L+1)/(d_j+1))]^power. Frequencies are
#'   counted across reference spectra, using the requested acquisition bins.
#'   Zero-intensity spectra are omitted. Nominal duplicates are summed before
#'   applying the relative threshold. Source intensities and totals are preserved.
#'
#'   Build weights once and reuse them in librarySearch(ion_weights=...). SQLite
#'   records are read in bounded batches; native frequency counts use OpenMP.
#'   For validation, estimate weights from the reference catalogue after excluding
#'   query records and near-duplicates. Uniform weights remain the search default.
#' @seealso [librarySearch()]
#' @examples
#' f <- tempfile(fileext=".msp")
#' writeLines(c("Name: A", "Num Peaks: 2", "43 100;57 20;",
#'              "Name: B", "Num Peaks: 2", "43 100;91 50;"), f)
#' lib <- readMSP(f)
#' w <- libraryIonWeights(lib, mass=c(43,57,91), ncores=1)
#' w
#' @export
libraryIonWeights <- function(library,mass=33:550,min_relative=.01,power=1,
                                ncores=NULL,batch_size=4096L) {
  begun<-proc.time()[['elapsed']]
  if(!inherits(library,"GCMSlibrary")) stop("Supply a GCMSlibrary.")
  if(!is.numeric(mass)||!length(mass)||any(!is.finite(mass)|mass<=0|mass!=round(mass))||anyDuplicated(mass))
    stop("Use distinct positive nominal mass bins.")
  if(!is.numeric(min_relative)||length(min_relative)!=1L||!is.finite(min_relative)||min_relative<0||min_relative>1||
      !is.numeric(power)||length(power)!=1L||!is.finite(power)||power<0)
    stop("Use min_relative between zero and one and nonnegative power.")
  if(!is.numeric(batch_size)||length(batch_size)!=1L||!is.finite(batch_size)||batch_size<1||
      batch_size!=floor(batch_size)||batch_size>.Machine$integer.max) stop("Use a positive integer batch_size.")
  ncores<-GCMS_threads(ncores);mass<-sort(as.numeric(mass))
  con<-if(library$backend=="sqlite") library_connect(library) else NULL
  if(!is.null(con)) on.exit(DBI::dbDisconnect(con),add=TRUE)
  frequency<-numeric(length(mass));n<-0L
  for(at in seq.int(1L,library$n,by=as.integer(batch_size))) {
    ids<-seq.int(at,min(library$n,at+batch_size-1L))
    a<-search_matrix(library,con,ids,mass)$matrix
    f<-ion_frequency_cpp(a,min_relative,ncores)
    frequency<-frequency+f$frequency;n<-n+f$spectra
  }
  if(!n) stop("No reference ions occur in the requested mass bins.")
  w<-(1+log((n+1)/(frequency+1)))^power;w<-w/mean(w);names(w)<-as.character(mass)
  attr(w,'frequency')<-stats::setNames(frequency,names(w));attr(w,'spectra')<-n
  attr(w,'min_relative')<-min_relative;attr(w,'power')<-power
  attr(w,'seconds')<-proc.time()[['elapsed']]-begun;w
}
