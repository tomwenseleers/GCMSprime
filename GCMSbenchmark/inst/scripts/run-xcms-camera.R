# CAMERA's legacy xcmsSet interface uses functions from attached dependencies.
# Keep those attachments inside this isolated R session and its workers.
load_tools <- function(libraries) {
  .libPaths(libraries)
  suppressPackageStartupMessages({
    library(xcms)
    library(CAMERA)
  })
  invisible(NULL)
}

xcms_detection <- function(file,clock,options=NULL) {
  params <- list(method="matchedFilter",fwhm=6*stats::median(diff(clock)),
    step=1,steps=1,max=1000L,snthresh=10)
  if (!is.null(options)) params <- utils::modifyList(params,options)
  do.call(xcms::xcmsSet,c(list(files=file,BPPARAM=BiocParallel::SerialParam()),params))
}

camera_fwhm <- function(xs,options=NULL) {
  if (!nrow(xs@peaks)) return(NULL)
  do.call(CAMERA::groupFWHM,c(list(object=CAMERA::xsAnnotate(xs,sample=1L)),options))
}

camera_correlation <- function(groups,raw,options=NULL) {
  if (is.null(groups) || identical(options,FALSE)) return(groups)
  params_corr <- list(cor_eic_th=.75,calcIso=FALSE,calcCiS=TRUE,calcCaS=FALSE)
  if (!is.null(options)) params_corr <- utils::modifyList(params_corr,options)
  do.call(CAMERA::groupCorr,c(list(object=groups,xraw=raw),params_corr))
}

camera_components <- function(groups,raw,min_ions,bin_width=1) {
  empty <- list(peaks=data.frame(id=integer(),rt=numeric(),area=numeric(),height=numeric()),
    spectra=list(),profiles=list())
  if (is.null(groups)) return(empty)
  ions <- groups@xcmsSet@peaks
  components <- lapply(groups@pspectra,function(i) {
    p <- ions[i,,drop=FALSE]
    p <- p[is.finite(p[,"into"]) & p[,"into"]>0,,drop=FALSE]
    p <- p[order(p[,"into"],decreasing=TRUE),,drop=FALSE]
    p[!duplicated(round(p[,"mz"])),,drop=FALSE]
  })
  components <- Filter(function(p) nrow(p)>=min_ions,components)
  if (!length(components)) return(empty)
  peaks <- data.frame(id=seq_along(components),
    rt=vapply(components,function(p) p[1L,"rt"]/60,0),
    area=vapply(components,function(p) sum(p[,"into"]),0),
    height=vapply(components,function(p) p[1L,"maxo"],0))
  spectra <- lapply(components,function(p) unname(p[,c("mz","into"),drop=FALSE]))
  profiles <- lapply(components,function(p) {
    base <- p[1L,]
    eic <- xcms::rawEIC(raw,mzrange=base[c("mzmin","mzmax")]+c(-.499,.499)*bin_width,
      rtrange=base[c("rtmin","rtmax")])
    cbind(raw@scantime[eic$scan]/60,eic$intensity)
  })
  list(peaks=peaks,spectra=spectra,profiles=profiles)
}

analyse <- function(file,config) {
  cat("Input:",basename(file),"\n")
  raw <- xcms::xcmsRaw(file,profstep=0)
  xs <- if (length(raw@scantime)>=3L) xcms_detection(file,raw@scantime,config$detection)
  coarse <- if (!is.null(xs)) camera_fwhm(xs,config$fwhm_grouping)
  groups <- camera_correlation(coarse,raw,config$grouping)
  step <- if (is.null(config$detection$step)) 1 else config$detection$step
  camera_components(groups,raw,config$min_ions,step)
}

run <- function(config) {
  if (config$workers==1L) {
    load_tools(config$libraries)
    return(lapply(config$files,analyse,config=config))
  }
  cluster <- parallel::makePSOCKcluster(config$workers,outfile="")
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster,load_tools,config$libraries)
  parallel::clusterExport(cluster,c("xcms_detection","camera_fwhm","camera_correlation",
    "camera_components"),envir=environment(run))
  parallel::parLapplyLB(cluster,config$files,analyse,config=config)
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  saveRDS(run(readRDS(args[1L])),args[2L])
}
