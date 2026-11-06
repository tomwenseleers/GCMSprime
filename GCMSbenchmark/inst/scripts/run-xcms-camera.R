args <- commandArgs(TRUE)
config <- readRDS(args[1L])

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

analyse <- function(file,config) {
  cat("Input:",basename(file),"\n")
  raw <- xcms::xcmsRaw(file,profstep=0)
  clock <- raw@scantime
  empty <- list(peaks=data.frame(id=integer(),rt=numeric(),area=numeric(),height=numeric()),
    spectra=list(),profiles=list())
  if (length(clock)<3L) return(empty)
  params <- list(method="matchedFilter",fwhm=6*stats::median(diff(clock)),
    step=1,steps=1,max=1000L,snthresh=10)
  if (!is.null(config$detection)) params <- utils::modifyList(params,config$detection)
  xs <- do.call(xcms::xcmsSet,c(list(files=file,BPPARAM=BiocParallel::SerialParam()),params))
  if (!nrow(xs@peaks)) return(empty)
  groups <- CAMERA::groupFWHM(CAMERA::xsAnnotate(xs,sample=1L))
  params_corr <- list(cor_eic_th=.75,calcIso=FALSE,calcCiS=TRUE,calcCaS=FALSE)
  if (!is.null(config$grouping)) params_corr <- utils::modifyList(params_corr,config$grouping)
  groups <- do.call(CAMERA::groupCorr,c(list(object=groups,xraw=raw),params_corr))
  ions <- groups@xcmsSet@peaks
  components <- lapply(groups@pspectra,function(i) {
    p <- ions[i,,drop=FALSE]
    p <- p[is.finite(p[,"into"]) & p[,"into"]>0,,drop=FALSE]
    p <- p[order(p[,"into"],decreasing=TRUE),,drop=FALSE]
    p[!duplicated(round(p[,"mz"])),,drop=FALSE]
  })
  components <- Filter(function(p) nrow(p)>=config$min_ions,components)
  if (!length(components)) return(empty)
  peaks <- data.frame(id=seq_along(components),
    rt=vapply(components,function(p) p[1L,"rt"]/60,0),
    area=vapply(components,function(p) sum(p[,"into"]),0),
    height=vapply(components,function(p) p[1L,"maxo"],0))
  spectra <- lapply(components,function(p) unname(p[,c("mz","into"),drop=FALSE]))
  profiles <- lapply(components,function(p) {
    base <- p[1L,]
    eic <- xcms::rawEIC(raw,mzrange=base[c("mzmin","mzmax")]+c(-.499,.499)*params$step,
      rtrange=base[c("rtmin","rtmax")])
    cbind(clock[eic$scan]/60,eic$intensity)
  })
  list(peaks=peaks,spectra=spectra,profiles=profiles)
}

run <- function(config) {
  if (config$workers==1L) {
    load_tools(config$libraries)
    return(lapply(config$files,analyse,config=config))
  }
  cluster <- parallel::makePSOCKcluster(config$workers,outfile="")
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster,load_tools,config$libraries)
  parallel::parLapplyLB(cluster,config$files,analyse,config=config)
}

saveRDS(run(config),args[2L])
