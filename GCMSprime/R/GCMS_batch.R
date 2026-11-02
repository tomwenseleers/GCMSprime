GCMS_cdf_metadata <- function(path) {
  nc<-ncdf4::nc_open(path);on.exit(ncdf4::nc_close(nc))
  attr<-function(key,var=0) {
    a<-ncdf4::ncatt_get(nc,var,key);if(isTRUE(a$hasatt))a$value else NULL
  }
  date<-as.POSIXct(NA_real_,origin="1970-01-01",tz="UTC")
  for(key in c("experiment_date_time_stamp","experiment_start_date_time","run_start_time","acquisition_date_time")) {
    value<-attr(key)
    if(is.null(value))next
    for(format in c("%Y%m%d%H%M%S%z","%Y-%m-%dT%H:%M:%S%z","%Y-%m-%d %H:%M:%S")) {
      candidate<-as.POSIXct(strptime(as.character(value),format,tz="UTC"))
      if(length(candidate)==1L&&!is.na(candidate)){date<-candidate;break}
    }
    if(!is.na(date))break
  }
  value<-function(key)if(key%in%names(nc$var))as.numeric(ncdf4::ncvar_get(nc,key))else NULL
  list(acquired_at=date,scan_duration=value("scan_duration"),
    intensity_units=attr("units","intensity_values"),
    mass_min=value("mass_range_min"),mass_max=value("mass_range_max"))
}

GCMS_batch_files <- function(directory,output_dir,sample_order=NULL,alkanes_file=NULL) {
  directory<-normalizePath(directory,winslash="/",mustWork=TRUE)
  paths<-list.files(directory,pattern="[.]cdf$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
  paths<-normalizePath(paths,winslash="/",mustWork=TRUE)
  out<-normalizePath(output_dir,winslash="/",mustWork=FALSE)
  paths<-paths[!startsWith(tolower(paths),paste0(tolower(out),"/"))]
  relative<-substring(paths,nchar(directory)+2L)
  is_standard<-grepl("(^|/)standards(/|$)",relative,ignore.case=TRUE)
  if(!is.null(alkanes_file))is_standard<-is_standard|tolower(paths)==tolower(normalizePath(alkanes_file,winslash="/",mustWork=TRUE))
  standards<-paths[is_standard];paths<-paths[!is_standard];relative<-relative[!is_standard]
  if(length(paths)<2L)stop("Directory input needs at least two sample CDFs outside standards/.")
  metadata<-lapply(paths,GCMS_cdf_metadata)
  dates<-vapply(metadata,function(x)as.numeric(x$acquired_at),0)
  if(!is.null(sample_order)) {
    index<-if(is.numeric(sample_order))sample_order else match(sample_order,relative)
    if(length(index)!=length(paths)||anyNA(index)||!setequal(index,seq_along(paths)))stop("sample_order must list every relative sample path once.")
  } else {
    if(anyNA(dates))stop("Acquisition dates missing: supply batch$sample_order with relative CDF paths. Export dates are not acquisition dates.")
    index<-order(dates,relative)
  }
  paths<-paths[index];relative<-relative[index];metadata<-metadata[index];dates<-dates[index]
  group<-dirname(relative);group[group=="."]<-basename(directory)
  name<-tools::file_path_sans_ext(basename(paths))
  id<-make.unique(gsub("[^[:alnum:]_.-]","_",paste(group,name,sep=".")))
  selected<-rep(NA_character_,length(paths))
  if(!is.null(alkanes_file))selected[]<-normalizePath(alkanes_file,winslash="/",mustWork=TRUE) else if(length(standards)) {
    st<-vapply(standards,function(p)as.numeric(GCMS_cdf_metadata(p)$acquired_at),0)
    if(length(standards)>1L&&(anyNA(st)||anyNA(dates)))stop("Multiple standards need acquisition dates or an explicit alkanes_file.")
    for(s in seq_along(paths))selected[s]<-if(length(standards)==1L)standards else standards[which.min(abs(st-dates[s]))]
  }
  list(manifest=data.frame(sample_id=id,sample_name=name,group=group,cdf_file=paths,
    acquired_at=as.POSIXct(dates,origin="1970-01-01",tz="UTC"),alkanes_file=selected,
    stringsAsFactors=FALSE),metadata=metadata,standards=standards)
}

#' Multisample GC/MS deconvolution, alignment and shared spectral recovery
#'
#' Supply a directory to [GCMSdeco()] to process an experiment. Sample group
#' subdirectories are retained, and `standards/` contains alkane calibrations.
#' The closest standard in acquisition time calibrates each sample. Results
#' share spectra through [L0nmf::coupledNMF()] while keeping original scan grids.
#'
#' @name GCMSbatch
#' @param x,object A multisample result returned by [GCMSdeco()].
#' @param sample One sample index or identifier. `NULL` in `coef()` returns the
#'   common spectral library. `NULL` in `plot()` plots the coupled area table.
#' @param normalize Spectrum normalization: unnormalized or unit L1.
#' @param ... Arguments passed to the corresponding sample or coupled method.
#' @details
#' Directory layout can be `experiment/standards/alkanes.cdf`,
#' `experiment/beer1/sample001.cdf`, and `experiment/beer2/sample002.cdf`.
#' Analytical samples are sorted by acquisition timestamp, including its time
#' zone. The exported CDF timestamp is never used as an acquisition date.
#' If acquisition dates are absent, pass explicit relative sample paths in
#' `batch=list(sample_order=...)`.
#'
#' `batch` supports `method="cosine"` or `"hellinger"`, `align=TRUE`,
#' `coupling=list()` for [L0nmf::coupledNMF()] controls, `sample_order=NULL`,
#' `exposure=NULL`, `resume=FALSE` and `plots=FALSE`. `exposure` may be a positive
#' scalar, a list with one scan-wise vector per acquisition-ordered sample, or
#' a function of a sample's metadata and time vector. The default preserves
#' vendor intensities. Known integration durations are appropriate when the
#' data are integrated ion counts. Sequential full-scan mass spectrometers
#' allocate acquisition time over m/z; exported peak counts and retained traces
#' do not identify the per-mass dwell time. No automatic division by the number
#' of exported channels is applied.
#'
#' Each sample has its own initial PSF and smooth background. Original-channel
#' means and the likelihood are used for coupled updates. Shared spectra are
#' exported on the union of recorded mass ranges; unmeasured channels receive
#' zero likelihood weight. Missing compounds retain `NA` areas. No gap filling
#' is performed. `coupled=FALSE` keeps the independently fitted results.
#'
#' The output directory contains per-sample MSP, Excel and fitted RDS files,
#' an acquisition-ordered `samples.csv`, `areas.csv`, a consensus MSP library,
#' and `experiment.xlsx` with samples, compounds, areas and packed spectra.
#' `batch=list(plots=TRUE)` also writes each sample's trace and heatmap PDFs.
#' Observations and factors are stored per sample, enabling streamed fitting;
#' the returned object contains paths, tables and compact common spectra.
#' Keep the output directory with the saved object. Explicit `resume=TRUE`
#' reuses independent fits only when the file metadata and fitting arguments
#' match; coupled refinement is recomputed.
#'
#' @return A `GCMSbatch` object with `samples`, optional `coupled` fit, sample
#'   RDS `paths`, consensus `areas`, output files and stage timings. Use
#'   `readRDS(x$paths[i])` for a full per-sample [GCMSdeco()] object.
#' @examples
#' \dontrun{
#' experiment <- GCMSdeco("chromatograms", output_dir="results",
#'   family=quasipoisson("identity"), precluster=32L,
#'   batch=list(method="hellinger",plots=TRUE))
#' experiment$areas
#' plot(experiment,sample=1,type="lines")
#' }
#' @seealso [GCMSdeco()], [L0nmf::coupledNMF()], [write_GCMS_cdf()]
#' @md
NULL

GCMSdeco_directory <- function(cdf_file,msp_file=NULL,family,lambda,nwindows,
    alkanes_file,time_range,calibration_range,baseline_df,precluster,dispersion,theta,
    correlation,ncores,control,psf_control,verbose,xlsx_file,sample_name,library,search,
    output_dir,coupled,batch) {
  begun<-proc.time()[["elapsed"]]
  if(!requireNamespace("ncdf4",quietly=TRUE)||!requireNamespace("openxlsx",quietly=TRUE))
    stop("Install ncdf4 and openxlsx for directory input and batch exports.")
  if(!is.null(msp_file)&&(!is.character(msp_file)||length(msp_file)!=1L||is.na(msp_file)||
      !grepl("[.]msp$",msp_file,ignore.case=TRUE)))stop("msp_file must end in .msp.")
  if(!is.null(xlsx_file))GCMS_xlsx_path(xlsx_file)
  defaults<-list(method="cosine",align=TRUE,coupling=list(),sample_order=NULL,
    exposure=NULL,resume=FALSE,plots=FALSE)
  if(!is.list(batch)||length(batch)&&(is.null(names(batch))||anyDuplicated(names(batch))||
      length(setdiff(names(batch),names(defaults)))))stop("Use named batch options; see ?GCMSbatch.")
  batch<-utils::modifyList(defaults,batch,keep.null=TRUE);ncores<-GCMS_threads(ncores)
  if(!is.logical(coupled)||length(coupled)!=1L||is.na(coupled))stop("coupled must be TRUE or FALSE.")
  if(is.null(output_dir))output_dir<-paste0(normalizePath(cdf_file,winslash="/",mustWork=TRUE),"_GCMSprime")
  found<-GCMS_batch_files(cdf_file,output_dir,batch$sample_order,alkanes_file)
  manifest<-found$manifest;metadata<-found$metadata
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE);output_dir<-normalizePath(output_dir,winslash="/")
  dir.create(file.path(output_dir,"samples"),showWarnings=FALSE)
  descriptors<-vector("list",nrow(manifest));paths<-character(nrow(manifest));mass_axes<-vector("list",length(paths))
  cache_hit<-logical(length(paths))
  initial_time<-system.time(for(s in seq_along(paths)) {
    id<-manifest$sample_id[s];d<-file.path(output_dir,"samples",id);dir.create(d,showWarnings=FALSE)
    cache<-file.path(d,"initial.rds");file<-manifest$cdf_file[s]
    standard<-manifest$alkanes_file[s];if(is.na(standard))standard<-NULL
    args<-list(cdf_file=file,msp_file=file.path(d,"initial.msp"),family=family,lambda=lambda,
      nwindows=nwindows,alkanes_file=standard,time_range=time_range,calibration_range=calibration_range,
      baseline_df=baseline_df,precluster=precluster,dispersion=dispersion,theta=theta,
      correlation=correlation,ncores=ncores,control=control,psf_control=psf_control,verbose=verbose,
      sample_name=id)
    signature<-GCMS_batch_signature(file,standard,args)
    saved<-if(isTRUE(batch$resume)&&file.exists(cache))readRDS(cache)else NULL
    cache_hit[s]<-!is.null(saved)&&identical(saved$signature,signature)
    fit<-if(cache_hit[s])saved$fit else do.call(GCMSdeco,args)
    saveRDS(list(signature=signature,fit=fit),cache)
    paths[s]<-file.path(d,"fit.rds");saveRDS(fit,paths[s])
    raw<-GCMS_read_cdf(file,time_range);yfile<-file.path(d,"observed.rds");saveRDS(raw$Y,yfile,compress=FALSE)
    initfile<-file.path(d,"initial_model.rds");saveRDS(fit$model,initfile,compress=FALSE)
    scan_meta<-metadata[[s]]
    for(key in c("scan_duration","mass_min","mass_max"))
      if(length(scan_meta[[key]])>=max(fit$scans))scan_meta[[key]]<-scan_meta[[key]][fit$scans]
    exposure<-if(is.null(batch$exposure))1 else if(is.function(batch$exposure))batch$exposure(scan_meta,fit$time)else
      if(is.list(batch$exposure))batch$exposure[[s]]else batch$exposure
    if(!is.null(batch$exposure)&&is.null(exposure))stop("The requested exposure is missing for sample ",id,".")
    small<-fit$psf_width_scans*sqrt(log(100)/log(2))*median(diff(fit$time))
    descriptors[[s]]<-list(y=yfile,time=fit$time,fit=initfile,exposure=exposure,
      sample_id=id,group=manifest$group[s],small_width=small)
    mass_axes[[s]]<-fit$mass
    if(verbose)cat("Sample",s,"of",length(paths),id,":",fit$model$rank,"initial compounds\n")
  })[["elapsed"]]
  axis<-sort(unique(unlist(mass_axes)));common<-NULL
  if(coupled) {
    # Pad only unmeasured mass ranges; absent measured ions remain genuine zeros.
    for(s in seq_along(paths)) {
      x<-descriptors[[s]];mass<-mass_axes[[s]];n<-length(x$time)
      if(!identical(mass,axis)) {
        ix<-match(mass,axis);Y<-matrix(0,n,length(axis));Y[,ix]<-readRDS(x$y);colnames(Y)<-as.character(axis)
        saveRDS(Y,x$y,compress=FALSE);f<-readRDS(x$fit)
        S<-matrix(0,nrow(f$S),length(axis));S[,ix]<-f$S;f$S<-S
        A<-matrix(0,nrow(f$alpha),length(axis));A[,ix]<-f$alpha;f$alpha<-A
        saveRDS(f,x$fit,compress=FALSE)
        x$weights<-file.path(dirname(x$y),"channel_weights.rds")
        saveRDS(matrix(rep(as.numeric(axis%in%mass),each=n),n,length(axis)),x$weights,compress=FALSE)
      }
      descriptors[[s]]<-x
    }
    common<-L0nmf::coupledNMF(descriptors,family=family,lambda=lambda,
      path=tempfile("coupled_",tmpdir=output_dir),precluster=precluster,
      method=batch$method,align=batch$align,ncores=ncores,control=batch$coupling,verbose=verbose)
    for(s in seq_along(paths)) {
      fit<-readRDS(paths[s]);fit$model<-L0nmf::sampleNMF(common,s);fit$mass<-axis
      fit$compounds<-GCMS_compounds(fit$model,fit$time,fit$sample_name)
      fit$compounds$coupled_compound<-fit$model$compound_id[fit$compounds$compound]
      fit$compounds$integrated_area<-fit$model$integrated_area[fit$compounds$compound]
      fit$compounds$apex_scan<-fit$scans[fit$compounds$apex_scan]
      saveRDS(fit,paths[s])
    }
  }
  for(s in seq_along(paths)) {
    fit<-readRDS(paths[s]);d<-dirname(paths[s]);id<-manifest$sample_id[s]
    fit$msp_file<-file.path(d,paste0(id,".msp"));fit$xlsx_file<-file.path(d,paste0(id,".xlsx"))
    if(!is.null(library))fit<-do.call(identify_compounds,c(list(x=fit,library=library),search))
    fit$compounds<-GCMS_write_msp(fit$model,fit$time,fit$mass,fit$msp_file,fit$compounds)
    GCMS_write_xlsx(fit$model,fit$time,fit$mass,fit$compounds,fit$xlsx_file,
      fit$cdf_file,fit$alkanes_file,fit$seconds,
      matches=if(is.null(fit$identification))NULL else fit$identification$matches,
      calibration=if(is.null(fit$identification))NULL else fit$identification$calibration,
      search_settings=if(is.null(fit$identification))NULL else fit$identification$settings)
    if(isTRUE(batch$plots)) {
      plot(fit,type="lines",file=file.path(d,paste0(id,"_traces.pdf")))
      plot(fit,type="heatmap",file=file.path(d,paste0(id,"_heatmaps.pdf")))
    }
    saveRDS(fit,paths[s])
    # Observations now live in the coupled store; plotting reloads the CDF.
    unlink(c(descriptors[[s]]$y,descriptors[[s]]$fit,
      if(is.character(descriptors[[s]]$weights))descriptors[[s]]$weights else NULL))
  }
  manifest$fit_file<-paths
  manifest$cache_hit<-cache_hit
  manifest$compounds<-vapply(paths,function(p)readRDS(p)$model$rank,0L)
  manifest$GIC<-vapply(paths,function(p)GIC(readRDS(p)),0)
  areas<-if(is.null(common))NULL else common$areas
  object<-structure(list(samples=manifest,paths=paths,coupled=common,areas=areas,
    standards=found$standards,output_dir=output_dir,
    seconds=c(initial=initial_time,total=proc.time()[["elapsed"]]-begun)),class="GCMSbatch")
  utils::write.csv(manifest,file.path(output_dir,"samples.csv"),row.names=FALSE)
  if(!is.null(areas)) {
    utils::write.csv(cbind(sample_id=rownames(areas),as.data.frame(areas)),file.path(output_dir,"areas.csv"),row.names=FALSE)
    GCMS_batch_consensus(object,msp_file,xlsx_file)
  }
  saveRDS(object,file.path(output_dir,"experiment.rds"));object
}

GCMS_batch_signature <- function(file,standard,args) {
  f<-args$family;if(is.function(f))f<-f(link="identity")
  # Family closures have distinct R environments on each construction.
  # Compare their model settings rather than environment identity.
  if(is.list(f))args$family<-list(family=f$family,link=f$link,scalar=f$scalar,
    theta=if(is.function(f$getTheta))f$getTheta(TRUE)else NULL,n.theta=f$n.theta)
  info<-file.info(file)
  list(size=info$size,mtime=as.numeric(info$mtime),args=args,
    standard=if(is.null(standard))NULL else list(size=file.info(standard)$size,
      mtime=as.numeric(file.info(standard)$mtime)),
    package=as.character(utils::packageVersion("GCMSprime")))
}

GCMS_batch_consensus <- function(object,msp_file=NULL,xlsx_file=NULL) {
  common<-object$coupled;G<-nrow(common$P)
  spectra<-coef(common);axis<-as.numeric(common$mass)
  if(is.null(msp_file))msp_file<-file.path(object$output_dir,"consensus.msp")
  if(is.null(xlsx_file))xlsx_file<-file.path(object$output_dir,"experiment.xlsx")
  occurrence<-common$compounds
  rt<-vapply(seq_len(G),function(g)weighted.mean(occurrence$anchor[occurrence$root==g],occurrence$area[occurrence$root==g]),0)
  table<-data.frame(sample=rep("consensus",G),compound=seq_len(G),name=paste0("consensus_",seq_len(G)),
    retention_time=rt,area=rowSums(spectra),n_samples=colSums(!is.na(common$areas)),stringsAsFactors=FALSE)
  GCMS_write_msp(list(H=Matrix::Diagonal(G),S=spectra),rt,axis,msp_file,table)
  if(!requireNamespace("openxlsx",quietly=TRUE))stop("Install openxlsx for batch Excel output.")
  table$Spectrum<-GCMS_spectrum_text(spectra,axis)
  wb<-openxlsx::createWorkbook()
  sheets<-list(samples=object$samples,compounds=table,
    areas=cbind(sample_id=rownames(common$areas),as.data.frame(common$areas)))
  for(name in names(sheets)) {
    openxlsx::addWorksheet(wb,name);openxlsx::writeData(wb,name,sheets[[name]])
    openxlsx::addStyle(wb,name,openxlsx::createStyle(wrapText=FALSE),rows=seq_len(nrow(sheets[[name]])+1L),
      cols=seq_len(ncol(sheets[[name]])),gridExpand=TRUE)
    openxlsx::freezePane(wb,name,firstRow=TRUE,firstCol=TRUE)
  }
  openxlsx::saveWorkbook(wb,xlsx_file,overwrite=TRUE)
  invisible(object)
}

GCMS_batch_sample <- function(x,sample) {
  index<-if(is.character(sample))match(sample,x$samples$sample_id)else sample
  if(length(index)!=1L||is.na(index)||index<1||index>length(x$paths)||index!=floor(index))stop("Choose one sample by index or sample_id.")
  readRDS(x$paths[index])
}
#' @rdname GCMSbatch
#' @export
fitted.GCMSbatch <- function(object,sample=1L,...)fitted(GCMS_batch_sample(object,sample),...)
#' @rdname GCMSbatch
#' @export
coef.GCMSbatch <- function(object,sample=NULL,normalize=c("none","l1"),...) {
  normalize<-match.arg(normalize)
  if(is.null(sample)) {
    if(is.null(object$coupled))stop("Select a sample for independent fits.")
    return(coef(object$coupled,normalize=normalize,...))
  }
  coef(GCMS_batch_sample(object,sample),normalize=normalize,...)
}
#' @rdname GCMSbatch
#' @export
plot.GCMSbatch <- function(x,sample=NULL,...) {
  if(is.null(sample)) {
    if(is.null(x$coupled))stop("Select a sample for independent fits.")
    return(plot(x$coupled,...))
  }
  plot(GCMS_batch_sample(x,sample),...)
}
#' @rdname GCMSbatch
#' @export
logLik.GCMSbatch <- function(object,...) {
  if(is.null(object$coupled))stop("Select a sample for independent fits.")
  logLik(object$coupled,...)
}
#' @rdname GCMSbatch
#' @export
GIC.GCMSbatch <- function(object,...) {
  if(is.null(object$coupled))stop("Select a sample for independent fits.")
  GIC(object$coupled,...)
}
#' @rdname GCMSbatch
#' @param k Information-criterion multiplier.
#' @export
AIC.GCMSbatch <- function(object,...,k=2) {
  if(is.null(object$coupled))stop("Select a sample for independent fits.")
  AIC(object$coupled,...,k=k)
}
#' @rdname GCMSbatch
#' @export
BIC.GCMSbatch <- function(object,...)BIC(object$coupled,...)
#' @rdname GCMSbatch
#' @export
EBIC.GCMSbatch <- function(object,...)EBIC(object$coupled,...)
#' @rdname GCMSbatch
#' @export
MBIC.GCMSbatch <- function(object,...)MBIC(object$coupled,...)
#' @rdname GCMSbatch
#' @export
nobs.GCMSbatch <- function(object,...) {
  if(!is.null(object$coupled))return(nobs(object$coupled,...))
  sum(vapply(object$paths,function(p)nobs(readRDS(p)),0L))
}
#' @rdname GCMSbatch
#' @export
print.GCMSbatch <- function(x,...) {
  cat("GC/MS experiment:",nrow(x$samples),"samples in",length(unique(x$samples$group)),"groups\n")
  if(!is.null(x$coupled))cat(nrow(x$coupled$P),"consensus compounds; GIC",format(GIC(x),digits=9),"\n")
  cat("Results:",x$output_dir,"\n");invisible(x)
}
