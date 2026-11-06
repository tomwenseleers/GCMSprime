bench_mshub <- function(files,out,ncores,control,timeout,dots) {
  if (length(files)<10L)
    stop("The MSHub comparison requires at least ten samples.",call.=FALSE)
  if (is.null(control$export)) {
    exe <- executable(control$exe,"python")
    source <- control$source
    if (is.null(source) || !file.exists(file.path(source,"proc/io/importmsdata.py")))
      stop("Supply the standalone MSHub source directory in control.",call.=FALSE)
    prep <- bench_time(stage_inputs(files,file.path(out,"input"),"cdf"))
    source <- normalizePath(source,winslash="/",mustWork=TRUE)
    pyenv <- c(OMP_NUM_THREADS=as.character(ncores),OPENBLAS_NUM_THREADS=as.character(ncores))
    if (!is.null(control$pythonpath)) pyenv <- c(pyenv,PYTHONPATH=control$pythonpath)
    directory <- file.path(out,"export")
    run_logged(exe,c(adapter_asset("run-mshub.py"),source,dirname(prep$value[1L]),directory),
      out,file.path(out,"run.log"),timeout,env=pyenv)
  } else {
    directory <- control$export
    prep <- list(seconds=NA_real_)
  }
  a <- utils::read.delim(file.path(directory,"peaks.tsv"),quote="",comment.char="")
  peaks <- peak_table(a$sample,a$id,a$rt,a$area,group=a$id)
  if (!all(peaks$sample%in%sample_id(files)))
    stop("MSHub sample names do not match the input files.",call.=FALSE)
  ans <- records(peaks,lapply(a$spectrum,numeric_pairs,pattern=","),
    area_kind="MSHub aligned NMF integral",shared=length(unique(peaks$group)))
  ans$rt_basis <- "aligned"
  ans$note <- "Aligned group retention times and spectra are shared across samples."
  ans$analysis <- scan(file.path(directory,"seconds.txt"),quiet=TRUE,nmax=1)
  ans$preparation <- prep$seconds
  ans
}
