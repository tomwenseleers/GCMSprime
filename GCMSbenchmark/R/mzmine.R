read_mgf <- function(path) {
  x <- trimws(readLines(path,warn=FALSE))
  start <- which(x=="BEGIN IONS")
  stop <- which(x=="END IONS")
  if (length(start)!=length(stop) || any(stop<=start))
    stop("Malformed mzmine MGF export.",call.=FALSE)
  ids <- vapply(seq_along(start),function(i) {
    rows <- x[seq.int(start[i]+1L,stop[i]-1L)]
    id <- sub("^FEATURE_ID=","",grep("^FEATURE_ID=",rows,value=TRUE))
    if (length(id)==1L) id else sub("^SCANS=","",grep("^SCANS=",rows,value=TRUE))[1L]
  },"")
  ans <- lapply(seq_along(start),function(i) {
    rows <- x[seq.int(start[i]+1L,stop[i]-1L)]
    rows <- rows[grepl("^[0-9.]+[[:blank:]]",rows)]
    numeric_pairs(paste(rows,collapse=";"),"[[:blank:]]+")
  })
  if (anyNA(ids) || anyDuplicated(ids)) stop("Ambiguous mzmine MGF identifiers.",call.=FALSE)
  stats::setNames(ans,ids)
}

bench_mzmine <- function(files,out,ncores,control,timeout,dots) {
  exe <- executable(control$exe,"mzmine_console.exe")
  batch <- if (is.null(control$batch))
    system.file("extdata","mzmine-gcei.mzbatch",package="GCMSbenchmark") else control$batch
  if (length(batch)!=1L || !file.exists(batch))
    stop("Supply an exported mzmine GC-EI .mzbatch in control.",call.=FALSE)
  batch <- normalizePath(batch,winslash="/",mustWork=TRUE)
  prep <- bench_time(stage_inputs(files,file.path(out,"input"),"mzXML"))
  input_list <- file.path(out,"files.txt")
  writeLines(normalizePath(prep$value,winslash="\\",mustWork=TRUE),input_list)
  temp <- tempfile("temp-",tmpdir=out)
  dir.create(temp)
  prefix <- file.path(out,"result")
  if (length(Sys.glob(paste0(prefix,"*"))))
    stop("mzmine exports already exist. Use a fresh output directory.",call.=FALSE)
  args <- if (is.null(control$user_file)) character() else
    c("-user",normalizePath(control$user_file,winslash="\\",mustWork=TRUE))
  log <- file.path(out,"run.log")
  run <- bench_time(run_logged(exe,c(args,"-batch",batch,"-input",input_list,
    "-output",prefix,"-threads",as.character(ncores),"-temp",temp),out,log,timeout))
  q <- list.files(out,"_quant\\.csv$",full.names=TRUE,ignore.case=TRUE)
  g <- list.files(out,"\\.mgf$",full.names=TRUE,ignore.case=TRUE)
  if (length(q)!=1L || length(g)!=1L)
    stop("Enable the GC-EI quantification table and spectral MGF in the mzmine batch. See ",
      log,call.=FALSE)
  a <- utils::read.csv(q,check.names=FALSE)
  ms <- read_mgf(g)
  parts <- lapply(seq_along(files),function(i) {
    id <- sample_id(files[i])
    column <- paste0(basename(prep$value[i])," Peak area")
    if (!column %in% names(a))
      stop("mzmine quantification table does not contain ",column,call.=FALSE)
    keep <- is.finite(a[[column]]) & a[[column]]>0
    rt_column <- paste0(basename(prep$value[i])," RT")
    if (!rt_column %in% names(a)) rt_column <- "row retention time"
    if (!rt_column %in% names(a)) stop("No mzmine retention time column.",call.=FALSE)
    ids <- as.character(a[["row ID"]][keep])
    peaks <- peak_table(id,ids,a[[rt_column]][keep],a[[column]][keep],group=ids)
    records(peaks,unname(ms[ids]))
  })
  ans <- combine_records(parts,"mzmine representative-ion area","unavailable")
  if (length(files)>1L) {
    ans$shared <- nrow(a)
    ans$note <- "Aligned group retention times and spectra are shared across samples."
  }
  ans$analysis <- run$seconds
  ans$preparation <- prep$seconds
  ans
}
