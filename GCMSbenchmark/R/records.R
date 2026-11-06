sample_id <- function(path) tools::file_path_sans_ext(basename(path))

peak_table <- function(sample=character(), id=character(), rt=numeric(),
                       area=NA_real_, height=NA_real_, group=NA_character_) {
  n <- length(rt)
  data.frame(sample=rep_len(sample,n), id=rep_len(as.character(id),n), rt=rt,
    area=rep_len(as.numeric(area),n), height=rep_len(as.numeric(height),n),
    group=rep_len(as.character(group),n), stringsAsFactors=FALSE)
}

records <- function(peaks=peak_table(), spectra=vector("list",nrow(peaks)),
                    profiles=vector("list",nrow(peaks)), area_kind="unavailable",
                    profile_kind="unavailable", shared=NA_integer_, native=NULL) {
  stopifnot(length(spectra)==nrow(peaks), length(profiles)==nrow(peaks))
  list(peaks=peaks, spectra=spectra, profiles=profiles, area_kind=area_kind,
    profile_kind=profile_kind, shared=shared, native=native)
}

positive_pairs <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & y > 0
  cbind(x[keep],y[keep])
}

numeric_pairs <- function(text, pattern) {
  if (length(text)!=1L || is.na(text) || !nzchar(text)) return(matrix(numeric(),0L,2L))
  pattern <- paste0("([+-]?[0-9]+(?:\\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)",pattern,
    "([+-]?[0-9]+(?:\\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)")
  rows <- regmatches(text,gregexpr(pattern,text,perl=TRUE))[[1L]]
  if (!length(rows)) return(matrix(numeric(),0L,2L))
  cbind(as.numeric(sub(pattern,"\\1",rows,perl=TRUE)),
        as.numeric(sub(pattern,"\\2",rows,perl=TRUE)))
}

need_package <- function(name) {
  if (!requireNamespace(name,quietly=TRUE))
    stop("Install ",name," to run this adapter.",call.=FALSE)
}

settings <- function(default, overrides) {
  if (is.null(overrides)) default else utils::modifyList(default,overrides)
}

executable <- function(path, name) {
  if (is.null(path)) path <- unname(Sys.which(name))
  if (length(path)!=1L || !nzchar(path) || !file.exists(path))
    stop("Set the executable path for ",name," in control.",call.=FALSE)
  normalizePath(path,winslash="/",mustWork=TRUE)
}

run_logged <- function(exe, args, work, log, timeout, env=NULL, expected=NULL) {
  if (!is.null(env)) env <- c("current",env)
  proc <- processx::run(exe,args,wd=work,env=env,stdout=log,stderr_to_stdout=TRUE,
    timeout=timeout,error_on_status=FALSE,echo=FALSE,windows_hide_window=TRUE)
  if (isTRUE(proc$timeout)) stop("Timed out. See ",log,call.=FALSE)
  if (is.null(expected) && proc$status!=0L) stop("Exit ",proc$status,". See ",log,call.=FALSE)
  if (!is.null(expected) && !isTRUE(expected())) stop("Expected exports missing. See ",log,call.=FALSE)
  invisible(proc$status)
}

cdf_clock <- function(path) {
  need_package("ncdf4")
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc))
  time <- as.numeric(ncdf4::ncvar_get(nc,"scan_acquisition_time"))
  units <- tolower(ncdf4::ncatt_get(nc,"scan_acquisition_time","units")$value)
  if (!identical(units,"minutes") && !identical(units,"minute")) time <- time/60
  time
}

write_open_ms <- function(cdf, output, format) {
  need_package("mzR")
  input <- mzR::openMSfile(cdf)
  on.exit(mzR::close(input))
  h <- mzR::header(input)
  p <- lapply(mzR::spectra(input),function(x) x[x[,2L] > 0,,drop=FALSE])
  if (any(lengths(p)==0L)) {
    p <- lapply(p, function(x) if(nrow(x)) x else matrix(c(33,0),1L))
  }
  h$peaksCount <- vapply(p,nrow,integer(1))
  h$basePeakMZ <- vapply(p,function(x) x[which.max(x[,2L]),1L],0)
  h$basePeakIntensity <- vapply(p,function(x) max(x[,2L]),0)
  h$lowMZ <- vapply(p,function(x) min(x[,1L]),0)
  h$highMZ <- vapply(p,function(x) max(x[,1L]),0)
  h$centroided <- TRUE
  h$polarity <- 1L
  mzR::writeMSData(p,output,header=h,outformat=tolower(format))
  invisible(output)
}

stage_inputs <- function(files, directory, format="cdf") {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  paths <- file.path(directory,paste0(sample_id(files),".",format))
  for (i in seq_along(paths)) {
    if (file.exists(paths[i])) next
    if (tolower(format)=="cdf") {
      if (!file.copy(files[i],paths[i])) stop("Could not stage ",files[i],call.=FALSE)
    } else write_open_ms(files[i],paths[i],format)
  }
  paths
}
