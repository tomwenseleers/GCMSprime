bench_methods <- c("GCMSprime","MS-DIAL","mzmine","eRah","AMDIS","xcms+CAMERA","MSHub")

GCMSbenchmark <- function(cdf_file, output_dir, method="GCMSprime", truth=NULL,
                          ncores=1L, control=list(), timeout=7200,
                          rt_tolerance=.12, min_cosine=.70,
                          matching=c("assignment","reciprocal"),ri_tolerance=30,reference=NULL, ...) {
  matching <- match.arg(matching)
  files <- cdf_file
  if (length(cdf_file)==1L && dir.exists(cdf_file)) {
    files <- list.files(cdf_file,pattern="\\.cdf$",recursive=TRUE,ignore.case=TRUE)
    files <- file.path(cdf_file,files[!grepl("(^|/)standards/",files,ignore.case=TRUE)])
  }
  if (!length(files) || any(!file.exists(files)) || any(dir.exists(files)) ||
      anyDuplicated(sample_id(files)))
    stop("Supply CDF files with unique sample names or a directory containing them.",call.=FALSE)
  files <- normalizePath(files,winslash="/",mustWork=TRUE)
  if (matching=="reciprocal")
    invisible(scoreGCMSbenchmark(structure(list(methods=list(),files=files),class="GCMSbenchmark"),
      truth,rt_tolerance,min_cosine,matching,ri_tolerance,reference))
  if (!is.numeric(ncores) || length(ncores)!=1L || !is.finite(ncores) ||
      ncores<1 || ncores!=floor(ncores) || ncores>.Machine$integer.max)
    stop("ncores must be a positive integer.",call.=FALSE)
  if (!is.numeric(timeout) || length(timeout)!=1L || !is.finite(timeout) || timeout<=0)
    stop("timeout must be positive seconds.",call.=FALSE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!names(control)%in%bench_methods))))
    stop("Use named backend-specific lists in control.",call.=FALSE)
  all_methods <- identical(tolower(method),"all")
  if (all_methods) {
    method <- if(length(files)<10L) setdiff(bench_methods,"MSHub") else bench_methods
    if (.Platform$OS.type!="windows" && is.null(control$AMDIS$elu))
      method <- setdiff(method,"AMDIS")
  }
  method <- bench_methods[match(tolower(method),tolower(bench_methods))]
  if (!length(method) || anyNA(method))
    stop("method must be one or more of ",paste(bench_methods,collapse=", "),".",call.=FALSE)
  method <- unique(method)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir <- normalizePath(output_dir,winslash="/",mustWork=TRUE)
  if (file.exists(file.path(output_dir,"benchmark.csv")) ||
      any(file.exists(file.path(output_dir,method))))
    stop("Choose a fresh output directory for this benchmark.",call.=FALSE)
  dots <- list(...)
  backends <- list(GCMSprime=bench_gcmsprime, "MS-DIAL"=bench_msdial, mzmine=bench_mzmine,
    eRah=bench_erah, AMDIS=bench_amdis, "xcms+CAMERA"=bench_xcms_camera, MSHub=bench_mshub)
  fits <- stats::setNames(lapply(method,function(m) {
    directory <- file.path(output_dir,m)
    dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(backends[[m]](files,directory,ncores,control[[m]],timeout,
      if(m=="GCMSprime")dots else list()),error=identity)
    total <- proc.time()[["elapsed"]]-started
    if (inherits(result,"error")) {
      return(list(status="error",message=conditionMessage(result),total=total))
    }
    result$status <- "ok"; result$message <- ""; result$total <- total
    saveRDS(result,file.path(directory,"standardized.rds"))
    result
  }),method)
  answer <- structure(list(methods=fits, files=files, directory=output_dir),
    class="GCMSbenchmark")
  answer <- scoreGCMSbenchmark(answer,truth,rt_tolerance,min_cosine,matching,ri_tolerance,reference)
  utils::write.csv(answer$summary,file.path(output_dir,"benchmark.csv"),row.names=FALSE)
  if (nrow(answer$matches)) utils::write.csv(answer$matches,
    file.path(output_dir,"matches.csv"),row.names=FALSE)
  answer
}

print.GCMSbenchmark <- function(x, ...) {
  keep <- intersect(c("method","status","analysis_s","components","shared","TP","FP",
    "FN","F1","spectral_cosine","profile_cosine"),names(x$summary))
  print(x$summary[,keep,drop=FALSE],row.names=FALSE,...)
  invisible(x)
}
