bench_xcms_camera <- function(files,out,ncores,control,timeout,dots) {
  need_package("xcms")
  need_package("CAMERA")
  min_ions <- if (is.null(control$min_ions)) 3L else control$min_ions
  if (length(min_ions)!=1L || !is.numeric(min_ions) || !is.finite(min_ions) ||
      min_ions<1L || min_ions!=floor(min_ions))
    stop("xcms+CAMERA min_ions must be a positive integer.",call.=FALSE)
  workers <- min(ncores,length(files))
  config <- file.path(out,"settings.rds")
  result <- file.path(out,"xcms-camera.rds")
  saveRDS(list(files=files,libraries=.libPaths(),workers=workers,
    min_ions=min_ions,detection=control$detection,grouping=control$grouping),config)
  rscript <- file.path(R.home("bin"),if (.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
  run <- bench_time(run_logged(rscript,c("--vanilla",adapter_asset("run-xcms-camera.R"),
    config,result),out,file.path(out,"run.log"),timeout))
  native <- readRDS(result)
  parts <- lapply(seq_along(files),function(i) {
    sample <- native[[i]]
    p <- sample$peaks
    records(peak_table(sample_id(files[i]),p$id,p$rt,p$area,p$height),
      sample$spectra,sample$profiles)
  })
  answer <- combine_records(parts,"xcms summed ion area (count-seconds)","xcms representative-ion trace")
  answer$analysis <- run$seconds
  answer$workers <- workers
  if (length(files)>1L) answer$note <- "CAMERA groups coeluting ions separately within each sample."
  answer
}
