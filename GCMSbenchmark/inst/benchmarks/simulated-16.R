# Source this script with GCMSbenchmark installed. Time is in minutes.
benchmark16Simulation <- function() {
  library <- new.env()
  utils::data("simulatedSpectra",package="L0glm",envir=library)
  spectra <- library$simulatedSpectra
  ri <- unname(attr(spectra,"RI"))
  recorded_time <- unname(attr(spectra,"RT"))
  # The library attributes follow matrix row order. Their names retain source IDs.
  scaled_time <- .6 + (recorded_time-min(recorded_time)) /
    diff(range(recorded_time)) * 18.8
  calibration <- data.frame(ri=ri,time=scaled_time)
  logconNMF::simulateGCMSSamples(samples=16L,n=1200L,k=48L,S=spectra,
    preset="time_area_acid",acid_fraction=1/48,presence=.5,irregular=0,
    drift=.5,rt_area=0,family="quasipoisson",dispersion=2,seed=9118L,
    retentionIndex=ri,retentionMap=calibration)
}

simulateBenchmark16 <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  sim <- benchmark16Simulation()
  ids <- vapply(sim$samples, `[[`, "", "sample_id")
  files <- file.path(directory, paste0(ids, ".cdf"))
  # MSHub reports RTs on its mean alignment clock. Match simulated apices
  # on that clock for this engine. The fits never receive the truth.
  reference <- Reduce(`+`, lapply(sim$samples, `[[`, "time")) / length(ids)
  mean_calibration <- Reduce(`+`,lapply(sim$samples,function(x)
    x$retentionCalibration$time))/length(ids)
  truth <- lapply(seq_along(ids), function(i) {
    x <- sim$samples[[i]]
    GCMSprime::writeGCMSCDF(x$Y, files[i], time = x$time,
      acquired_at = x$acquired_at)
    present <- which(x$present)
    H <- Matrix::Matrix(x$H_true[, present, drop = FALSE] * x$exposure, sparse = TRUE)
    apex <- max.col(t(as.matrix(H)), ties.method = "first")
    list(id = present, time = x$time, H = H, P = x$P_true[present, , drop = FALSE],
      area = x$area_true[present], apex = x$time[apex], reference_apex = reference[apex],
      acid = present%in%x$acid, mass = as.integer(colnames(x$Y)),
      ri=x$retentionIndex[present],ri_calibration=x$retentionCalibration,
      reference_ri_calibration=data.frame(ri=sim$retentionMap$ri,time=mean_calibration))
  })
  known <- sim$samples[[1L]]
  reference_library <- list(id=seq_len(nrow(sim$P)),P=sim$P,
    mass=as.integer(colnames(sim$P)),ri=known$retentionIndex)
  selected <- which(vapply(sim$samples,function(x)any(x$present[x$acid]),TRUE))[1L]
  x <- sim$samples[[selected]]
  saveRDS(list(sample=x$sample_id,time=x$time,mass=as.numeric(colnames(x$Y)),
    observed=x$Y,mean=x$mean,apex=x$time[max.col(t(x$H_true),ties.method="first")],
    present=x$present,acid=x$acid),file.path(dirname(directory),"sample-heatmap.rds"))
  list(files = files, truth = stats::setNames(truth, ids),reference=reference_library)
}

runBenchmark16 <- function(output_dir, method = c("GCMSprime", "xcms+CAMERA"),
                           control = list(), ncores = 16L) {
  methods <- utils::getFromNamespace("bench_methods", "GCMSbenchmark")
  if (identical(tolower(method), "all")) {
    method <- methods
    if (.Platform$OS.type!="windows" && is.null(control$AMDIS$elu))
      method <- setdiff(method,"AMDIS")
  }
  method <- unique(methods[match(tolower(method), tolower(methods))])
  if (!length(method) || anyNA(method)) stop("Choose one or more GCMSbenchmark methods.")
  if (file.exists(output_dir)) stop("Choose a fresh output directory.")
  dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  input <- simulateBenchmark16(file.path(output_dir, "input"))
  saveRDS(input$truth, file.path(output_dir, "truth.rds"))
  saveRDS(input$reference,file.path(output_dir,"reference-library.rds"))

  settings <- list(
    "MS-DIAL" = list(method_file = system.file("benchmarks", "msdial-simulated-16.txt",
      package = "GCMSbenchmark", mustWork = TRUE)),
    eRah = list(deconvolution = list(min.peak.width = 4, min.peak.height = 20,
      noise.threshold = 3), alignment = list(mz.range = 30:550)),
    "xcms+CAMERA" = list(detection = list(fwhm = 8.6, snthresh = 3),
      fwhm_grouping = list(perfwhm = 0.8), grouping = FALSE, min_ions = 6L))
  control <- utils::modifyList(settings, control)
  rt_tolerance <- 0.12
  fits <- lapply(method, function(engine) {
    GCMSbenchmark::GCMSbenchmark(input$files,
      output_dir = file.path(output_dir, "runs", engine), method = engine,
      truth = input$truth, ncores = if (engine == "xcms+CAMERA") min(4L, ncores) else ncores,
      control = control[intersect(names(control), engine)], rt_tolerance = rt_tolerance,
      min_cosine = 0.70,matching="reciprocal",ri_tolerance=30,reference=input$reference,
      family = stats::quasipoisson("identity"), dispersion = 2,
      precluster = 32L, nwindows = 2L, baseline_df = 4L, correlation = "none",
      psf_control = list(maxiter = 3L),
      batch = list(coupling = list(cycles = 2L, initial_cycles = 1L, rounds = 1L, q = 32L)))
  })
  fit <- structure(list(methods = do.call(c, lapply(fits, `[[`, "methods")),
    files = input$files, directory = output_dir,
    reference=input$reference,matching="reciprocal",
    summary = do.call(rbind, lapply(fits, `[[`, "summary")),
    matches = do.call(rbind, lapply(fits, `[[`, "matches")),
    identifications=do.call(rbind,lapply(fits,`[[`,"identifications"))), class = "GCMSbenchmark")
  hits <- fit$identifications[fit$identifications$reciprocal,,drop=FALSE]
  truth_peaks <- do.call(rbind, lapply(names(input$truth), function(sample) {
    z <- input$truth[[sample]]
    data.frame(sample = sample, truth_id = z$id,rt_min = z$apex,area_true = z$area,ri=z$ri,acid=z$acid)
  }))
  utils::write.csv(fit$summary, file.path(output_dir, "benchmark.csv"), row.names = FALSE)
  utils::write.csv(fit$matches, file.path(output_dir, "matches.csv"), row.names = FALSE)
  utils::write.csv(hits, file.path(output_dir, "spectral-hits.csv"), row.names = FALSE)
  utils::write.csv(fit$identifications,file.path(output_dir,"identifications.csv"),row.names=FALSE)
  utils::write.csv(truth_peaks, file.path(output_dir, "truth-peaks.csv"), row.names = FALSE)
  saveRDS(fit,file.path(output_dir,"benchmark.rds"))
  benchmark16AcidProfiles(fit,input$truth,output_dir)
  fit
}

benchmark16AcidProfiles <- function(fit,truth,output_dir) {
  acid <- unique(unlist(lapply(truth,function(z)z$id[z$acid])))
  profiles <- list()
  for (sample in names(truth)) {
    z <- truth[[sample]]
    j <- match(acid,z$id)
    if (is.na(j)) next
    add <- function(method,component,time,intensity,reciprocal) {
      keep <- is.finite(intensity)&intensity>max(intensity,na.rm=TRUE)*.001
      if(!any(keep)) return(NULL)
      ri <- utils::getFromNamespace("ri_from_time","GCMSbenchmark")(time[keep],z$ri_calibration)
      data.frame(method=method,sample=sample,component=as.character(component),ri=ri,
        intensity=intensity[keep]/max(intensity),reciprocal=reciprocal)
    }
    profiles[[length(profiles)+1L]] <- add("Ground truth",acid,z$time,z$H[,j],TRUE)
    for (method in names(fit$methods)) {
      calls <- fit$identifications[fit$identifications$method==method &
        fit$identifications$sample==sample & fit$identifications$reference_id%in%as.character(acid),]
      for (k in seq_len(nrow(calls))) {
        p <- fit$methods[[method]]$profiles[[calls$peak_row[k]]]
        if(length(p)) profiles[[length(profiles)+1L]] <-
          add(method,calls$estimate_id[k],p[,1L],p[,2L],calls$reciprocal[k])
      }
    }
  }
  utils::write.csv(do.call(rbind,profiles),file.path(output_dir,"acid-profiles.csv"),row.names=FALSE)
}
