# Source this script with GCMSbenchmark installed. Time is in minutes.
benchmark16Simulation <- function() {
  rng<-RNGkind()
  on.exit(do.call(RNGkind,as.list(rng)))
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  library <- new.env()
  utils::data("simulatedSpectra",package="L0glm",envir=library)
  spectra <- library$simulatedSpectra
  ri <- unname(attr(spectra,"RI"))
  recorded_time <- unname(attr(spectra,"RT"))
  # The library attributes follow matrix row order. Their names retain source IDs.
  scaled_time <- .6 + (recorded_time-min(recorded_time)) /
    diff(range(recorded_time)) * 18.8
  calibration <- data.frame(ri=ri,time=scaled_time)
  template <- readRDS(system.file("benchmarks","gcms-data2-background.rds",
    package="GCMSbenchmark",mustWork=TRUE))
  logconNMF::simulateGCMSSamples(samples=16L,n=1200L,k=100L,S=spectra,
    preset="time_area_acid",acid_fraction=1/100,presence=.5,irregular=0,
    drift=.5,rt_area=0,family="quasipoisson",dispersion=template$dispersion,seed=9118L,
    area=template$area,background=template$basis%*%template$alpha,
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
  list(files = files, truth = stats::setNames(truth, ids),reference=reference_library,sim=sim)
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
    "xcms+CAMERA" = list(detection = list(fwhm = 8.6, snthresh = 1.5),
      fwhm_grouping = list(perfwhm = 0.8), grouping = FALSE, min_ions = 4L))
  control <- utils::modifyList(settings, control)
  rt_tolerance <- 0.12
  fits <- lapply(method, function(engine) {
    GCMSbenchmark::GCMSbenchmark(input$files,
      output_dir = file.path(output_dir, "runs", engine), method = engine,
      truth = input$truth, ncores = if (engine == "xcms+CAMERA") min(4L, ncores) else ncores,
      control = control[intersect(names(control), engine)], rt_tolerance = rt_tolerance,
      min_cosine = 0.70,matching="reciprocal",ri_tolerance=30,reference=input$reference,
      ion_threshold=.001,
      family = stats::quasipoisson("identity"),
      precluster = 32L, nwindows = 2L, baseline_df = 4L, correlation = "none",
      psf_control = list(maxiter = 3L),
      batch = list(coupling = list(cycles = 2L, initial_cycles = 1L, rounds = 1L, q = 32L)))
  })
  fit <- structure(list(methods = do.call(c, lapply(fits, `[[`, "methods")),
    files = input$files, directory = output_dir,
    reference=input$reference,matching="reciprocal",ion_threshold=.001,
    simulation=list(seed=input$sim$seed,compounds=nrow(input$sim$P),
      samples=length(input$sim$samples),occurrences=sum(input$sim$present),
      dispersion=input$sim$settings$dispersion),
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
  benchmark16Profiles(fit,input$truth,input$sim,output_dir)
  fit
}

benchmark16Profiles <- function(fit,truth,sim,output_dir) {
  acid <- sim$samples[[1L]]$acid
  regular <- setdiff(seq_len(nrow(sim$P)),acid)
  regular <- regular[which.max(rowSums(sim$areas)[regular])]
  targets <- c(Regular=regular,Acid=acid)
  available <- names(fit$methods)[vapply(fit$methods,function(f)
    identical(f$status,"ok")&&any(lengths(f$profiles)>0),TRUE)]
  profiles <- list()
  for(kind in names(targets)) {
    id<-targets[[kind]]
    abundance<-sim$areas[id,]
    ranks<-rank(abundance,ties.method="first")
    for(s in seq_along(sim$samples)) {
      x<-sim$samples[[s]];sample<-x$sample_id;z<-truth[[sample]]
      add<-function(method,component,time,intensity,reciprocal,called) {
        ri<-utils::getFromNamespace("ri_from_time","GCMSbenchmark")(time,z$ri_calibration)
        data.frame(kind=kind,truth_id=id,method=method,sample=sample,sample_rank=ranks[s],
          true_area=abundance[s],component=as.character(component),ri=ri,
          intensity=as.numeric(intensity),reciprocal=reciprocal,called=called)
      }
      profiles[[length(profiles)+1L]]<-add("Ground truth",id,x$time,
        x$H_true[,id]*x$exposure*x$area_true[id],TRUE,x$present[id])
      for(method in available) {
        calls<-fit$identifications[fit$identifications$method==method &
          fit$identifications$sample==sample &
          fit$identifications$reference_id%in%as.character(id),,drop=FALSE]
        nonempty<-FALSE
        for(k in seq_len(nrow(calls))) {
          p<-fit$methods[[method]]$profiles[[calls$peak_row[k]]]
          if(length(p)) {
            profiles[[length(profiles)+1L]]<-add(method,calls$estimate_id[k],
              p[,1L],p[,2L],calls$reciprocal[k],TRUE)
            nonempty<-TRUE
          }
        }
        if(!nonempty)profiles[[length(profiles)+1L]]<-add(method,"none",x$time,
          rep(0,length(x$time)),FALSE,FALSE)
      }
    }
  }
  saveRDS(do.call(rbind,profiles),file.path(output_dir,"profiles.rds"))
  target_info<-data.frame(kind=names(targets),truth_id=unname(targets),
    total_area=rowSums(sim$areas)[targets],present=rowSums(sim$present)[targets])
  utils::write.csv(target_info,file.path(output_dir,"profile-targets.csv"),row.names=FALSE)
  quality<-do.call(rbind,lapply(names(targets),function(kind) {
    id<-targets[[kind]];n<-sum(sim$present[id,])
    do.call(rbind,lapply(available,function(method) {
      q<-fit$matches[fit$matches$method==method &fit$matches$truth_id==id,,drop=FALSE]
      data.frame(kind=kind,method=method,present=n,recovered=nrow(q),
        mean_profile_cosine=sum(q$profile_cosine,na.rm=TRUE)/n,
        median_matched_profile_cosine=if(nrow(q))stats::median(q$profile_cosine,na.rm=TRUE)else NA_real_)
    }))
  }))
  utils::write.csv(quality,file.path(output_dir,"profile-quality.csv"),row.names=FALSE)
  areas<-do.call(rbind,lapply(seq_along(sim$samples),function(s) {
    x<-sim$samples[[s]]
    known<-data.frame(sample=x$sample_id,truth_id=seq_len(nrow(sim$P)),
      area_true=x$area_true,present=x$present)
    rows<-lapply(c("Ground truth",names(fit$methods)),function(method) {
      r<-known;r$method<-method;r$area_estimated<-if(method=="Ground truth")r$area_true else 0
      r$called<-if(method=="Ground truth")r$present else FALSE
      if(method!="Ground truth") {
        hits<-fit$identifications[fit$identifications$method==method &
          fit$identifications$sample==x$sample_id &fit$identifications$reciprocal,,drop=FALSE]
        i<-match(hits$reference_id,as.character(r$truth_id))
        r$area_estimated[i]<-fit$methods[[method]]$peaks$area[hits$peak_row]
        r$called[i]<-TRUE
      }
      r$false_negative<-r$present & !r$called
      r
    })
    do.call(rbind,rows)
  }))
  utils::write.csv(areas,file.path(output_dir,"area-grid.csv"),row.names=FALSE)
}
