# Generate the common CDF inputs and a separate scoring key. Times are in minutes.
benchmarkSignalHash <- function(file, mass, scans) {
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc))
  n <- as.integer(ncdf4::ncvar_get(nc, "point_count"))
  time <- ncdf4::ncvar_get(nc, "scan_acquisition_time")
  first <- ncdf4::ncvar_get(nc, "mass_values", start = 1, count = length(mass))
  stopifnot(length(n) == scans, all(n == length(mass)), identical(as.integer(first), as.integer(mass)),
            length(time) == scans, all(diff(time) > 0))
  # The CDF header includes its creation date. Hash the numerical inputs
  # independently, with a specified byte order, for reproducible checks.
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  con <- file(tmp, "wb")
  tryCatch({
    writeBin(as.double(time), con, size = 8L, endian = "little")
    for (name in c("mass_values", "intensity_values")) {
      values <- as.double(ncdf4::ncvar_get(nc, name))
      stopifnot(length(values) == sum(n), all(is.finite(values)))
      writeBin(values, con, size = 8L, endian = "little")
    }
  }, finally = close(con))
  unname(tools::md5sum(tmp))
}

prepareFullBenchmark <- function(directory, background_file = system.file(
    "benchmarks", "gcms-data2-background.rds", package = "GCMSbenchmark", mustWork = TRUE)) {
  if (file.exists(file.path(directory, "inputs.csv"))) stop("This dataset already exists.")
  env <- new.env()
  utils::data("simulatedSpectra", package = "L0glm", envir = env)
  library_spectra <- env$simulatedSpectra
  mass <- 30:579
  ions <- match(as.integer(colnames(library_spectra)), mass)
  spectra <- matrix(0, nrow(library_spectra), length(mass),
                    dimnames = list(rownames(library_spectra), as.character(mass)))
  spectra[, ions] <- as.matrix(library_spectra)
  template <- readRDS(background_file)
  observed_background <- template$basis %*% template$alpha
  background <- matrix(.2 * median(observed_background), nrow(observed_background), length(mass))
  background[, ions] <- observed_background
  n <- 10000L
  k <- 1000L
  samples <- 16L
  dir.create(file.path(directory, "input"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(directory, "standard"), recursive = TRUE, showWarnings = FALSE)
  rng <- RNGkind()
  on.exit(do.call(RNGkind, as.list(rng)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(171116L)
  absent <- t(replicate(k, !(seq_len(samples) %in% sample.int(samples, samples / 2))))
  sim <- logconNMF::simulateGCMSSamples(samples = samples, n = n, k = k, S = spectra,
    preset = "time_area_acid", acid_fraction = .01, absent = absent, presence = .5,
    irregular = .02, drift = .5, rt_area = 0, width_area = .45, tail_area = 1.5,
    family = "quasipoisson", dispersion = template$dispersion, area = template$area,
    background = background, seed = 171117L)
  ids <- vapply(sim$samples, `[[`, "", "sample_id")
  acid <- sim$samples[[1L]]$acid
  stopifnot(identical(dim(sim$present), c(k, samples)), all(rowSums(sim$present) == 8),
            length(acid) == 10L, sum(sim$present[acid, ]) == 80)
  # A specified synthetic RI clock allows every engine's reported times to be
  # scored on the same scale. Pure-compound library RIs are not assumed.
  clock <- function(minutes) 650 + 34 * minutes
  nominal_time <- sim$samples[[1L]]$centers * .3
  ref <- list(id = seq_len(k), P = sim$P, mass = mass, ri = clock(nominal_time))
  calibration <- lapply(sim$samples, function(x) {
    at <- seq(1L, n, length.out = 101L)
    data.frame(ri = clock(x$reference_time[round(at)] * .3), time = x$time[round(at)] * .3)
  })
  reference_time <- Reduce(`+`, lapply(sim$samples, function(x) x$time * .3)) / samples
  mean_calibration <- Reduce(`+`, calibration) / samples
  settings <- list(seed = sim$seed, presence_seed = 171116L, mass = mass, compounds = k,
    scans = n, sim_settings = sim$settings, source_ids = sim$compound_ids, present = sim$present,
    areas = sim$areas, acid = acid, dispersion = template$dispersion,
    representative = c(Regular = setdiff(seq_len(k), acid)[which.max(rowSums(sim$areas)[-acid])],
                       Acid = acid[which.max(rowSums(sim$areas)[acid])]))
  truth <- stats::setNames(vector("list", samples), ids)
  files <- file.path(directory, "input", paste0(ids, ".cdf"))
  for (i in seq_along(ids)) {
    x <- sim$samples[[i]]
    stopifnot(identical(dim(x$Y), c(n, length(mass))))
    time <- x$time * .3
    GCMSprime::writeGCMSCDF(x$Y, files[i], time = time, acquired_at = x$acquired_at)
    present <- which(x$present)
    h <- x$H_true[, present, drop = FALSE] * x$exposure
    apex <- max.col(t(h), ties.method = "first")
    h[h < rep(h[cbind(apex, seq_along(present))] * 1e-12, each = n)] <- 0
    h <- Matrix::Matrix(h, sparse = TRUE)
    h <- h %*% Matrix::Diagonal(x = 1 / Matrix::colSums(h))
    truth[[i]] <- list(id = present, time = time, H = h, P = x$P_true[present, , drop = FALSE],
      area = x$area_true[present], apex = time[apex], reference_apex = reference_time[apex],
      acid = present %in% acid, mass = mass, ri = ref$ri[present], ri_calibration = calibration[[i]],
      reference_ri_calibration = mean_calibration)
    sim$samples[i] <- list(NULL)
    message("Wrote ", ids[i], " with ", length(present), " occurrences.")
  }
  saveRDS(settings, file.path(directory, "settings.rds"))
  saveRDS(truth, file.path(directory, "truth.rds"))
  saveRDS(ref, file.path(directory, "reference-library.rds"))
  rm(sim, truth, x, h)
  gc()
  standard <- logconNMF::simulateGCMSSamples(samples = 1L, n = n, k = 20L, S = spectra,
    preset = "time", acid_fraction = 0, presence = 1, irregular = .02, drift = 0, rt_area = 0,
    family = "quasipoisson", dispersion = template$dispersion, area = template$area,
    background = background, seed = 171118L)$samples[[1L]]
  GCMSprime::writeGCMSCDF(standard$Y, file.path(directory, "standard", "psf-standard.cdf"),
    time = standard$time * .3, acquired_at = as.POSIXct("2026-11-01 07:00:00", tz = "UTC"))
  manifest <- data.frame(sample = ids, file = basename(files), md5 = unname(tools::md5sum(files)),
    scans = n, channels = length(mass), occurrences = colSums(settings$present))
  manifest$signal_md5 <- unname(vapply(files, benchmarkSignalHash, "", mass = mass, scans = n))
  write.csv(manifest, file.path(directory, "inputs.csv"), row.names = FALSE)
  invisible(manifest)
}

# Run once per program. No scoring information enters this function.
fitFullBenchmark <- function(directory, method, control = list(), ncores = 16L) {
  methods <- c(GCMSprime = "gcmsprime", `MS-DIAL` = "msdial", mzmine = "mzmine",
    eRah = "erah", AMDIS = "amdis", `xcms+CAMERA` = "xcms_camera", MSHub = "mshub")
  if (is.null(control)) control <- list()
  stopifnot(length(method) == 1L, method %in% names(methods), is.list(control),
            length(ncores) == 1L, is.finite(ncores), ncores >= 1L)
  directory <- normalizePath(directory, winslash = "/")
  files <- sort(list.files(file.path(directory, "input"), "[.]cdf$", full.names = TRUE))
  stopifnot(length(files) == 16L, !anyDuplicated(basename(files)))
  out <- file.path(directory, "runs", method)
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
    stop("A run directory already exists for ", method, ". Choose a fresh location for a new fit.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  raw_file <- file.path(out, "raw.rds")
  start <- proc.time()[["elapsed"]]
  set.seed(171117L)
  if (method == "GCMSprime") {
    standard_path <- file.path(directory, "standard", "psf-standard.cdf")
    raw <- utils::getFromNamespace("GCMS_read_cdf", "GCMSprime")(standard_path, NULL)
    standard <- logconNMF::estimatePSF(raw$Y, time = raw$time, family = quasipoisson("identity"),
      method = "standard", nwindows = 10L, precluster = 1L, ncores = min(8L, ncores))
    rm(raw)
    calibration_s <- proc.time()[["elapsed"]] - start
    fit_start <- proc.time()[["elapsed"]]
    fit <- GCMSprime::gcmsDeco(dirname(files[1L]), output_dir = file.path(out, "analysis"),
      ncores = ncores, psf_control = list(standard = standard),
      batch = list(workers = min(8L, ncores), coupling = list(penalty_search = TRUE)), verbose = TRUE)
    elapsed <- proc.time()[["elapsed"]] - fit_start
    saveRDS(list(seconds = elapsed, calibration_s = calibration_s, result = fit),
            file.path(out, "native.rds"), compress = FALSE)
    result <- utils::getFromNamespace("gcmsprime_batch", "GCMSbenchmark")(fit,
      tools::file_path_sans_ext(basename(files)))
    result$analysis <- elapsed + calibration_s
    result$calibration <- calibration_s
  } else {
    defaults <- switch(method,
      "MS-DIAL" = list(method_file = system.file("benchmarks", "msdial-simulated-16.txt",
         package = "GCMSbenchmark", mustWork = TRUE)),
      eRah = list(deconvolution = list(min.peak.width = 2, min.peak.height = 20, noise.threshold = 3),
                  alignment = list(mz.range = 30:579)),
      "xcms+CAMERA" = list(detection = list(fwhm = 2.6, snthresh = 1.5),
                  fwhm_grouping = list(perfwhm = .8), grouping = FALSE, min_ions = 4L), list())
    engine <- utils::getFromNamespace(paste0("bench_", methods[[method]]), "GCMSbenchmark")
    cores <- switch(method, eRah = min(8L, ncores), "xcms+CAMERA" = min(4L, ncores), ncores)
    result <- engine(files, out, cores, utils::modifyList(defaults, control), 14400, list())
  }
  result$total <- proc.time()[["elapsed"]] - start
  result$status <- "ok"
  result$message <- ""
  saveRDS(result, raw_file)
  write.csv(data.frame(method = method, seconds_per_sample = result$analysis / length(files),
    occurrences_called = nrow(result$peaks)), file.path(out, "timing.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
  invisible(result)
}

# Link sample-local matches to the engine's original shared component identifiers.
benchmarkFactorMatches <- function(record, matches) {
  peak <- record$peaks[matches$peak_row, , drop = FALSE]
  stopifnot(all(peak$sample == matches$sample), all(peak$id == matches$estimate_id))
  used <- !is.na(peak$group) & nzchar(peak$group)
  data.frame(method = matches$method[used], sample = matches$sample[used],
    estimate_id = matches$estimate_id[used], factor_id = peak$group[used], truth_id = matches$truth_id[used])
}

# Score saved fits independently of their timing. Curves retain total-ion counts per scan.
scoreFullBenchmark <- function(directory, method, output_dir = file.path(directory, "scores")) {
  input <- file.path(directory, "input")
  files <- sort(list.files(input, "[.]cdf$", full.names = TRUE))
  truth <- readRDS(file.path(directory, "truth.rds"))
  reference <- readRDS(file.path(directory, "reference-library.rds"))
  settings <- readRDS(file.path(directory, "settings.rds"))
  run <- file.path(directory, "runs", method)
  record <- readRDS(file.path(run, "raw.rds"))
  record <- utils::getFromNamespace("tic_units", "GCMSbenchmark")(record, method, files)
  fit <- structure(list(methods = setNames(list(record), method), files = files),
                   class = "GCMSbenchmark")
  fit <- GCMSbenchmark::scoreGCMSbenchmark(fit, truth, reference = reference,
    matching = "reciprocal", ri_tolerance = 30, min_cosine = .70, ion_threshold = .001)
  z <- fit$summary
  stopifnot(z$TP + z$FN == sum(settings$present), z$TP + z$FP == nrow(record$peaks),
    z$true_rank == nrow(reference$P), z$acid_truth == sum(settings$present[settings$acid, ]),
    !anyDuplicated(fit$matches[c("sample", "truth_id")]), !anyDuplicated(fit$matches[c("sample", "estimate_id")]))
  count <- length(files)
  z$analysis_s_per_sample <- z$analysis_s / count
  z$conversion_s_per_sample <- z$input_conversion_s / count
  z$analysis_and_conversion_s_per_sample <- z$analysis_s_per_sample +
    ifelse(is.na(z$conversion_s_per_sample), 0, z$conversion_s_per_sample)
  z$adapter_s_per_sample <- z$adapter_elapsed_s / count
  z$calibration_s_per_sample <- if (length(record$calibration)) record$calibration / count else NA_real_
  z[c("analysis_s", "input_conversion_s", "adapter_elapsed_s")] <- NULL
  out <- file.path(output_dir, method)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  put <- function(d, n) utils::write.csv(d, file.path(out, paste0(n, ".csv")), row.names = FALSE)
  put(z, "benchmark")
  put(fit$matches, "matches")
  put(fit$identifications, "identifications")
  if (is.finite(z$estimated_rank)) put(benchmarkFactorMatches(record, fit$matches), "factor-matches")
  peaks <- do.call(rbind, lapply(names(truth), function(id) {
    x <- truth[[id]]
    data.frame(sample = id, truth_id = x$id, rt_min = x$apex, area_true = x$area, ri = x$ri, acid = x$acid)
  }))
  keys <- function(sample, id) paste(sample, id, sep = "/")
  matched <- keys(peaks$sample, peaks$truth_id) %in% keys(fit$matches$sample, fit$matches$truth_id)
  quartile <- cut(peaks$area_true, quantile(peaks$area_true, seq(0, 1, .25)), include.lowest = TRUE,
                  labels = c("Smallest", "Lower-middle", "Upper-middle", "Largest"))
  tab <- do.call(rbind, lapply(levels(quartile), function(q) {
    at <- which(quartile == q)
    data.frame(method, quartile = q, true = length(at), TP = sum(matched[at]),
               sensitivity = mean(matched[at]), min_area = min(peaks$area_true[at]), max_area = max(peaks$area_true[at]))
  }))
  put(tab, "sensitivity-by-area")
  per_sample <- do.call(rbind, lapply(names(truth), function(s) {
    tp <- sum(fit$matches$sample == s)
    data.frame(method, sample = s, true = length(truth[[s]]$id), TP = tp,
               FP = sum(record$peaks$sample == s) - tp, FN = length(truth[[s]]$id) - tp)
  }))
  put(per_sample, "by-sample")
  if (method == "GCMSprime") {
    put(peaks, "truth-peaks")
    # A relative time grid for each compound is available separately for every acid.
    put(data.frame(kind = names(settings$representative), truth_id = unname(settings$representative),
      total_area = rowSums(settings$areas)[settings$representative],
      present = rowSums(settings$present)[settings$representative]), "profile-targets")
  }
  acid_quality <- do.call(rbind, lapply(settings$acid, function(id) {
    q <- subset(fit$matches, truth_id == id)
    finite <- is.finite(q$profile_cosine)
    data.frame(method, truth_id = id, present = sum(settings$present[id, ]),
      recovered = nrow(q), mean_matched_profile_cosine = if(any(finite)) mean(q$profile_cosine[finite]) else NA_real_,
      median_area_ratio = if(nrow(q)) median(q$area_estimated/q$area_true, na.rm = TRUE) else NA_real_)
  }))
  put(acid_quality, "acid-by-compound")
  profiles <- list()
  quality <- list()
  if (any(lengths(record$profiles) > 0)) for (kind in names(settings$representative)) {
    id <- unname(settings$representative[kind])
    abundance <- settings$areas[id, ]
    ranks <- rank(abundance, ties.method = "first")
    q <- fit$matches[fit$matches$truth_id == id, , drop = FALSE]
    quality[[kind]] <- data.frame(kind, method, present = sum(abundance > 0), recovered = nrow(q),
      mean_profile_cosine = sum(q$profile_cosine, na.rm = TRUE) / sum(abundance > 0))
    for (i in seq_along(truth)) {
      s <- names(truth)[i]
      x <- truth[[i]]
      center <- reference$ri[id]
      ri <- function(t) stats::approx(x$ri_calibration$time, x$ri_calibration$ri, t, rule = 2)$y
      add <- function(engine, comp, t, h, reciprocal, called) {
        coordinate <- ri(t)
        keep <- which(abs(coordinate - center) <= 30)
        if (!length(keep)) return()
        profiles[[length(profiles) + 1L]] <<- data.frame(kind, truth_id = id, method = engine,
          sample = s, sample_rank = ranks[i], true_area = abundance[i], component = as.character(comp),
          ri = coordinate[keep], intensity = h[keep], reciprocal, called)
      }
      if (method == "GCMSprime") {
        j <- match(id, x$id)
        h <- if (is.na(j)) numeric(length(x$time)) else as.numeric(x$H[, j]) * abundance[i]
        add("Ground truth", id, x$time, h, TRUE, !is.na(j))
      }
      calls <- subset(fit$identifications, sample == s & reference_id %in% id)
      have <- FALSE
      for (j in seq_len(nrow(calls))) {
        p <- record$profiles[[calls$peak_row[j]]]
        if (length(p)) {
          add(method, calls$estimate_id[j], p[, 1], p[, 2], calls$reciprocal[j], TRUE)
          have <- TRUE
        }
      }
      if (!have) add(method, "none", x$time, numeric(length(x$time)), FALSE, FALSE)
    }
  }
  if (length(quality)) put(do.call(rbind, quality), "profile-quality")
  saveRDS(do.call(rbind, profiles), file.path(out, "profiles.rds"))
  invisible(z)
}

collectFullBenchmark <- function(directory, methods, output_dir) {
  scores <- file.path(directory, "scores")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (key in c("benchmark", "matches", "identifications", "factor-matches", "by-sample", "sensitivity-by-area",
                "truth-peaks", "profile-targets", "profile-quality", "acid-by-compound")) {
    paths <- file.path(scores, methods, paste0(key, ".csv"))
    if (key == "benchmark") stopifnot(all(file.exists(paths)))
    d <- do.call(rbind, lapply(paths[file.exists(paths)], utils::read.csv, check.names = FALSE))
    utils::write.csv(d, file.path(output_dir, paste0(key, ".csv")), row.names = FALSE)
    if (key == "identifications")
      utils::write.csv(d[d$reciprocal, ], file.path(output_dir, "spectral-hits.csv"), row.names = FALSE)
  }
  curves <- do.call(rbind, lapply(file.path(scores, methods, "profiles.rds"), readRDS))
  saveRDS(curves, file.path(output_dir, "profiles.rds"))
  manifest <- c("inputs.csv", "settings.rds", "reference-library.rds")
  stopifnot(all(file.copy(file.path(directory, manifest), output_dir, overwrite = TRUE)))
}
