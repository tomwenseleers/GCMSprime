# Source this script with GCMSbenchmark installed. Time is in minutes.
simulateBenchmark16 <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  sim <- logconNMF::simulateGCMSSamples(samples = 16L, n = 1200L, k = 48L,
    preset = "time_area", presence = 0.5, irregular = 0, drift = 2,
    family = "quasipoisson", dispersion = 2, seed = 9118L)
  ids <- vapply(sim$samples, `[[`, "", "sample_id")
  files <- file.path(directory, paste0(ids, ".cdf"))
  # MSHub reports RTs on its mean alignment clock. Match simulated apices
  # on that clock for this engine. The fits never receive the truth.
  reference <- Reduce(`+`, lapply(sim$samples, `[[`, "time")) / length(ids)
  truth <- lapply(seq_along(ids), function(i) {
    x <- sim$samples[[i]]
    GCMSprime::writeGCMSCDF(x$Y, files[i], time = x$time,
      acquired_at = x$acquired_at)
    present <- which(x$present)
    H <- Matrix::Matrix(x$H_true[, present, drop = FALSE] * x$exposure, sparse = TRUE)
    apex <- max.col(t(as.matrix(H)), ties.method = "first")
    list(id = present, time = x$time, H = H, P = x$P_true[present, , drop = FALSE],
      area = x$area_true[present], apex = x$time[apex], reference_apex = reference[apex],
      acid = rep(FALSE, length(present)), mass = as.integer(colnames(x$Y)))
  })
  list(files = files, truth = stats::setNames(truth, ids))
}

benchmark16Hits <- function(fit, truth, rt_tolerance) {
  norm_rows <- utils::getFromNamespace("norm_rows", "GCMSbenchmark")
  spectrum_matrix <- utils::getFromNamespace("spectrum_matrix", "GCMSbenchmark")
  do.call(rbind, lapply(names(fit$methods), function(method) {
    result <- fit$methods[[method]]
    if (!identical(result$status, "ok")) return(NULL)
    do.call(rbind, lapply(names(truth), function(sample) {
      z <- truth[[sample]]
      i <- which(result$peaks$sample == sample)
      spectra <- result$spectra[i]
      present <- vapply(spectra, function(s) length(s) > 0L &&
        any(is.finite(s[, 2L]) & s[, 2L] > 0), TRUE)
      spectra <- spectra[present]
      rt <- result$peaks$rt[i[present]]
      target <- if (identical(result$rt_basis, "aligned")) z$reference_apex else z$apex
      mass <- sort(unique(round(c(z$mass, unlist(lapply(spectra, function(s) s[, 1L]))))))
      known <- matrix(0, nrow(z$P), length(mass))
      known[, match(round(z$mass), mass)] <- z$P
      cosine <- pmin(pmax(tcrossprod(norm_rows(known),
        norm_rows(spectrum_matrix(spectra, mass))), 0), 1)
      best <- function(x) if (length(x)) max(x) else NA_real_
      nearest <- lapply(target, function(t) which(abs(rt - t) <= rt_tolerance))
      data.frame(method = method, sample = sample, truth_id = z$id,
        best_cosine = vapply(seq_along(target), function(j) best(cosine[j, ]), 0),
        near_cosine = vapply(seq_along(target), function(j) best(cosine[j, nearest[[j]]]), 0))
    }))
  }))
}

runBenchmark16 <- function(output_dir, method = c("GCMSprime", "xcms+CAMERA"),
                           control = list(), ncores = 16L) {
  methods <- utils::getFromNamespace("bench_methods", "GCMSbenchmark")
  if (identical(tolower(method), "all")) method <- methods
  method <- unique(methods[match(tolower(method), tolower(methods))])
  if (!length(method) || anyNA(method)) stop("Choose one or more GCMSbenchmark methods.")
  if (file.exists(output_dir)) stop("Choose a fresh output directory.")
  dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  input <- simulateBenchmark16(file.path(output_dir, "input"))
  saveRDS(input$truth, file.path(output_dir, "truth.rds"))

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
      min_cosine = 0.70, family = stats::quasipoisson("identity"), dispersion = 2,
      precluster = 32L, nwindows = 2L, baseline_df = 4L, correlation = "none",
      psf_control = list(maxiter = 3L),
      batch = list(coupling = list(cycles = 2L, initial_cycles = 1L, rounds = 1L, q = 32L)))
  })
  fit <- structure(list(methods = do.call(c, lapply(fits, `[[`, "methods")),
    files = input$files, directory = output_dir,
    summary = do.call(rbind, lapply(fits, `[[`, "summary")),
    matches = do.call(rbind, lapply(fits, `[[`, "matches"))), class = "GCMSbenchmark")
  hits <- benchmark16Hits(fit, input$truth, rt_tolerance)
  if (is.null(hits)) {
    hits <- data.frame(method = character(), sample = character(), truth_id = integer(),
      best_cosine = numeric(), near_cosine = numeric())
  }
  truth_peaks <- do.call(rbind, lapply(names(input$truth), function(sample) {
    z <- input$truth[[sample]]
    data.frame(sample = sample, truth_id = z$id, rt_min = z$apex, area_true = z$area)
  }))
  utils::write.csv(fit$summary, file.path(output_dir, "benchmark.csv"), row.names = FALSE)
  utils::write.csv(fit$matches, file.path(output_dir, "matches.csv"), row.names = FALSE)
  utils::write.csv(hits, file.path(output_dir, "spectral-hits.csv"), row.names = FALSE)
  utils::write.csv(truth_peaks, file.path(output_dir, "truth-peaks.csv"), row.names = FALSE)
  fit
}
