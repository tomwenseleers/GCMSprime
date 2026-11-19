test_that("CDF initialization exposes iterative solvers and retries workspace errors", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("logconNMF", "0.1.0.9003")
  set.seed(472);n <- 60L;m <- 4L
  kernel <- logconNMF::convKernel("gaussian", 9L, list(sigma = 1))
  X <- logconNMF::convMatrix(kernel, n)
  mu <- as.matrix(X[, c(15L, 42L)] %*%
    rbind(c(8000, 4000, 2000, 1000), c(1000, 3000, 6000, 4000))) + 5
  Y <- matrix(rpois(n * m, mu), n, m)
  input <- tempfile(fileext = ".cdf")
  outputs <- c(tempfile(fileext = ".msp"), tempfile(fileext = ".msp"))
  on.exit(unlink(c(input, outputs)))
  writeGCMSCDF(Y, input, time = seq_len(n) / 60, mass = 40:43)
  controls <- list(initial_solver = list(max.memory = 0))
  args <- list(cdf_file = input, msp_file = outputs[1L],
    family = quasipoisson("identity"), dispersion = 2, precluster = 2L,
    baseline_df = 4L, nwindows = 1L, correlation = "none", ncores = 1L,
    psf_control = list(kernel_est = kernel, maxiter = 0L), control = controls)
  expect_message(fallback <- do.call(gcmsDeco, args),
    "Eigen CG and diagonal preconditioning")
  args$msp_file <- outputs[2L]
  args$control$initial_solver <- list(solver = "eigen", method = 3L,
    precond = 1L, max.memory = 0)
  args$control$initial_fallback <- FALSE
  explicit <- do.call(gcmsDeco, args)
  expect_s3_class(fallback, "GCMSdeco")
  expect_gt(fallback$model$rank, 0L)
  expect_equal(coef(fallback), coef(explicit), tolerance = 1e-8)
  expect_equal(fitted(fallback), fitted(explicit), tolerance = 1e-8)
  expect_true(all(file.exists(outputs)))
  expect_identical(readLines(outputs[1L]), readLines(outputs[2L]))
  expect_equal(dim(fitted(fallback)), dim(Y))
  expect_true(all(is.finite(fitted(fallback))))
})

test_that("coupled and independent sample series finish with a zero Gram budget", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("logconNMF", "0.1.0.9003")
  set.seed(479);n <- 75L;m <- 5L
  kernel <- logconNMF::convKernel("gaussian", 9L, list(sigma = 1))
  X <- logconNMF::convMatrix(kernel, n)
  S <- rbind(c(8000, 4000, 2000, 1000, 3000), c(1000, 3000, 6000, 4000, 2000))
  directory <- tempfile();dir.create(directory)
  outputs <- replicate(3L, tempfile())
  on.exit(unlink(c(directory, outputs), recursive = TRUE))
  for (sample in 1:2) {
    mu <- as.matrix(X[, c(20L, 52L) + sample - 1L] %*% S) + 5
    writeGCMSCDF(matrix(rpois(n * m, mu), n, m),
      file.path(directory, paste0("sample", sample, ".cdf")),
      time = seq_len(n) / 60, mass = 40:44,
      acquired_at = as.POSIXct("2026-11-01", tz = "UTC") + sample * 3600)
  }
  controls <- list(initial_solver = list(method = 3L, precond = 1L, max.memory = 0),
    initial_fallback = FALSE)
  args <- list(cdf_file = directory, family = quasipoisson("identity"), dispersion = 2,
    precluster = 2L, baseline_df = 4L, nwindows = 1L, correlation = "none", ncores = 1L,
    psf_control = list(kernel_est = kernel, maxiter = 0L), control = controls,
    batch = list(workers = 1L, cache = TRUE,
      coupling = list(initial_cycles = 1L, cycles = 1L, rounds = 1L, q = 12L)))
  reference <- do.call(gcmsDeco, c(args, list(output_dir = outputs[1L])))
  args$ncores <- 2L
  args$batch$cache <- FALSE
  args$control <- list(initial_solver = list(max.memory = 0))
  expect_message(fallback <- do.call(gcmsDeco, c(args, list(output_dir = outputs[2L]))),
    "Eigen CG and diagonal preconditioning")
  expect_s3_class(fallback, "GCMSbatch")
  expect_equal(coef(fallback), coef(reference), tolerance = 1e-8)
  expect_equal(fallback$areas, reference$areas, tolerance = 1e-8)
  expect_equal(GIC(fallback), GIC(reference), tolerance = 1e-8)
  expect_true(all(is.finite(fallback$areas)))
  expect_equal(fallback$parallel, list(ncores = 2L, workers = 1L, threads = 2L))
  expect_equal(fallback$coupled$parallel, list(ncores = 2L, workers = 1L, threads = 2L))
  expect_null(fallback$coupled$store$cache)
  expect_true(any(fallback$coupled$history$stage == "simultaneous_cycle"))
  for (sample in 1:2) {
    fit <- readRDS(fallback$paths[sample])
    expect_equal(dim(fitted(fallback, sample = sample)), c(n, m))
    expect_true(all(is.finite(fitted(fallback, sample = sample))))
    expect_true(file.exists(fit$msp_file) && file.exists(fit$xlsx_file))
    initial <- readRDS(file.path(dirname(fallback$paths[sample]), "initial.rds"))$fit
    expect_equal(initial$model$control$initial_solver$max.memory, 0)
    expect_equal(initial$model$ncores, 2L)
  }
  args$coupled <- FALSE
  expect_message(independent <- do.call(gcmsDeco, c(args, list(output_dir = outputs[3L]))),
    "Eigen CG and diagonal preconditioning")
  expect_null(independent$coupled)
  expect_true(all(file.exists(independent$paths)))
  for (path in independent$paths) {
    fit <- readRDS(path)
    expect_true(all(is.finite(fitted(fit))))
    expect_true(file.exists(fit$msp_file) && file.exists(fit$xlsx_file))
  }
})
