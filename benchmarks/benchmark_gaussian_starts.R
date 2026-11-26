# Fit the public 16-sample example at a fixed Gaussian full width at half maximum.
# Run from the GCMSprime repository. The two benchmark widths are 3.64015 and
# 14.56059 scans, relative to the 7.28030-scan median for small regular peaks.
benchmarkGaussianStart <- function(input, output_dir, fwhm_scans, ncores = 8L) {
  stopifnot(length(fwhm_scans) == 1L, is.finite(fwhm_scans), fwhm_scans > 0)
  if (dir.exists(output_dir)) stop("Choose a fresh output directory.")
  time <- lapply(input$files, function(path) {
    con <- ncdf4::nc_open(path)
    on.exit(ncdf4::nc_close(con))
    stopifnot(tolower(ncdf4::ncatt_get(con, "scan_acquisition_time", "units")$value) %in% c("seconds", "second"))
    as.numeric(ncdf4::ncvar_get(con, "scan_acquisition_time")) / 60
  })
  sigma <- fwhm_scans / sqrt(8 * log(2))
  radius <- ceiling(6 * sigma)
  offset <- seq.int(-radius, radius)
  kernel <- exp(-0.5 * (offset / sigma)^2)
  kernel <- kernel / sum(kernel)
  spacing <- median(vapply(time, function(t) median(diff(t)), numeric(1)))
  standard <- list(calibration = list(psf_offsets = spacing * offset,
    peak_timepoints = range(unlist(time)), psf_kernels = rbind(kernel, kernel)))
  set.seed(9118L)
  GCMSbenchmark::GCMSbenchmark(input$files, output_dir = output_dir, method = "GCMSprime",
    truth = input$truth, reference = input$reference, matching = "reciprocal",
    ri_tolerance = 30, rt_tolerance = 0.12, min_cosine = 0.70, ion_threshold = 0.001,
    ncores = ncores, precluster = 32L, nwindows = 2L, baseline_df = 4L,
    family = stats::quasipoisson("identity"), correlation = "none",
    psf_control = list(standard = standard, maxiter = 0L),
    batch = list(workers = 4L, coupling = list(alignment_window = 0.5, penalty_search = TRUE)))
}
# source(system.file("benchmarks", "simulated-16.R", package = "GCMSbenchmark"))
# input <- simulateBenchmark16("gaussian-start-input")
# half <- benchmarkGaussianStart(input, "gaussian-half", 7.280297 / 2)
# double <- benchmarkGaussianStart(input, "gaussian-double", 7.280297 * 2)
