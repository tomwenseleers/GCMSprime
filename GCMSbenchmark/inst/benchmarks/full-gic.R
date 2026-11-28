# Re-run the calibrated full benchmark with both objective traces.
runFullBenchmarkGIC <- function(directory, output_dir, ncores = 16L) {
  if (utils::packageVersion("logconNMF") < "0.1.0.9015")
    stop("The matched-sweep diagnostic requires logconNMF 0.1.0.9015 or later.")
  directory <- normalizePath(directory, winslash = "/")
  truth <- readRDS(file.path(directory, "truth.rds"))
  profiles <- lapply(truth, function(x) x[c("id", "time", "H")])
  rm(truth)
  files <- sort(list.files(file.path(directory, "input"), "[.]cdf$", full.names = TRUE))
  stopifnot(length(files) == 16L, identical(names(profiles), tools::file_path_sans_ext(basename(files))))
  set.seed(171117L)
  started <- proc.time()[["elapsed"]]
  raw <- utils::getFromNamespace("GCMS_read_cdf", "GCMSprime")(file.path(directory, "standard", "psf-standard.cdf"), NULL)
  standard <- logconNMF::estimatePSF(raw$Y, time = raw$time, family = quasipoisson("identity"),
    method = "standard", nwindows = 10L, precluster = 1L, ncores = min(8L, ncores))
  rm(raw)
  fit <- GCMSprime::gcmsDeco(dirname(files[1L]), output_dir = file.path(output_dir, "analysis"),
    ncores = ncores, psf_control = list(standard = standard),
    batch = list(workers = min(8L, ncores), coupling = list(penalty_search = TRUE)),
    truth.profiles = profiles, objective.trace = TRUE, verbose = TRUE)
  elapsed <- proc.time()[["elapsed"]] - started
  gap <- logconNMF::gicGap(fit)
  saveRDS(gap, file.path(output_dir, "gic-gap.rds"))
  utils::write.csv(gap$trace, file.path(output_dir, "gic-trace.csv"), row.names = FALSE)
  summary <- data.frame(true_rank = length(unique(unlist(lapply(profiles, "[[", "id")))),
    estimated_rank = fit$coupled$rank, fitted_gic = gap$fitted_gic,
    reference_gic = gap$reference_gic, difference = gap$difference, percent = gap$percent,
    diagnostic_seconds = elapsed, reference_seconds = unname(gap$reference_seconds))
  utils::write.csv(summary, file.path(output_dir, "gic-gap.csv"), row.names = FALSE)
  print(summary, digits = 12)
  invisible(gap)
}

# Plot the blind fit at its recorded checkpoints in their original order.
plotFullBenchmarkGIC <- function(trace, file) {
  blind <- trace[trace$start == "standard", c("sweep", "GIC"), drop = FALSE]
  stopifnot(nrow(blind) > 0L, all(is.finite(blind$GIC)),
            all(is.finite(blind$sweep)), !is.unsorted(blind$sweep))
  gic <- blind$GIC / 1e6
  ink <- "#19354a"
  blue <- "#17649a"

  grDevices::png(file, width = 1460, height = 760, res = 155)
  on.exit(grDevices::dev.off())
  graphics::par(mar = c(4.0, 4.6, 4.2, 1.1), mgp = c(2.6, 0.7, 0),
                las = 1, family = "sans", col.lab = ink, col.axis = ink)
  graphics::plot(blind$sweep, gic, type = "n", axes = FALSE,
                 xlab = "Profile and spectrum sweeps", ylab = "Original-ion GIC (millions)")
  xticks <- pretty(blind$sweep, n = 8)
  yticks <- pretty(gic)
  graphics::abline(v = xticks, h = yticks, col = "#e5ebef", lwd = 0.8)
  graphics::axis(1, at = xticks, tck = -0.012, col = "#a9b7bf")
  graphics::axis(2, at = yticks, tck = -0.012, col = "#a9b7bf")
  graphics::box(col = "#afbbc3")
  graphics::lines(blind$sweep, gic, col = blue, lwd = 2.4)
  graphics::points(blind$sweep, gic, pch = 16, cex = 0.62, col = blue)
  graphics::title(main = "GIC convergence of the blind fit", line = 2.25,
                  adj = 0, cex.main = 1.45, col.main = ink)
  graphics::mtext(sprintf("16 full chromatograms  |  %d recorded checkpoints", nrow(blind)),
                  side = 3, line = 0.7, adj = 0, col = "#475e70", cex = 0.96)
  invisible(file)
}
