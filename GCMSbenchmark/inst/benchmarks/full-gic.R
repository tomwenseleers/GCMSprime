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

# Draw both measured trajectories from the block-level CSV returned by the diagnostic.
plotFullBenchmarkGIC <- function(trace, file) {
  arms <- lapply(c("standard", "true profiles"), function(s) trace[trace$start == s, , drop = FALSE])
  stopifnot(all(vapply(arms, nrow, 0L) > 0L), all(is.finite(trace$GIC)))
  reference <- tail(arms[[2]]$GIC, 1L)
  blind <- tail(arms[[1]]$GIC, 1L)
  stopifnot(reference != 0)
  delta <- function(gic) 100 * (gic - reference) / abs(reference)
  colours <- c("#17649a", "#a8521b")
  baseline <- "#c32d35"
  early <- arms[[2]][arms[[2]]$stage == "matched reference continuation" & arms[[2]]$block == "start", , drop = FALSE]
  initial <- range(trace$GIC / 1e6)
  grDevices::png(file, width = 1700, height = 790, res = 160)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, 2), mar = c(3.8, 4.7, 2.0, 1.0),
    oma = c(2.2, 0.2, 2.5, 0.3), mgp = c(2.9, .75, 0), las = 1,
    family = "sans", fg = "#34485a", col.lab = "#23384a", col.axis = "#34485a")
  for (zoom in c(FALSE, TRUE)) {
    lines <- if (zoom) lapply(arms, function(z) z[z$sweep >= max(z$sweep)/2, , drop = FALSE]) else arms
    y <- lapply(lines, function(z) if (zoom) delta(z$GIC) else z$GIC / 1e6)
    x <- lapply(lines, "[[", "sweep")
    levels <- if (zoom) range(c(0, unlist(y))) else initial
    if (diff(levels) == 0) levels <- levels + c(-1, 1) * 0.01
    graphics::plot(range(unlist(x)), levels, type = "n", axes = FALSE,
      xlab = "Profile and spectrum sweeps",
      ylab = if (zoom) "GIC difference from reference (%)" else "Original-ion GIC (millions)",
      main = if (zoom) "Final half of each fit" else "All recorded updates")
    graphics::grid(col = "#e1e9ed", lty = 1)
    graphics::axis(1)
    graphics::axis(2)
    graphics::box(col = "#9cabb6")
    graphics::abline(h = if (zoom) 0 else reference / 1e6, col = baseline, lty = 3, lwd = 2)
    for (i in seq_along(lines)) {
      graphics::lines(x[[i]], y[[i]], col = colours[i], lty = i, lwd = 2.2)
      graphics::points(tail(x[[i]], 1L), tail(y[[i]], 1L), col = colours[i], pch = 16, cex = .72)
    }
    if (nrow(early)) graphics::points(early$sweep[1L],
      if (zoom) delta(early$GIC[1L]) else early$GIC[1L] / 1e6,
      pch = 21, col = colours[2L], bg = "white", lwd = 1.5, cex = .85)
    if (!zoom) {
      labels <- c("Blind start", "Ground truth-initialised", "Reference endpoint")
      symbols <- rep(NA_integer_, 3L)
      types <- 1:3
      inks <- c(colours, baseline)
      if (nrow(early)) {
        labels <- c(labels, "Reference met stopping rule")
        symbols <- c(symbols, 21L)
        types <- c(types, NA_integer_)
        inks <- c(inks, colours[2L])
      }
      graphics::legend("topright", bty = "n", cex = .8, legend = labels,
        pch = symbols, pt.bg = "white", lty = types, lwd = 2, col = inks)
    }
  }
  graphics::mtext("GIC convergence on 16 full chromatograms", side = 3, outer = TRUE,
    font = 2, cex = 1.32, col = "#14273b", line = .6)
  graphics::mtext(sprintf("Blind GIC %s    |    Reference GIC %s    |    Final signed gap %+.4f%%",
    format(round(blind), big.mark = ","), format(round(reference), big.mark = ","), delta(blind)),
    side = 1, outer = TRUE, line = .55, cex = .92, col = "#354e63")
  invisible(file)
}
