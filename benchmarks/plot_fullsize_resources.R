# Plot a monitored CDF workflow. Published tables can be used without the raw data.
plotFullsize16Resources <- function(directory = "benchmarks/results/fullsize16",
                                    file = NULL, physical_cores = NULL, logical_cores = NULL) {
  x <- utils::read.csv(file.path(directory, "resources.csv"))
  times <- utils::read.csv(file.path(directory, "timings.csv"))
  end <- times$total
  # The published monitor ends when the fit is saved, just after the timed export.
  # Align elapsed stage times to that end point to account for process startup.
  stopifnot(nrow(x) > 1L, all(diff(x$seconds) > 0))
  startup <- max(0, max(x$seconds) - end)
  cuts <- c(0, startup + c(times$initial, times$initial + times$psf_regressions,
                            times$initial + times$coupled, end)) / 60
  stages <- c("CDFs and calibration", "Grouped PSF", "Alignment and joint fit", "Export")
  xl <- range(cuts)
  if (is.null(file)) file <- file.path(directory, "fullsize-16-resources.png")
  grDevices::png(file, width = 1900, height = 1260, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(1:3, 3, 1), heights = c(.88, 2.9, 3.1))
  graphics::par(mar = c(0, 4.8, 2.6, 1), xaxs = "i", family = "sans")
  graphics::plot.new()
  graphics::plot.window(xl, c(0, 1))
  graphics::title(main = "CPU and memory during deconvolution of 16 full-size chromatograms",
                  line = .8, cex.main = 1.18)
  left <- head(cuts, -1)
  right <- tail(cuts, -1)
  graphics::rect(left, .02, right, .8, col = c("#d6e6f3", "#b0d5da", "#d4e6db", "#f0d9b8"), border = "white")
  graphics::text((left + right) / 2, .4, stages, cex = .76)
  graphics::par(mar = c(.6, 4.8, 1.15, 1), las = 1)
  t <- x$seconds
  graphics::plot(0, type = "n", xlim = xl,
                 ylim = c(0, max(x$available_GiB, x$working_GiB, x$private_GiB, na.rm = TRUE) * 1.08),
                 xlab = "", ylab = "Memory (GiB)", xaxt = "n", bty = "l")
  graphics::abline(h = graphics::axTicks(2), col = "#e5eaee")
  graphics::abline(v = cuts[2:4], col = "#a3b1be", lty = 3)
  available <- "available_GiB" %in% names(x) && any(is.finite(x$available_GiB))
  if (available) graphics::lines(t / 60, x$available_GiB, col = "#71818c", lwd = 1.2, lty = 2)
  graphics::lines(t / 60, x$working_GiB, col = "#008b91", lwd = 1.9)
  graphics::lines(t / 60, x$private_GiB, col = "#163e63", lwd = 1.15)
  series <- c(available, TRUE, TRUE)
  graphics::legend("right", c("System RAM available", "Process-tree working set", "Process-tree private memory")[series],
                   col = c("#71818c", "#008b91", "#163e63")[series], lty = c(2, 1, 1)[series], lwd = c(1.2, 1.9, 1.15)[series],
                   bty = "n", cex = .83)
  graphics::par(mar = c(4.6, 4.8, 1.15, 1), mgp = c(2.6, .7, 0))
  # Interpolate cumulative CPU time so irregular monitor intervals get their actual duration.
  cpu <- (x$cpu_seconds - stats::approx(t, x$cpu_seconds, t - 15)$y) / 15
  hardware <- c(physical_cores, logical_cores)
  graphics::plot(t / 60, cpu, type = "n", xlim = xl, ylim = c(0, 1.07 * max(1, cpu, hardware, na.rm = TRUE)),
                 xlab = "Elapsed time (minutes)", ylab = "Logical-core equivalents", bty = "l")
  graphics::abline(h = graphics::axTicks(2), col = "#e5eaee")
  if (length(hardware)) graphics::abline(h = hardware, col = "#697780", lty = 2)
  graphics::abline(v = cuts[2:4], col = "#a3b1be", lty = 3)
  graphics::lines(t / 60, cpu, col = "#164d78", lwd = 1.65)
  model <- c(if (length(physical_cores)) paste(physical_cores, "physical cores"),
             if (length(logical_cores)) paste(logical_cores, "logical cores"))
  graphics::legend("topright", c("CPU, 15-second average", model),
                   col = c("#164d78", rep("#697780", length(model))),
                   lty = c(1L, rep(2L, length(model))), lwd = c(1.65, rep(1, length(model))),
                   bty = "o", box.lty = 0, bg = "white", cex = .83)
  invisible(file)
}
