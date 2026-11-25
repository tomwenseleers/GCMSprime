# Recreate the two diagnostic panels from the nine fitted models.
plotPSFStarts <- function(results, file) {
  d <- utils::read.csv(file.path(results, "psf-starts.csv"))
  starts <- unique(d$initial_psf)
  stages <- unique(d$search)
  cols <- c("#1968a6", "#d17920", "#704a86")
  offsets <- c(-0.06, 0, 0.06)
  grDevices::png(file, width = 1650, height = 720, res = 145)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, 2), mar = c(4.8, 4.2, 2.2, 1.4),
                oma = c(1.0, 0.2, 1.4, 0.2), mgp = c(2.3, 0.65, 0), las = 1)
  for (metric in c("FN", "FP")) {
    ylim <- if (metric == "FN") c(-0.4, 7) else c(-7, 140)
    graphics::plot(NA, xlim = c(0.75, 3.25), ylim = ylim, xaxt = "n",
      bty = "l", xlab = "", ylab = "Sample-level occurrences",
      main = if (metric == "FN") "Missed true peaks" else "False positive calls")
    graphics::abline(h = graphics::axTicks(2), col = "#e3e9ee")
    graphics::axis(1, at = 1:3, labels = c("Full\nsweeps", "Fourfold\nsearch", "Search + blind\nPSF exchange"),
                   tick = FALSE, line = 0.45)
    for (i in seq_along(starts)) {
      arm <- d[d$initial_psf == starts[i], ]
      values <- arm[[metric]][match(stages, arm$search)]
      graphics::lines(seq_along(stages) + offsets[i], values, type = "o",
                      col = cols[i], lwd = 2.2, pch = 14 + i)
    }
    if (metric == "FN")
      graphics::legend("topright", starts, col = cols, pch = 15:17,
                       lwd = 2, bty = "n", cex = 0.81)
  }
  exact <- d[d$search==utils::tail(stages,1),"rank"]
  graphics::mtext(sprintf("True rank %d   |   Final estimated ranks %s%s",
                  unique(d$true_rank),paste(exact,collapse=", "),
                  if (all(exact==unique(d$true_rank))) " (all exact)" else ""),
                  outer = TRUE, side = 3, line = 0.1, font = 2)
  invisible(file)
}
