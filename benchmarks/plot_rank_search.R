# Recreate the two-size rank and timing comparison from the scored CDF runs.
plotRankSearch <- function(data = "benchmarks/results/rank-search.csv",
                           file = "rank-search-comparison.png") {
  x <- utils::read.csv(data)
  stopifnot(nrow(x) == 4L, all(is.finite(x$seconds)), all(x$seconds > 0))
  grDevices::png(file, width = 1720, height = 790, res = 165)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.4, 5.5, 2.5, .6), oma = c(.3, 0, .8, .4),
                las = 1, family = "sans", xaxs = "i")
  for (size in c("small", "full")) {
    z <- x[x$size == size, ]
    z <- z[match(c("Default", "Local search"), z$method), ]
    y <- c(2.1, 1.05)
    minutes <- z$seconds / 60
    ticks <- pretty(c(0, max(minutes)), n = 4)
    limit <- max(ticks, minutes) * 1.07
    graphics::plot.new()
    graphics::plot.window(c(0, limit), c(.45, 2.65))
    graphics::title(main = sprintf("%s scans per sample", format(z$scans[1], big.mark = ",")),
                    xlab = "CDF input through export (minutes)", cex.main = 1.12)
    graphics::axis(1, at = ticks[ticks >= 0 & ticks <= limit])
    graphics::axis(2, at = y, labels = z$method, tick = FALSE, line = -.2, cex.axis = .85)
    graphics::abline(v = ticks, col = "#e7eaed", lwd = .8)
    cols <- c("#74818e", "#087e83")
    graphics::segments(0, y, minutes, y, col = cols, lwd = 2)
    graphics::points(minutes, y, pch = 16, col = cols, cex = 1.25)
    rank <- sprintf("Estimated rank %s / true %s%s", format(z$estimated_rank, trim=TRUE, big.mark=","),
                    format(z$true_rank, trim=TRUE, big.mark=","),
                    ifelse(z$estimated_rank == z$true_rank, " (exact)", ""))
    graphics::text(0, y + .23, rank, adj = 0, col = "#203447", cex = .86)
    graphics::text(0, y - .25, sprintf("%.1f s  |  %s TP, %s FP", z$seconds,
                  format(z$TP, trim = TRUE, big.mark = ","), z$FP), adj = 0, cex = .82)
  }
  invisible(file)
}
