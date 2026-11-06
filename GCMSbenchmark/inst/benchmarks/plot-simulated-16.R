plotBenchmark16 <- function(results, output_dir = results) {
  read_table <- function(file) utils::read.csv(file.path(results, file), check.names = FALSE)
  summary <- read_table("benchmark.csv")
  matches <- read_table("matches.csv")
  hits <- read_table("spectral-hits.csv")
  truth <- read_table("truth-peaks.csv")
  methods <- summary$method[summary$status == "ok"]
  if (!length(methods)) stop("No successful benchmark runs to plot.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  labels <- c(GCMSprime = "GCMSprime", "MS-DIAL" = "MS-DIAL", mzmine = "mzmine GC-EI",
    eRah = "eRah", AMDIS = "AMDIS", "xcms+CAMERA" = "xcms + CAMERA", MSHub = "MSHub")
  units <- c(GCMSprime = "counts", "MS-DIAL" = "model-ion area",
    mzmine = "representative-ion area", eRah = "eRah area", AMDIS = "counts",
    "xcms+CAMERA" = "count-seconds", MSHub = "aligned NMF integral")
  colour <- stats::setNames(scales::hue_pal(l = 48, c = 82)(length(labels)), names(labels))
  fill <- stats::setNames(scales::hue_pal(l = 88, c = 42)(length(labels)), names(labels))
  reference_colour <- "#26364b"
  fmt <- function(x, digits) ifelse(is.na(x), "–", formatC(x, format = "f", digits = digits))
  draw <- function(name, width, height, expr) {
    grDevices::png(file.path(output_dir, name), width = width, height = height, res = 155)
    on.exit(grDevices::dev.off())
    force(expr)
  }
  area_regression <- utils::getFromNamespace("area_regression", "GCMSbenchmark")
  models <- stats::setNames(lapply(methods, function(m) {
    area_regression(matches[matches$method == m, , drop = FALSE])
  }), methods)

  draw("simulated-16-area-predictions.png", 2390, 1440, {
    graphics::par(mfrow = c(ceiling((length(methods) + 1L) / 4L), 4L),
      mar = c(3, 3.7, 3.8, 0.85), oma = c(2.3, 0.5, 0.7, 0.5),
      mgp = c(2.15, 0.57, 0), las = 1, tcl = -0.25, cex.axis = 0.9)
    xlim <- range(log10(truth$area_true)) + c(-0.06, 0.06)
    for (m in methods) {
      z <- models[[m]]
      x <- z$data$log_true
      y <- z$data$log_estimated
      interval <- if (!is.null(z$fit)) {
        grid <- seq(min(x), max(x), length.out = 160L)
        stats::predict(z$fit, newdata = data.frame(log_true = grid), interval = "prediction")
      } else NULL
      counts <- identical(summary$native_area[match(m, summary$method)], "integrated counts")
      ylim <- range(c(y, interval, if (counts) xlim, if (!length(y)) c(0, 1)))
      graphics::plot(x, y, type = "n", xlim = xlim, ylim = ylim,
        xlab = "", ylab = bquote(log[10]~.(unname(units[m]))))
      graphics::grid(col = "#e5ebee", lty = 1)
      if (!is.null(interval)) {
        graphics::polygon(c(grid, rev(grid)), c(interval[, "lwr"], rev(interval[, "upr"])),
          col = fill[m], border = NA)
      }
      graphics::points(x, y, pch = 16, cex = 0.64,
        col = grDevices::adjustcolor(colour[m], alpha.f = 0.66))
      if (!is.null(interval))
        graphics::lines(grid, interval[, "fit"], col = colour[m], lwd = 2.3)
      if (counts) graphics::abline(0, 1, col = reference_colour, lty = 2, lwd = 1.65)
      graphics::mtext(labels[m], 3, 2.05, font = 2)
      digits <- if (isTRUE(z$adjusted_R2 > 0.9995)) 5L else 3L
      graphics::mtext(bquote("Adjusted"~R^2==.(fmt(z$adjusted_R2, digits))~~
        .(sprintf("(%d / %d peaks)", length(x), nrow(truth)))), 3, 0.65, cex = 0.8)
    }
    graphics::plot.new()
    graphics::legend("center", bty = "n", cex = 1.27, y.intersp = 1.45,
      legend = c("Matched peak area", "Log–log fit", "95% prediction interval",
        "1:1 line in count units"),
      col = c(colour[1L], colour[1L], fill[1L], reference_colour),
      pch = c(16, NA, 15, NA), lty = c(NA, 1, NA, 2), lwd = c(NA, 2.3, NA, 1.65),
      pt.cex = c(1.06, 1, 1.45, 1))
    graphics::mtext(expression(log[10]*"(True peak area, integrated counts)"),
      1, outer = TRUE, line = 0.6, cex = 1.1)
  })

  if (!any(is.finite(c(hits$best_cosine, hits$near_cosine))))
    return(invisible(NULL))
  plot_hits <- rbind(
    data.frame(method = hits$method, window = "Best spectrum anywhere in the sample",
      cosine = hits$best_cosine),
    data.frame(method = hits$method, window = "Best spectrum within ±7.2 seconds",
      cosine = hits$near_cosine))
  plot_hits$method <- factor(plot_hits$method, levels = rev(methods))
  plot_hits$window <- factor(plot_hits$window, levels = unique(plot_hits$window))
  stats <- do.call(rbind, lapply(split(plot_hits, plot_hits[c("method", "window")], drop = TRUE),
    function(z) {
      q <- stats::quantile(z$cosine, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
      data.frame(method = z$method[1L], window = z$window[1L], q25 = q[1L],
        median = q[2L], q75 = q[3L], compared = sum(is.finite(z$cosine)), total = nrow(z))
    }))
  stats$label <- sprintf("%s · %d/%d", fmt(stats$median, 6), stats$compared, stats$total)
  plot_hits <- plot_hits[is.finite(plot_hits$cosine), ]
  ticks <- c(0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999)
  span <- range(stats::qlogis(c(plot_hits$cosine, ticks)), finite = TRUE)
  logit <- function(x) {
    z <- stats::qlogis(x)
    # Display exact endpoints just beyond the finite logits, if they occur.
    z[which(x == 0)] <- span[1L] - 1.4
    z[which(x == 1)] <- span[2L] + 1.4
    z
  }
  ticks <- sort(unique(c(ticks, intersect(plot_hits$cosine, c(0, 1)))))
  plot_hits$position <- logit(plot_hits$cosine)
  limits <- range(c(plot_hits$position, logit(ticks)), finite = TRUE)
  stats$label_position <- limits[2L] + 0.7
  violin <- ggplot2::ggplot(plot_hits, ggplot2::aes(x = method, y = position)) +
    ggplot2::geom_violin(ggplot2::aes(fill = method, colour = method),
      scale = "width", width = 0.78, adjust = 0.75, linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(colour = method), size = 0.84, alpha = 0.46,
      position = ggplot2::position_jitter(width = 0.23, height = 0, seed = 23)) +
    ggplot2::geom_hline(yintercept = logit(0.70), linetype = "dotted",
      colour = "#607387", linewidth = 0.45) +
    ggplot2::geom_linerange(data = stats[stats$compared > 0L, ],
      ggplot2::aes(x = method, ymin = logit(q25), ymax = logit(q75), colour = method),
      inherit.aes = FALSE, linewidth = 0.85) +
    ggplot2::geom_point(data = stats[stats$compared > 0L, ],
      ggplot2::aes(x = method, y = logit(median), colour = method),
      inherit.aes = FALSE, shape = 21, fill = "white", stroke = 0.74, size = 2.6) +
    ggplot2::geom_text(data = stats,
      ggplot2::aes(x = method, y = label_position, label = label),
      inherit.aes = FALSE, colour = "#526579", hjust = 0, size = 3.35) +
    ggplot2::facet_wrap(~window, nrow = 1) +
    ggplot2::scale_fill_manual(values = fill, guide = "none") +
    ggplot2::scale_colour_manual(values = colour, guide = "none") +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_y_continuous(breaks = logit(ticks),
      labels = sub("\\.?0+$", "", format(ticks, trim = TRUE, scientific = FALSE)),
      limits = limits + c(-0.4, 6), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Cosine similarity (logit scale)",
      title = "Best spectral hits across the 16 simulations",
      subtitle = "Small dots show individual best hits. The open dot and bar show the median and interquartile range.",
      caption = "Labels show median cosine and candidate coverage. The dotted line marks the 0.70 cosine cutoff.") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(), panel.spacing.x = grid::unit(1.4, "lines"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#14273b",
        margin = ggplot2::margin(5, 0, 8, 0)),
      plot.title = ggplot2::element_text(face = "bold", colour = "#14273b", size = 18),
      plot.subtitle = ggplot2::element_text(colour = "#526579", margin = ggplot2::margin(3, 0, 7, 0)),
      plot.caption = ggplot2::element_text(hjust = 0, colour = "#526579"),
      plot.margin = ggplot2::margin(13, 18, 9, 12))
  draw("simulated-16-spectral-violin.png", 2480, 1100, print(violin))
  invisible(list(area = models, spectral = violin))
}
