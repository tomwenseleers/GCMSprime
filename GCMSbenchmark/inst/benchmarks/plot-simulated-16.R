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

  if (!any(is.finite(hits$spectral_cosine)))
    return(invisible(NULL))
  plot_hits <- data.frame(method=hits$method,cosine=hits$spectral_cosine,
    presence=factor(ifelse(hits$present,"Present","Absent"),levels=c("Present","Absent")))
  plot_hits$method <- factor(plot_hits$method, levels = rev(methods))
  hit_stats <- do.call(rbind, lapply(methods,function(method) {
      z<-plot_hits[plot_hits$method==method,]
      q <- stats::quantile(z$cosine, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
      data.frame(method=method,q25=q[1L],median=q[2L],q75=q[3L],
        compared=sum(is.finite(z$cosine)),true=sum(z$presence=="Present"),total=nrow(truth))
    }))
  hit_stats$label <- sprintf("%s · %d/%d true",fmt(hit_stats$median,6),
    hit_stats$true,hit_stats$total)
  plot_hits <- plot_hits[is.finite(plot_hits$cosine), ]
  ticks <- c(0.5,0.7,0.9,0.99,0.999,0.9999,0.99999,0.999999)
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
  hit_stats$label_position <- limits[2L] + 0.7
  violin <- ggplot2::ggplot(plot_hits, ggplot2::aes(x = method, y = position)) +
    ggplot2::geom_violin(ggplot2::aes(fill = method, colour = method),
      scale = "width", width = 0.78, adjust = 0.75, linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(colour=method,shape=presence),size=1.1,alpha=.5,
      position = ggplot2::position_jitter(width = 0.23, height = 0, seed = 23)) +
    ggplot2::geom_hline(yintercept = logit(0.70), linetype = "dotted",
      colour = "#607387", linewidth = 0.45) +
    ggplot2::geom_linerange(data = hit_stats[hit_stats$compared > 0L, ],
      ggplot2::aes(x = method, ymin = logit(q25), ymax = logit(q75), colour = method),
      inherit.aes = FALSE, linewidth = 0.85) +
    ggplot2::geom_point(data = hit_stats[hit_stats$compared > 0L, ],
      ggplot2::aes(x = method, y = logit(median), colour = method),
      inherit.aes = FALSE, shape = 21, fill = "white", stroke = 0.74, size = 2.6) +
    ggplot2::geom_text(data = hit_stats,
      ggplot2::aes(x = method, y = label_position, label = label),
      inherit.aes = FALSE, colour = "#526579", hjust = 0, size = 3.35) +
    ggplot2::scale_shape_manual(values=c(Present=16,Absent=4),
      name="Reference compound in the sample") +
    ggplot2::scale_fill_manual(values = fill, guide = "none") +
    ggplot2::scale_colour_manual(values = colour, guide = "none") +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_y_continuous(breaks = logit(ticks),
      labels = sub("\\.?0+$", "", format(ticks, trim = TRUE, scientific = FALSE)),
      limits = limits + c(-0.4, 6), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Cosine similarity (logit scale)",
      title="Reciprocal spectral hits across 16 simulations",
      subtitle="RI window ±30 · open dot and bar show the median and interquartile range",
      caption="Each point is a reciprocal library hit. Crosses identify compounds absent from the sample.") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(), panel.spacing.x = grid::unit(1.4, "lines"),
      strip.text = ggplot2::element_text(face = "bold", colour = "#14273b",
        margin = ggplot2::margin(5, 0, 8, 0)),
      plot.title = ggplot2::element_text(face = "bold", colour = "#14273b", size = 18),
      plot.subtitle = ggplot2::element_text(colour = "#526579", margin = ggplot2::margin(3, 0, 7, 0)),
      legend.position="bottom",
      plot.caption = ggplot2::element_text(hjust = 0, colour = "#526579"),
      plot.margin = ggplot2::margin(13, 18, 9, 12))
  draw("simulated-16-spectral-violin.png",1800,1050,print(violin))
  benchmark16Heatmap(results,output_dir,draw)
  benchmark16AcidPlot(results,methods,labels,draw)
  invisible(list(area = models, spectral = violin))
}

benchmark16Heatmap <- function(results,output_dir,draw) {
  file <- file.path(results,"sample-heatmap.rds")
  if(!file.exists(file)) return(invisible(NULL))
  x <- readRDS(file)
  grid <- seq(min(x$time),max(x$time),length.out=length(x$time))
  matrices <- lapply(list(x$observed,x$mean),function(z)
    apply(z,2L,function(v)stats::approx(x$time,v,grid)$y))
  palette <- grDevices::colorRampPalette(
    c("black","black","blue","green","yellow","orange","red"))(1000L)
  lim <- range(sqrt(unlist(matrices)))
  draw("simulated-16-sample-heatmap.png",2200,1250,{
    graphics::par(mfrow=c(2,1),mar=c(3.3,4.2,2.5,1),mgp=c(2.3,.7,0),las=1)
    titles <- c("Observed counts","True mean, including smooth background")
    for(i in 1:2) {
      graphics::image(grid,x$mass,sqrt(matrices[[i]]),col=palette,zlim=lim,
        useRaster=TRUE,xlab="Retention time (min)",ylab="Mass-to-charge ratio",
        main=paste(x$sample,titles[i],sep=" · "))
      graphics::abline(v=x$apex[x$present],col=grDevices::adjustcolor("red",.55),lwd=.7)
      graphics::abline(v=x$apex[intersect(which(x$present),x$acid)],col="red",lwd=1.8)
    }
  })
}

benchmark16AcidPlot <- function(results,methods,labels,draw) {
  file <- file.path(results,"acid-profiles.csv")
  if(!file.exists(file)) return(invisible(NULL))
  x <- utils::read.csv(file)
  available <- intersect(methods,unique(x$method))
  if(!length(available)) return(invisible(NULL))
  known <- x[x$method=="Ground truth",]
  truth <- do.call(rbind,lapply(available,function(m) { z<-known;z$method<-m;z }))
  estimated <- x[x$method%in%available,]
  estimated$kind <- "Recovered"
  truth$kind <- "True"
  z <- rbind(estimated,truth)
  z$method <- factor(z$method,levels=available)
  z$curve <- interaction(z$method,z$sample,z$component,z$kind,drop=TRUE)
  z$emphasis <- ifelse(z$kind=="True"|z$reciprocal,.75,.35)
  plot <- ggplot2::ggplot(z,ggplot2::aes(x=ri,y=intensity,colour=sample,
    linetype=kind,group=curve,linewidth=emphasis)) +
    ggplot2::geom_line(alpha=.8) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_linetype_manual(values=c(True="dashed",Recovered="solid"),name=NULL) +
    ggplot2::facet_wrap(~method,ncol=2,labeller=ggplot2::as_labeller(labels)) +
    ggplot2::labs(x="Retention index",y="Relative peak height",
      title="One broad acid-like compound across samples",
      subtitle="Dashed curves show truth. Solid curves show each compatible returned component.",
      colour="Sample",
      caption="Thin curves are additional calls to the same reference. Profiles use each engine's native export.") +
    ggplot2::theme_minimal(base_size=13) +
    ggplot2::theme(legend.position="bottom",strip.text=ggplot2::element_text(face="bold"),
      panel.grid.minor=ggplot2::element_blank())
  draw("simulated-16-acid-profiles.png",2100,1500,print(plot))
}
