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
  colour <- stats::setNames(scales::hue_pal(l = 48, c = 82)(length(labels)), names(labels))
  fill <- stats::setNames(scales::hue_pal(l = 88, c = 42)(length(labels)), names(labels))
  area_band_colour <- "#d9d9d9"
  reference_colour <- "black"
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

  columns <- min(4L, length(methods) + 1L)
  rows <- ceiling((length(methods) + 1L) / columns)
  draw("simulated-16-area-predictions.png", 600L * columns, 700L * rows, {
    graphics::par(mfrow = c(rows, columns), pty = "s", xaxs = "i", yaxs = "i",
      mar = c(2.9, 1.6, 3.8, 1.0), oma = c(2.7, 3.4, 0.7, 0.5),
      mgp = c(2.15, 0.57, 0), las = 1, tcl = -0.25, cex.axis = 1.2)
    limits <- c(7.5, 10.5)
    ticks <- seq(limits[1L], limits[2L], by = 0.5)
    bands <- lapply(models, function(z) {
      if (is.null(z$fit)) return(NULL)
      grid <- data.frame(log_true = seq(min(z$data$log_true), max(z$data$log_true), length.out = 160L))
      cbind(grid, stats::predict(z$fit, grid, interval = "prediction", level = 0.95))
    })
    for (i in seq_along(methods)) {
      m <- methods[i]
      z <- models[[m]]
      x <- z$data$log_true
      y <- z$data$log_estimated
      interval <- bands[[m]]
      graphics::plot(x, y, type = "n", xlim = limits, ylim = limits,
        xlab = "", ylab = "", axes = FALSE)
      graphics::abline(h = ticks, v = ticks, col = "#e5ebee", lty = 1)
      if (!is.null(interval)) {
        graphics::polygon(c(interval$log_true, rev(interval$log_true)),
          c(interval$lwr, rev(interval$upr)), col = area_band_colour, border = NA)
        graphics::lines(interval$log_true, interval$fit, col = "black", lwd = 2.3)
      }
      graphics::points(x, y, pch = 16, cex = 0.78,
        col = grDevices::adjustcolor(colour[m], alpha.f = 0.5))
      graphics::abline(0, 1, col = reference_colour, lty = 2, lwd = 1.65)
      graphics::axis(1, at = ticks)
      graphics::axis(2, at = ticks, labels = (i - 1L) %% columns == 0L)
      graphics::box()
      outside <- sum(x < limits[1L] | x > limits[2L] | y < limits[1L] | y > limits[2L])
      if (outside > 0L) {
        graphics::legend("topleft", legend = sprintf("%d outside axes", outside),
          bty = "n", text.col = "#485767", cex = 1.0)
      }
      fit_row <- summary[summary$method==m,,drop=FALSE]
      panel_title <- labels[m]
      if (all(c("true_rank","estimated_rank") %in% names(fit_row)) &&
          is.finite(fit_row$estimated_rank) && is.finite(fit_row$true_rank)) {
        panel_title <- sprintf("%s  |  rank %d/%d%s",labels[m],fit_row$estimated_rank,
          fit_row$true_rank,if (isTRUE(fit_row$rank_exact)) " (exact)" else "")
      }
      graphics::mtext(panel_title, 3, 2.05, font = 2)
      digits <- if (isTRUE(z$adjusted_R2 > 0.9995)) 6L else 3L
      graphics::mtext(bquote("Adjusted"~R^2==.(fmt(z$adjusted_R2, digits))~~
        .(sprintf("(%d / %d peaks)", length(x), nrow(truth)))), 3, 0.65, cex = 0.8)
    }
    graphics::plot.new()
    graphics::legend("center", bty = "n", cex = 1.27, y.intersp = 1.45,
      legend = c("Matched peak area", "Least-squares log–log fit", "95% prediction interval",
        "1:1 reference"),
      col = c(grDevices::adjustcolor(colour[1L], alpha.f = 0.5),
        "black", area_band_colour, reference_colour),
      pch = c(16, NA, 15, NA), lty = c(NA, 1, NA, 2), lwd = c(NA, 2.3, NA, 1.65),
      pt.cex = c(1.30, 1, 1.45, 1))
    graphics::mtext(expression(log[10]*"(true peak area, integrated counts)"),
      1, outer = TRUE, line = 0.55, cex = 1.2)
    graphics::mtext(expression(log[10]*"(estimated peak area, integrated counts)"),
      2, outer = TRUE, line = 1.55, cex = 1.2, las = 0)
  })


  benchmark16Similarity(hits,"spectral_cosine",methods,labels,colour,fill,draw,
    "simulated-16-spectral-violin.png","Cosine similarity (logit scale)",
    "Reciprocal spectral hits across 16 samples",stats::median,"Median",
    c(.5,.7,.9,.99,.999,.9999,.99999,.999999),nrow(truth),.7,cap=.999999)
  benchmark16Similarity(hits[hits$present,,drop=FALSE],"spectral_F1",methods,
    labels,colour,fill,draw,"simulated-16-spectral-support.png",
    "Within-spectrum F1 (logit scale)","Correctness of spectral ion support",
    mean,"Mean",c(0,.25,.5,.75,.9,.95,.99,1),nrow(truth))
  benchmark16ProfilePlots(results,methods,labels,draw)
  benchmark16Performance(summary,labels,colour,draw)
  invisible(list(area=models))
}

benchmark16Similarity <- function(hits,metric,methods,labels,colour,fill,draw,
    filename,axis,title,statistic,stat_label,ticks,total,cutoff=NULL,cap=NULL) {
  x<-data.frame(method=hits$method,score=hits[[metric]],
    presence=factor(ifelse(hits$present,"Present","Absent"),levels=c("Present","Absent")))
  x<-x[is.finite(x$score),,drop=FALSE]
  if(!nrow(x))return(invisible(NULL))
  x$method<-factor(x$method,levels=rev(methods))
  stats<-do.call(rbind,lapply(methods,function(m) {
    z<-x[x$method==m,,drop=FALSE]
    q<-stats::quantile(z$score,c(.25,.75),na.rm=TRUE,names=FALSE)
    data.frame(method=m,q25=q[1],q75=q[2],value=if(nrow(z))statistic(z$score)else NA_real_,
      compared=nrow(z),true=sum(z$presence=="Present"))
  }))
  stats$label<-sprintf("%s %.7f · %d/%d true",stat_label,stats$value,stats$true,total)
  display<-function(v)if(is.null(cap))v else pmin(v,cap)
  span<-range(stats::qlogis(display(c(x$score,ticks))),finite=TRUE)
  logit<-function(v) {
    v<-display(v)
    z<-stats::qlogis(v)
    z[which(v==0)]<-span[1]-1.4
    z[which(v==1)]<-span[2]+1.4
    z
  }
  x$position<-logit(x$score)
  limits<-range(c(x$position,logit(ticks)),finite=TRUE)
  plot_limits<-limits+c(-.4,.4)
  tick_labels<-vapply(ticks,function(v)if(!is.null(cap)&&v>=cap)"1"else
    format(v,trim=TRUE,scientific=FALSE),"")
  stats$label_position<-limits[2]+.7
  varying<-methods[vapply(methods,function(m)length(unique(x$score[x$method==m]))>1L,TRUE)]
  plot<-ggplot2::ggplot(x,ggplot2::aes(x=method,y=position))+
    ggplot2::geom_violin(data=x[x$method%in%varying,,drop=FALSE],
      ggplot2::aes(fill=method,colour=method),scale="width",width=.78,adjust=.75,linewidth=.45)+
    ggplot2::geom_point(ggplot2::aes(colour=method,shape=presence),size=1.05,alpha=.5,
      position=ggplot2::position_jitter(width=.23,height=0,seed=23))+
    ggplot2::geom_linerange(data=stats[stats$compared>0L,],
      ggplot2::aes(x=method,ymin=logit(q25),ymax=logit(q75),colour=method),
      inherit.aes=FALSE,linewidth=.85)+
    ggplot2::geom_point(data=stats[stats$compared>0L,],
      ggplot2::aes(x=method,y=logit(value),colour=method),inherit.aes=FALSE,
      shape=21,fill="white",stroke=.74,size=2.6)+
    ggplot2::geom_text(data=stats,ggplot2::aes(x=method,y=label_position,label=label),
      inherit.aes=FALSE,colour="#526579",hjust=0,size=3.2)+
    ggplot2::scale_shape_manual(values=c(Present=16,Absent=4),
      name="Reference compound in the sample",guide=if(any(x$presence=="Absent"))"legend"else"none")+
    ggplot2::scale_fill_manual(values=fill,guide="none")+
    ggplot2::scale_colour_manual(values=colour,guide="none")+
    ggplot2::scale_x_discrete(labels=labels)+
    ggplot2::scale_y_continuous(breaks=logit(ticks),
      labels=tick_labels,limits=plot_limits,oob=scales::oob_keep,
      expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_flip(clip="off")+
    ggplot2::labs(x=NULL,y=axis,title=title,
      subtitle=paste("RI window ±30 · open dot:",tolower(stat_label),"· bar: interquartile range"),
      caption=if(metric=="spectral_F1")
        "Correctly matched true spectra · shared cutoff: 0.1% of base-peak intensity in reference and estimated spectra."
      else "Each point is a reciprocal library hit. Crosses: absent compounds. Display values ≥0.999999 are shown at 1; statistics use original values.")+
    ggplot2::theme_minimal(base_size=13)+
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major.y=ggplot2::element_blank(),legend.position="bottom",
      plot.title=ggplot2::element_text(face="bold",colour="#14273b",size=18),
      plot.subtitle=ggplot2::element_text(colour="#526579",margin=ggplot2::margin(3,0,7,0)),
      plot.caption=ggplot2::element_text(hjust=0,colour="#526579"),
      axis.text.x=ggplot2::element_text(angle=0,hjust=.5),
      axis.title.x=ggplot2::element_text(hjust=.5),
      plot.margin=ggplot2::margin(13,190,9,12))
  if(!is.null(cutoff))plot<-plot+ggplot2::geom_hline(yintercept=logit(cutoff),
    linetype="dotted",colour="#607387",linewidth=.45)
  draw(filename,1900,1050,print(plot))
  invisible(stats)
}

benchmark16ProfilePlots <- function(results,methods,labels,draw) {
  file<-file.path(results,"profiles.rds")
  if(!file.exists(file))return(invisible(NULL))
  x<-readRDS(file)
  quality<-utils::read.csv(file.path(results,"profile-quality.csv"))
  labels<-c("Ground truth"="Ground truth",labels)
  available<-c("Ground truth",methods)
  colours<-grDevices::colorRampPalette(c("#144bb8","#189ad5","#f5c343","#d82b28"))(16)
  for(kind in unique(x$kind)) {
    z<-x[x$kind==kind,,drop=FALSE]
    known<-z[z$method=="Ground truth",,drop=FALSE]
    window<-range(known$ri[known$intensity>max(known$intensity)*1e-4])
    window<-window+c(-.1,.1)*diff(window)
    z<-z[z$ri>=window[1]&z$ri<=window[2],,drop=FALSE]
    known<-z[z$method=="Ground truth",,drop=FALSE]
    ymax<-1.06*max(z$intensity)
    filename<-paste0("simulated-16-",tolower(kind),"-profiles.png")
    draw(filename,2400,1450,{
      graphics::par(mfrow=c(ceiling(length(available)/4),4L))
      graphics::par(mar=c(3.6,4.3,3,.6),oma=c(2.8,0,2.6,0),mgp=c(2.5,.7,0),cex=.95)
      curves<-function(d,lty,lwd) {
        for(rank in seq_len(16)) {
          pieces<-split(d[d$sample_rank==rank,,drop=FALSE],d$component[d$sample_rank==rank])
          for(curve in pieces) {
            curve<-curve[order(curve$ri),,drop=FALSE]
            graphics::lines(curve$ri,curve$intensity,col=colours[rank],lty=lty,lwd=lwd)
          }
        }
      }
      for(method in available) {
        d<-z[z$method==method,,drop=FALSE]
        graphics::plot(NA,xlim=window,ylim=c(0,ymax),xlab="Retention index",
          ylab="Total-ion intensity",main=labels[[method]])
        graphics::grid(col="#e5ebee",lty=1)
        if(method!="Ground truth" && !nrow(d)) {
          graphics::text(mean(window), ymax*.55, "No sample profile exported", col="#526579")
          next
        }
        if(method!="Ground truth")curves(d,1,1.5)
        curves(known,2,1.2)
        q<-quality[quality$kind==kind &quality$method==method,,drop=FALSE]
        if(nrow(q))graphics::mtext(sprintf("Profile cosine %.6f · %d/%d true calls",
          q$mean_profile_cosine,q$recovered,q$present),3,line=.2,cex=.75)
      }
      graphics::mtext(paste("Full chromatograms |",
        if(kind=="Regular")"Highest-abundance regular compound"else
          "Shark-fin-shaped acid compound"),3,outer=TRUE,line=.9,font=2,cex=1.15)
      graphics::mtext("Solid: inferred profiles · Dashed: ground truth · All 16 samples, blue (smallest area) to red (largest area)",
        1,outer=TRUE,line=1.1,cex=.95)
      present<-length(unique(known$sample[known$true_area>0]))
      graphics::mtext(sprintf("1,000 compounds, 10,000 scans per sample. Present in %d of 16 samples. Zero lines show absence. Heights retain total-ion units.",present),
        1,outer=TRUE,line=.05,cex=.9)
    })
  }
}


benchmark16Performance <- function(summary, labels, colour, draw) {
  z <- summary[rev(order(summary$analysis_and_conversion_s_per_sample)), , drop = FALSE]
  draw("simulated-16-performance.png", 1850L, 820L, {
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 8, 3.2, 1.4), las = 1)
    times <- z$analysis_and_conversion_s_per_sample
    ticks <- graphics::barplot(times, horiz = TRUE, names.arg = labels[z$method],
      col = colour[z$method], border = NA, xlim = c(0, max(times) * 1.22),
      xlab = "Wall seconds per sample", main = "Computing time", cex.names = .9)
    graphics::text(times, ticks, sprintf("%.2f", times), pos = 4, cex = .85)
    recall <- 100 * z$TP / (z$TP + z$FN)
    precision <- 100 * z$TP / (z$TP + z$FP)
    y <- seq_len(nrow(z))
    graphics::plot(NA, xlim = c(0, 102), ylim = range(y) + c(-.6, .6), yaxt = "n",
      xlab = "Percent", ylab = "", main = "Occurrence recovery")
    graphics::grid(nx = NULL, ny = NA, col = "#e5ebee", lty = 1)
    graphics::axis(2, at = y, labels = labels[z$method], tick = FALSE, cex.axis = .9)
    graphics::segments(recall, y, precision, y, col = colour[z$method], lwd = 2)
    graphics::points(recall, y, pch = 16, col = colour[z$method], cex = 1.35)
    graphics::points(precision, y, pch = 1, col = "#14273b", cex = 1.35)
    graphics::legend("bottomleft", c("Sensitivity", "Precision"), pch = c(16, 1),
      col = c("#526579", "#14273b"), bty = "n", cex = .85)
  })
}
