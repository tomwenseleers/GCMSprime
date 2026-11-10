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
  units <- stats::setNames(rep("Total-ion counts",length(labels)),names(labels))
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
      ylim <- range(c(y, interval, xlim))
      graphics::plot(x, y, type = "n", xlim = xlim, ylim = ylim,
        xlab = "", ylab = bquote(log[10]~.(unname(units[m]))))
      graphics::grid(col = "#e5ebee", lty = 1)
      if (!is.null(interval)) {
        graphics::polygon(c(grid, rev(grid)), c(interval[, "lwr"], rev(interval[, "upr"])),
          col = area_band_colour, border = NA)
      }
      if (!is.null(interval))
        graphics::lines(grid, interval[, "fit"], col = "black", lwd = 2.3)
      graphics::abline(0, 1, col = reference_colour, lty = 2, lwd = 1.65)
      graphics::points(x, y, pch = 16, cex = 0.78,
        col = grDevices::adjustcolor(colour[m], alpha.f = 0.5))
      graphics::mtext(labels[m], 3, 2.05, font = 2)
      digits <- if (isTRUE(z$adjusted_R2 > 0.9995)) 5L else 3L
      graphics::mtext(bquote("Robust adjusted"~R^2==.(fmt(z$adjusted_R2, digits))~~
        .(sprintf("(%d / %d peaks)", length(x), nrow(truth)))), 3, 0.65, cex = 0.8)
    }
    graphics::plot.new()
    graphics::legend("center", bty = "n", cex = 1.27, y.intersp = 1.45,
      legend = c("Matched peak area", "Robust log–log fit", "95% robust prediction band",
        "1:1 reference"),
      col = c(grDevices::adjustcolor(colour[1L], alpha.f = 0.5),
        "black", area_band_colour, reference_colour),
      pch = c(16, NA, 15, NA), lty = c(NA, 1, NA, 2), lwd = c(NA, 2.3, NA, 1.65),
      pt.cex = c(1.30, 1, 1.45, 1))
    graphics::mtext(expression(log[10]*"(True peak area, integrated counts)"),
      1, outer = TRUE, line = 0.6, cex = 1.1)
  })


  benchmark16Similarity(hits,"spectral_cosine",methods,labels,colour,fill,draw,
    "simulated-16-spectral-violin.png","Cosine similarity (logit scale)",
    "Reciprocal spectral hits across 16 samples",stats::median,"Median",
    c(.5,.7,.9,.99,.999,.9999,.99999,.999999),nrow(truth),.7,cap=.999999)
  benchmark16Similarity(hits[hits$present,,drop=FALSE],"spectral_F1",methods,
    labels,colour,fill,draw,"simulated-16-spectral-support.png",
    "Within-spectrum F1 (logit scale)","Correctness of spectral ion support",
    mean,"Mean",c(0,.25,.5,.75,.9,.95,.99,1),nrow(truth))
  benchmark16Heatmap(results,draw)
  benchmark16ProfilePlots(results,methods,labels,draw)
  benchmark16AreaHeatmaps(results,output_dir,methods,labels,draw)
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
  stats$label<-sprintf("%s %.6f · %d/%d true",stat_label,stats$value,stats$true,total)
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

benchmark16Palette <- function(n=1000L)
  grDevices::colorRampPalette(c("black","blue","green","yellow","orange","red"))(n)

benchmark16Heatmap <- function(results,draw) {
  file<-file.path(results,"sample-heatmap.rds")
  if(!file.exists(file))return(invisible(NULL))
  x<-readRDS(file)
  grid<-seq(min(x$time),max(x$time),length.out=length(x$time))
  matrices<-lapply(list(x$observed,x$mean),function(z)
    apply(z,2L,function(v)stats::approx(x$time,v,grid)$y))
  lim<-range(c(0,sqrt(unlist(matrices))))
  draw("simulated-16-sample-heatmap.png",2200,1250,{
    graphics::layout(matrix(c(1,3,2,3),2,2,byrow=TRUE),widths=c(12,1.6))
    graphics::par(mar=c(3.3,4.2,2.5,1),mgp=c(2.3,.7,0),las=1)
    titles<-c("Observed counts","True mean, including smooth background")
    for(i in 1:2) {
      graphics::image(grid,x$mass,sqrt(matrices[[i]]),col=benchmark16Palette(),zlim=lim,
        useRaster=TRUE,xlab="Retention time (min)",ylab="Mass-to-charge ratio",
        main=paste(x$sample,titles[i],sep=" · "))
      graphics::abline(v=x$apex[x$present],col=grDevices::adjustcolor("red",.4),lwd=.65)
      graphics::abline(v=x$apex[intersect(which(x$present),x$acid)],col="red",lwd=1.8)
    }
    graphics::par(mar=c(3.3,.2,2.5,4.7))
    values<-seq(lim[1],lim[2],length.out=1000L)
    graphics::image(c(0,1),values,matrix(rep(values,each=2),2),axes=FALSE,
      col=benchmark16Palette(),xlab="",ylab="",useRaster=TRUE)
    at<-pretty(lim,4);at<-at[at>=lim[1]&at<=lim[2]]
    graphics::axis(4,at=at,labels=format(at^2,trim=TRUE,scientific=TRUE),cex.axis=.8)
    graphics::mtext("Intensity (square-root colours)",4,line=3.6,cex=.8,las=0)
  })
}

benchmark16ProfilePlots <- function(results,methods,labels,draw) {
  file<-file.path(results,"profiles.rds")
  if(!file.exists(file))return(invisible(NULL))
  x<-readRDS(file)
  quality<-utils::read.csv(file.path(results,"profile-quality.csv"))
  labels<-c("Ground truth"="Ground truth",labels)
  available<-c("Ground truth",intersect(methods,unique(x$method)))
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
      graphics::layout(matrix(seq_along(available),2L,3L,byrow=TRUE))
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
        if(method!="Ground truth")curves(d,1,1.5)
        curves(known,2,1.2)
        q<-quality[quality$kind==kind &quality$method==method,,drop=FALSE]
        if(nrow(q))graphics::mtext(sprintf("Profile cosine %.3f · %d/%d true calls",
          q$mean_profile_cosine,q$recovered,q$present),3,line=.2,cex=.75)
      }
      graphics::mtext(if(kind=="Regular")"Highest-abundance regular compound"else
        "Shark-fin-shaped acid compound",3,outer=TRUE,line=.9,font=2,cex=1.3)
      graphics::mtext("Solid: inferred profiles · Dashed: ground truth · All 16 samples, blue (smallest area) to red (largest area)",
        1,outer=TRUE,line=1.1,cex=.95)
      graphics::mtext("Original total-ion amplitudes retained. Zero-area samples are included; no unit-area or peak-height normalization.",
        1,outer=TRUE,line=.05,cex=.9)
    })
  }
}

benchmark16AreaHeatmaps <- function(results,output_dir,methods,labels,draw) {
  file<-file.path(results,"area-grid.csv")
  if(!file.exists(file))return(invisible(NULL))
  x<-utils::read.csv(file)
  order<-rev(c("Ground truth",methods))
  x<-x[x$method%in%order,,drop=FALSE]
  labels<-c("Ground truth"="Ground truth",labels)
  k<-max(x$truth_id)
  lim<-range(log10(x$area_estimated[x$area_estimated>0]),finite=TRUE)
  palette<-grDevices::colorRampPalette(c("blue2","red2"))(1000L)
  panel<-function(z,title) {
    graphics::layout(matrix(1:2,1,2),widths=c(15,1.7))
    graphics::par(mar=c(4.5,8,3,1),mgp=c(2.7,.7,0),las=1)
    values<-matrix(NA_real_,k,length(order));fn<-matrix(FALSE,k,length(order))
    i<-cbind(z$truth_id,match(z$method,order))
    values[i]<-ifelse(z$area_estimated>0,log10(z$area_estimated),NA_real_)
    fn[i]<-z$false_negative
    graphics::plot.new();graphics::plot.window(c(.5,k+.5),c(.5,length(order)+.5),xaxs="i",yaxs="i")
    graphics::rect(.5,.5,k+.5,length(order)+.5,col="#edf0f2",border=NA)
    graphics::image(seq_len(k),seq_along(order),values,add=TRUE,col=palette,
      zlim=lim,useRaster=TRUE)
    j<-which(fn,arr.ind=TRUE)
    if(nrow(j))graphics::rect(j[,1]-.5,j[,2]-.5,j[,1]+.5,j[,2]+.5,col="yellow",border=NA)
    graphics::axis(1,at=unique(c(1,seq(10,k,10))))
    graphics::axis(2,at=seq_along(order),labels=unname(labels[order]),cex.axis=.9)
    graphics::title(main=title,xlab="Compound (ordered by retention index)")
    graphics::box()
    graphics::mtext("Yellow: false negative · Grey: no signal/call · Positive cells show total-ion count equivalents",1,line=3.6,cex=.83)
    graphics::par(mar=c(4.5,.5,3,4.5))
    levels<-seq(lim[1],lim[2],length.out=1000)
    graphics::image(c(0,1),levels,matrix(rep(levels,each=2),2),col=palette,
      axes=FALSE,xlab="",ylab="",useRaster=TRUE)
    graphics::axis(4,at=pretty(lim,4))
    graphics::mtext(expression(log[10]*"(peak area)"),4,line=2.8,las=0)
  }
  selected<-readRDS(file.path(results,"sample-heatmap.rds"))$sample
  draw("simulated-16-area-heatmap.png",2300,700,
    panel(x[x$sample==selected,,drop=FALSE],paste(selected,"· true and inferred peak areas")))
  pooled<-do.call(rbind,lapply(split(x,interaction(x$method,x$truth_id,drop=TRUE)),function(z)
    data.frame(method=z$method[1],truth_id=z$truth_id[1],area_estimated=sum(z$area_estimated),
      false_negative=any(z$present)&!any(z$called &z$present))))
  draw("simulated-16-area-heatmap-pooled.png",2300,700,
    panel(pooled,"True and inferred peak areas pooled across 16 samples"))
  grDevices::pdf(file.path(output_dir,"simulated-16-area-heatmaps-all-samples.pdf"),width=15,height=4.6)
  on.exit(grDevices::dev.off())
  for(sample in unique(x$sample))panel(x[x$sample==sample,,drop=FALSE],
    paste(sample,"· true and inferred peak areas"))
}
