# Use the published occurrence counts and matched areas without refitting.
plotFullsize16Recovery <- function(directory="benchmarks/results/fullsize16",
                                   file="fullsize-16-recovery.png") {
  metrics <- utils::read.csv(file.path(directory,"metrics.csv"))
  pairs <- utils::read.csv(file.path(directory,"matches.csv"))
  rank_file <- file.path(directory,"rank.csv")
  fit <- utils::getFromNamespace("area_regression","GCMSbenchmark")(pairs)
  if (is.null(fit$fit)) stop("At least three positive area pairs and two distinct true areas are required.")
  x <- fit$data$log_true
  y <- fit$data$log_estimated
  grid <- seq(min(x),max(x),length.out=200L)
  band <- stats::predict(fit$fit,data.frame(log_true=grid),interval="prediction",level=0.95)
  acid <- with(pairs,acid[is.finite(area_true) & area_true>0 & is.finite(area_estimated) & area_estimated>0])
  grDevices::png(file,width=2200L,height=1000L,res=180)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow=c(1,2),mar=c(5,5,4,1))
  heights <- rbind(metrics$TP,metrics$FP)
  mid <- graphics::barplot(heights,col=c("#1876C2","#DF4747"),border=NA,
    names.arg=rep("",nrow(metrics)),xlab="Sample",ylab="Compound occurrences",
    ylim=c(0,1.23*max(colSums(heights),metrics$true_compounds)),
    main=sprintf("%s / %s true occurrences recovered",format(sum(metrics$TP),big.mark=","),
                 format(sum(metrics$true_compounds),big.mark=",")))
  graphics::points(mid,metrics$true_compounds,pch=1,cex=.85)
  if (file.exists(rank_file)) {
    ranks <- utils::read.csv(rank_file)
    graphics::mtext(sprintf("True rank %s   |   Estimated rank %s%s",
      format(ranks$true_rank,big.mark=","),format(ranks$estimated_rank,big.mark=","),
      if (isTRUE(ranks$rank_exact)) " (exact)" else ""),3,line=.35,cex=.83)
  }
  graphics::axis(1,at=mid,labels=seq_len(nrow(metrics)),cex.axis=.8,tick=FALSE)
  graphics::legend("topleft",c("True positives","False positives","True count"),
    fill=c("#1876C2","#DF4747",NA),pch=c(NA,NA,1),border=NA,bty="n",cex=.85,ncol=3L)
  lim <- range(x,y,band)
  graphics::plot(x,y,type="n",xlim=lim,ylim=lim,xlab=expression(log[10]~"true total-ion area"),
    ylab=expression(log[10]~"estimated total-ion area"),main="Peak-area accuracy")
  graphics::mtext(sprintf("Adjusted R² %.6f   |   Median absolute area error %.3f%%",
    fit$adjusted_R2,100*stats::median(abs(10^(y-x)-1))),3,line=.35,cex=.83)
  graphics::polygon(c(grid,rev(grid)),c(band[,"lwr"],rev(band[,"upr"])),col="grey83",border=NA)
  graphics::lines(grid,band[,"fit"],lwd=2)
  graphics::abline(0,1,lty=2,lwd=1.2,col="#23364D")
  blue <- grDevices::adjustcolor("#1876C2",alpha.f=.5)
  violet <- grDevices::adjustcolor("#69418e",alpha.f=.5)
  graphics::points(x[!acid],y[!acid],pch=16,cex=.37,col=blue)
  graphics::points(x[acid],y[acid],pch=23,cex=.95,col="#563470",bg=violet)
  graphics::legend("topleft",c("Matched peak","Matched acid peak","Least-squares fit",
                              "95% prediction interval","1:1 reference"),
    col=c(blue,"#563470","black","grey83","#23364D"),pt.bg=c(NA,violet,NA,NA,NA),
    lty=c(NA,NA,1,NA,2),lwd=c(NA,NA,2,NA,1.2),pch=c(16,23,NA,15,NA),bty="n",cex=.77)
  invisible(file)
}
