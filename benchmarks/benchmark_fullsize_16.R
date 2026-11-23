# Reproducible full-size resource/quality benchmark, separate from the
# 100-compound seven-engine comparison. Run from the repository root.
# Requires L0glm, logconNMF, GCMSprime, Matrix, and the adjacent GCMSbenchmark.
# CDF simulation and PSF-standard calibration are timed separately from fitting.
# Source this file, then call benchmarkFullsize16(stage="all").

benchmarkFullsize16 <- function(directory="codex_sandbox/fullsize16",
  stage=c("all","prepare","fit","score"),ncores=NULL,
  background_file="GCMSbenchmark/inst/benchmarks/gcms-data2-background.rds") {
  stage <- match.arg(stage)
  library(L0glm);library(logconNMF);library(GCMSprime);library(Matrix)
  root <- normalizePath(directory,winslash="/",mustWork=FALSE)
  input <- file.path(root,"input")
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  if(stage %in% c("all","prepare")) {
    if(file.exists(file.path(root,"truth.rds")))
      stop("Use a new directory for simulation, or stage='fit'/'score'.")
    dir.create(file.path(input,'samples'),recursive=TRUE,showWarnings=FALSE)
    dir.create(file.path(input,'standards'),showWarnings=FALSE)
    data(simulatedSpectra,package='L0glm')
    mass <- 30:579
    S <- matrix(0,nrow(simulatedSpectra),length(mass),dimnames=list(rownames(simulatedSpectra),as.character(mass)))
    S[,match(as.integer(colnames(simulatedSpectra)),mass)] <- as.matrix(simulatedSpectra)
    template <- readRDS(background_file)
    oldbg <- template$basis%*%template$alpha
    background <- matrix(.2*stats::median(oldbg),nrow(oldbg),length(mass))
    background[,match(as.integer(colnames(simulatedSpectra)),mass)] <- oldbg
    cat('16 samples; 10,000 scans; 550 channels; 1,000 reference compounds, 50% presence; one acid.\n')
    cat('Dispersion',template$dispersion,'typical area',template$area,'\n');flush.console()
    simulation <- simulateGCMSSamples(samples=16L,n=10000L,k=1000L,S=S,
      preset='time_area_acid',acid_fraction=1/1000,presence=.5,irregular=.02,
      drift=.5,rt_area=0,family='quasipoisson',dispersion=template$dispersion,
      area=template$area,background=background,seed=171117L)
    truth <- vector('list',16L)
    for(s in 1:16) {
      x <- simulation$samples[[s]];x$time <- x$time*.3
      file <- file.path(input,'samples',paste0(x$sample_id,'.cdf'))
      writeGCMSCDF(x$Y,file,time=x$time,acquired_at=x$acquired_at)
      present <- which(x$present)
      H <- x$H_true[,present,drop=FALSE]*x$exposure
      maximum <- apply(H,2,max)
      H[H<rep(maximum*1e-12,each=nrow(H))] <- 0
      H <- as(Matrix(H,sparse=TRUE),'dgCMatrix');H <- H%*%Diagonal(x=1/colSums(H))
      truth[[s]] <- list(sample=x$sample_id,file=file,time=x$time,id=present,H=H,
        S=x$P_true[present,,drop=FALSE],area=x$area_true[present],acid=present%in%x$acid,
        mass=mass,dispersion=x$dispersion,apex=x$time[max.col(t(as.matrix(H)),ties.method='first')])
      cat('Exported',s,'compounds',length(present),'acid present',any(truth[[s]]$acid),'\n');flush.console()
      simulation$samples[s] <- list(NULL)
    }
    saveRDS(truth,file.path(root,'truth.rds'))
    saveRDS(list(P=simulation$P,compound_ids=simulation$compound_ids,present=simulation$present,
      areas=simulation$areas,seed=simulation$seed,settings=simulation$settings,
      description='1000 library spectra in 550 recorded mass bins (258 bins carry simulated compound ions), approximately 500 peaks per sample; one acid; concentration-dependent EMG broadening, drift and a real-window-derived smooth background'),
      file.path(root,'settings.rds'))
    rm(simulation,truth);gc()
    standard <- simulateGCMSSamples(samples=1L,n=10000L,k=20L,S=S,preset='time',
      acid_fraction=0,presence=1,irregular=.02,drift=0,rt_area=0,
      family='quasipoisson',dispersion=template$dispersion,area=template$area,
      background=background,seed=171118L)$samples[[1L]]
    writeGCMSCDF(standard$Y,file.path(input,'standards','psf-standard.cdf'),time=standard$time*.3,
      acquired_at=as.POSIXct('2026-11-01 07:00:00',tz='UTC'))
    cat('Full-size CDF benchmark inputs complete.\n')
  }
  if(stage %in% c("all","fit")) {
    standard_file <- file.path(root,"standard-calibration.rds")
    if(!file.exists(standard_file)) {
      raw <- getFromNamespace("GCMS_read_cdf","GCMSprime")(
        file.path(input,"standards","psf-standard.cdf"),NULL)
      standard <- estimatePSF(raw$Y,time=raw$time,
        family=quasipoisson("identity"),method="standard",nwindows=10L,
        precluster=1L,ncores=8L)
      saveRDS(standard,standard_file)
      rm(raw);gc()
    } else standard <- readRDS(standard_file)
    # Omit ncores or use NULL to test the automatic logical-core budget.
    object <- gcmsDeco(file.path(input,"samples"),
      output_dir=file.path(root,"analysis"),ncores=ncores,
      psf_control=list(standard=standard),verbose=TRUE)
    saveRDS(object,file.path(root,"analysis.rds"))
    capture.output(sessionInfo(),file=file.path(root,"sessionInfo.txt"))
    print(object$parallel);print(object$seconds)
  }
  if(stage %in% c("all","score")) {
    metrics <- new.env(parent=globalenv())
    source("GCMSbenchmark/R/metrics.R",local=metrics)
    object <- readRDS(file.path(root,"analysis.rds"))
    truth <- readRDS(file.path(root,"truth.rds"))
    names(truth) <- vapply(truth,`[[`,"","sample")
    assessed <- lapply(object$paths,function(path) {
      o <- readRDS(path);t <- truth[[o$sample_name]]
      if(is.null(t))t <- truth[[tools::file_path_sans_ext(basename(o$cdf_file))]]
      stopifnot(!is.null(t));f <- o$model
      H <- as.matrix(f$H);S <- as.matrix(f$S);area <- rowSums(S)
      record <- list(peaks=data.frame(id=seq_len(ncol(H)),
        rt=o$time[max.col(t(H),ties.method="first")],area=area),
        spectra=lapply(seq_len(ncol(H)),function(i) {
          keep <- S[i,]>=.001*max(S[i,]) & S[i,]>0
          cbind(o$mass[keep],S[i,keep])
        }),profiles=lapply(seq_len(ncol(H)),function(i)cbind(o$time,H[,i]*area[i])))
      t$P <- t$S;t$H <- as.matrix(t$H)%*%Diagonal(x=t$area)
      pairs <- metrics$match_sample(record,seq_len(ncol(H)),t,o$sample_name,
        "GCMSprime",rt_tolerance=.25,min_cosine=.9,ion_threshold=.001)
      TP <- nrow(pairs);FP <- ncol(H)-TP;FN <- length(t$id)-TP
      list(summary=data.frame(sample=o$sample_name,called=ncol(H),
        true_compounds=length(t$id),TP=TP,FP=FP,FN=FN,
        F1=2*TP/(2*TP+FP+FN),spectral_cosine=median(pairs$spectral_cosine),
        spectral_F1=mean(pairs$spectral_F1),profile_cosine=mean(pairs$profile_cosine),
        median_area_ratio=median(pairs$area_estimated/pairs$area_true),
        acid_called=sum(pairs$acid),acid_present=sum(t$acid)),matches=pairs)
    })
    summary <- do.call(rbind,lapply(assessed,`[[`,"summary"))
    pairs <- do.call(rbind,lapply(assessed,`[[`,"matches"))
    write.csv(summary,file.path(root,"quality.csv"),row.names=FALSE)
    write.csv(pairs,file.path(root,"matches.csv"),row.names=FALSE)
    print(summary);cat("Joint GIC:",GIC(object),"\n")
  }
  invisible(root)
}

# Recreate the large-benchmark positive-case acid display without refitting.
# All matched acid-positive samples are shown; acid-negative samples are omitted.
plotFullsize16Acid <- function(directory="codex_sandbox/fullsize16", file=NULL) {
  object <- readRDS(file.path(directory,"analysis.rds"))
  truth <- readRDS(file.path(directory,"truth.rds"))
  names(truth) <- vapply(truth,`[[`,"","sample")
  pairs <- utils::read.csv(file.path(directory,"matches.csv"))
  acid <- pairs[pairs$acid,,drop=FALSE]
  acid <- acid[order(acid$area_true),,drop=FALSE]
  if(!nrow(acid)) stop("No matched acid occurrences to plot.")
  positive <- sum(vapply(truth,function(x)sum(x$acid),0L))
  curves <- lapply(seq_len(nrow(acid)),function(i) {
    row <- acid[i,]
    t <- truth[[sub("^samples[.]","",row$sample)]]
    o <- readRDS(object$paths[match(row$sample,object$samples$sample_id)])
    k <- match(as.integer(row$truth_id),t$id)
    h <- as.numeric(t$H[,k])*t$area[k]
    inferred <- as.numeric(o$model$H[,row$estimate_id])*sum(o$model$S[row$estimate_id,])
    onset <- t$time[which(h>=max(h)*.01)[1L]]
    keep <- which(pmax(h,inferred)>max(h,inferred)*1e-5)
    list(time=(t$time-onset)*60,true=h,estimate=inferred,keep=keep)
  })
  if(is.null(file)) file <- file.path(directory,"fullsize-16-acid.png")
  grDevices::png(file,width=1900,height=1000,res=180)
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::par(mar=c(5,5,4,1))
  xr <- range(unlist(lapply(curves,function(z)z$time[z$keep])))
  yr <- max(unlist(lapply(curves,function(z)c(z$true,z$estimate))))
  graphics::plot(xr,c(0,yr*1.05),type="n",
    xlab="Seconds relative to the true 1%-height onset",
    ylab="Total-ion intensity per scan",
    main=sprintf("Large-scale %d-sample benchmark | acid-positive profiles",length(truth)))
  graphics::mtext(sprintf("10,000 scans; 1,000 reference compounds | %d matched / %d positive samples; %d acid-negative samples omitted",
    nrow(acid),positive,length(truth)-positive),3,line=.6,cex=.8)
  colours <- grDevices::colorRampPalette(c("blue2","red2"))(length(curves))
  for(i in seq_along(curves)) {
    z <- curves[[i]];at <- z$keep
    graphics::lines(z$time[at],z$true[at],lty=2,col=colours[i],lwd=1)
    graphics::lines(z$time[at],z$estimate[at],col=colours[i],lwd=1.5)
  }
  graphics::legend("topleft",c("Inferred","Ground truth"),
    lty=c(1,2),lwd=c(1.5,1),bty="n")
  graphics::mtext("Original amplitudes | blue to red: increasing true area among the displayed positive samples",
    1,line=3.8,cex=.8)
  invisible(file)
}


# Command-line stages use the same function as the interactive recipe.
if(sys.nframe()==0L && !interactive()) {
  args <- commandArgs(TRUE)
  benchmarkFullsize16(stage=if(length(args))args[1L]else"all",
    directory=if(length(args)>1L)args[2L]else"codex_sandbox/fullsize16")
}
