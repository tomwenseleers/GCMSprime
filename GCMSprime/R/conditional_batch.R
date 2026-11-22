# Rebin cumulative observed counts to calibrate one PSF on a common clock.
# The factorization itself always uses the measured counts and times.
GCMS_scan_edges <- function(t) {
  d <- diff(t)
  c(t[1L]-d[1L]/2, utils::head(t,-1L)+d/2, utils::tail(t,1L)+utils::tail(d,1L)/2)
}

GCMS_background <- function(time, width, df) {
  n <- length(time)
  if (is.null(df)) df <- max(4L, min(40L, 3L + ceiling(n/(50*width))))
  if (length(df)!=1L || !is.finite(df) || df<1L || df!=floor(df) || df>n || df %in% 2:3)
    stop("Use baseline_df=1 for a constant or an integer from four to the number of scans.")
  if (df==1L) matrix(1,n,1L) else as.matrix(splines::bs(time,df=df,intercept=TRUE))
}

GCMS_pooled_starts <- function(manifest, args, output, precluster, ncores) {
  slices <- lapply(seq_len(nrow(manifest)), function(i) {
    id <- manifest$sample_id[i]
    d <- file.path(output, "samples", id)
    dir.create(d, showWarnings = FALSE)
    raw <- GCMS_read_cdf(manifest$cdf_file[i], args$time_range)
    noise <- estimateSignalNoise(raw$Y, args$family, args$dispersion, args$theta,
      correlation=args$correlation, ncores=ncores)
    observed <- file.path(d, "observed.rds")
    saveRDS(raw$Y, observed, compress=FALSE)
    list(y=observed, time=raw$time, scans=raw$scans, mass=raw$mass, total=rowSums(raw$Y),
      spec=noise$spec, working_weight=noise$working_weight, d=d)
  })
  mass <- slices[[1L]]$mass
  if (any(!vapply(slices, function(x) identical(x$mass, mass), TRUE))) return(NULL)
  edges <- lapply(slices, function(x) GCMS_scan_edges(x$time))
  n <- as.integer(round(stats::median(lengths(lapply(slices, "[[", "time")))))
  if (n<9L || args$nwindows>n) stop("Use at least nine scans and at most one PSF window per scan.")
  bins <- seq(min(vapply(edges,min,0)), max(vapply(edges,max,0)), length.out=n+1L)
  time <- utils::head(bins,-1L)+diff(bins)/2
  totals <- vapply(seq_along(slices), function(i) {
    diff(stats::approx(edges[[i]], c(0,cumsum(slices[[i]]$total)), bins, rule=2)$y)
  }, numeric(n))
  weight <- colSums(totals)
  phi <- vapply(slices,function(x) x$spec$dispersion,0)
  ctrl <- utils::modifyList(list(fit=list(IC=if(is.character(args$lambda)) args$lambda else "gic")),
    args$psf_control[setdiff(names(args$psf_control),"standard")])
  standard <- args$psf_control$standard
  standard_file <- manifest$alkanes_file[1L]
  if (is.null(standard) && !is.na(standard_file))
    standard <- GCMS_read_cdf(standard_file,args$calibration_range)[c("Y","time")]
  common <- estimatePSF(matrix(rowSums(totals),n,1L), time=time, family=stats::quasipoisson("identity"),
    method="auto", standard=standard, dispersion=if(sum(weight)>0) stats::weighted.mean(phi,weight) else mean(phi),
    precluster=FALSE, nwindows=args$nwindows, ncores=ncores, control=ctrl)
  k <- if (identical(precluster, TRUE)) 32L else precluster
  groups <- if (inherits(k,"nmf_channel_groups")) k else NULL
  if (is.numeric(k) && (length(k)!=1L || !is.finite(k) || k!=floor(k) || k<1L))
    stop("precluster must be a positive integer or a channel map.")
  if (is.numeric(k) && k < length(mass)) {
    width <- max(1L, round(common$width_scans/2))
    summary <- do.call(rbind,lapply(slices,function(x) {
      y <- readRDS(x$y)
      rowsum(y,(seq_len(nrow(y))-1L)%/%width)
    }))
    groups <- logconNMF::channelGroups(summary,k=k,channels=as.character(mass),ncores=ncores,
      reserve_unknown=FALSE)
  }
  results <- lapply(seq_along(slices),function(i) {
    z <- slices[[i]]
    psf <- estimatePSF(matrix(z$total,ncol=1L),time=z$time,standard=common,method="transfer",ncores=ncores)
    width <- max(1,psf$width_scans)
    B <- GCMS_background(z$time,width,args$baseline_df)
    id <- manifest$sample_id[i]
    template <- structure(list(psf=psf,psf_standard=psf$calibration,psf_width_scans=width,background_df=ncol(B),
      alkanes_file=if(!is.na(standard_file)) standard_file else NULL,
      psf_working_weight=z$working_weight,time=z$time,scans=z$scans,mass=mass,
      cdf_file=manifest$cdf_file[i],sample_name=id,ion_threshold=args$ion_threshold,
      seconds=c(total=0)),class="GCMSdeco")
    path <- file.path(z$d,"fit.rds")
    saveRDS(template,path)
    desc <- list(y=z$y, time=z$time, X=psf$X, background=B, dispersion=z$spec$dispersion,
      noise_spec=z$spec, sample_id=id, group=manifest$group[i],
      small_width=width*sqrt(log(100)/log(2))*stats::median(diff(z$time)))
    list(path=path,descriptor=desc,mass=mass,cache_hit=FALSE)
  })
  list(results=results, groups=groups)
}
