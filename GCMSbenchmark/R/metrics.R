empty_matches <- function() data.frame(method=character(),sample=character(),
  truth_id=character(), estimate_id=character(), peak_row=integer(),
  rt_true=numeric(),rt_estimated=numeric(),rt_error_min=numeric(),acid=logical(),
  spectral_cosine=numeric(),profile_cosine=numeric(),profile_width_ratio=numeric(),
  area_true=numeric(),area_estimated=numeric())

cosine <- function(a,b) {
  z <- sqrt(sum(a*a)*sum(b*b))
  if (z>0) min(1,sum(a*b)/z) else NA_real_
}

spectrum_matrix <- function(spectra, mass) {
  m <- matrix(0,length(spectra),length(mass))
  for (j in seq_along(spectra)) {
    z <- spectra[[j]]
    if (!length(z)) next
    k <- match(round(z[,1L]),mass)
    valid <- is.finite(z[,2L]) & z[,2L]>0 & !is.na(k)
    if (!any(valid)) next
    s <- rowsum(z[valid,2L],k[valid],reorder=FALSE)
    m[j,as.integer(rownames(s))] <- s[,1L]
  }
  m
}

norm_rows <- function(x) {
  n <- sqrt(rowSums(x*x))
  x/pmax(n,.Machine$double.xmin)
}

profile_quality <- function(profile, time, truth) {
  if (is.null(profile) || nrow(profile)<2L) return(c(NA_real_,NA_real_))
  profile <- profile[is.finite(profile[,1L]) & is.finite(profile[,2L]),,drop=FALSE]
  if (length(unique(profile[,1L]))<2L) return(c(NA_real_,NA_real_))
  fit <- stats::approx(profile[,1L],profile[,2L],time,yleft=0,yright=0,ties=mean)$y
  width <- function(z) {
    i <- which(z>=.5*max(z) & z>0)
    if (length(i)>1L) diff(range(time[i])) else 0
  }
  w <- width(truth)
  c(cosine(fit,truth),if(w>0)width(fit)/w else NA_real_)
}

match_sample <- function(fit, ii, truth, sample, method, rt_tolerance, min_cosine) {
  if (!length(ii) || !length(truth$id)) return(empty_matches())
  ids <- fit$peaks[ii,,drop=FALSE]
  masses <- sort(unique(round(c(truth$mass,
    unlist(lapply(fit$spectra[ii],function(z) if(length(z)) z[,1L]))))))
  specs <- spectrum_matrix(fit$spectra[ii],masses)
  true_spec <- matrix(0,nrow(truth$P),length(masses))
  true_spec[,match(round(truth$mass),masses)] <- truth$P
  sim <- tcrossprod(norm_rows(true_spec),norm_rows(specs))
  comparison_rt <- if (identical(fit$rt_basis,"aligned") &&
    !is.null(truth$reference_apex)) truth$reference_apex else truth$apex
  rt_gap <- abs(outer(comparison_rt,ids$rt,"-"))
  eligible <- is.finite(rt_gap) & rt_gap<=rt_tolerance & sim>=min_cosine
  # Maximise cardinality first. The sum of all similarity bonuses is
  # smaller than the reward for one additional admissible pair.
  k <- length(truth$id)
  reward <- matrix(0,k,nrow(ids)+k)
  reward[,seq_len(nrow(ids))][eligible] <- k+1+sim[eligible]+
    1e-6*(1-rt_gap[eligible]/rt_tolerance)
  j <- as.integer(clue::solve_LSAP(reward,maximum=TRUE))
  pair <- cbind(seq_len(k),j)
  i <- which(j<=nrow(ids) & reward[pair]>0)
  j <- j[i]
  if (!length(i)) return(empty_matches())
  prof <- vapply(seq_along(i),function(l)
    profile_quality(fit$profiles[[ii[j[l]]]],truth$time,truth$H[,i[l]]),c(0,0))
  data.frame(method=method,sample=sample,truth_id=as.character(truth$id[i]),
    estimate_id=ids$id[j],peak_row=ii[j],rt_true=comparison_rt[i],
    rt_estimated=ids$rt[j],rt_error_min=ids$rt[j]-comparison_rt[i],
    acid=as.logical(truth$acid[i]),spectral_cosine=sim[cbind(i,j)],
    profile_cosine=prof[1L,],profile_width_ratio=prof[2L,],
    area_true=truth$area[i],area_estimated=ids$area[j])
}

metric_median <- function(x) {
  if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
}

area_regression <- function(pairs) {
  true <- pairs$area_true
  estimated <- pairs$area_estimated
  keep <- is.finite(true) & true>0 & is.finite(estimated) & estimated>0
  data <- data.frame(log_true=log10(true[keep]),log_estimated=log10(estimated[keep]))
  fit <- NULL
  r2 <- NA_real_
  if (nrow(data)>=3L && diff(range(data$log_true))>0) {
    fit <- stats::lm(log_estimated ~ log_true,data=data)
    if (diff(range(data$log_estimated))>0)
      r2 <- suppressWarnings(summary(fit)$adj.r.squared)
  }
  list(data=data,fit=fit,adjusted_R2=r2)
}

score_row <- function(method,fit,truth,pairs) {
  ok <- identical(fit$status,"ok")
  optional <- function(x, default=NA_real_) if (length(x)==1L) x else default
  n <- if(ok) nrow(fit$peaks) else NA_integer_
  ntruth <- if (is.null(truth)) NA_integer_ else sum(lengths(lapply(truth,`[[`,"id")))
  scored <- ok && !is.na(ntruth)
  tp <- if(scored) nrow(pairs) else NA_integer_
  fp <- n-tp
  fn <- ntruth-tp
  f1 <- if(scored && (n+ntruth)>0) 2*tp/(n+ntruth) else NA_real_
  acid_true <- if(is.null(truth)) NA_integer_ else
    sum(vapply(truth,function(z) sum(z$acid),0))
  acid <- pairs[pairs$acid,,drop=FALSE]
  area_model <- area_regression(pairs)
  ar <- log(pairs$area_estimated/pairs$area_true)
  ar[!is.finite(ar)] <- NA_real_
  ratio <- if (identical(fit$area_kind,"integrated counts"))
    exp(metric_median(abs(ar)))-1 else NA_real_
  ac_ratio <- if (identical(fit$area_kind,"integrated counts"))
    exp(metric_median(abs(ar[pairs$acid])))-1 else NA_real_
  across <- split(ar,pairs$truth_id)
  fold <- unlist(lapply(across,function(x) {
    x <- x[is.finite(x)]
    if(length(x)>1L) abs(x-stats::median(x)) else numeric()
  }),use.names=FALSE)
  data.frame(method=method,status=fit$status,message=fit$message,
    note=optional(fit$note,""),analysis_s=optional(fit$analysis),
    input_conversion_s=optional(fit$preparation),adapter_elapsed_s=optional(fit$total),
    components=n,shared=optional(fit$shared),TP=tp,FP=fp,FN=fn,F1=f1,
    TP_cosine_90=if(scored)sum(pairs$spectral_cosine>=.90) else NA_integer_,
    median_abs_rt_error_s=60*metric_median(abs(pairs$rt_error_min)),
    spectral_cosine=metric_median(pairs$spectral_cosine),
    profile_cosine=metric_median(pairs$profile_cosine),
    profiles_scored=sum(is.finite(pairs$profile_cosine)),
    native_area=optional(fit$area_kind,NA_character_),
    native_profile=optional(fit$profile_kind,NA_character_),
    area_pairs=nrow(area_model$data),area_adjusted_R2=area_model$adjusted_R2,
    area_median_error_factor=ratio,
    area_median_scaled_error_factor=exp(metric_median(abs(ar-metric_median(ar))))-1,
    within_compound_area_error_factor=exp(metric_median(fold))-1,
    acid_truth=acid_true,acid_TP=if(scored) nrow(acid) else NA_integer_,
    acid_TP_cosine_90=if(scored)sum(acid$spectral_cosine>=.90) else NA_integer_,
    acid_spectral_cosine=metric_median(acid$spectral_cosine),
    acid_profile_cosine=metric_median(acid$profile_cosine),
    acid_profile_width_ratio=metric_median(acid$profile_width_ratio),
    acid_area_error_factor=ac_ratio)
}

scoreGCMSbenchmark <- function(x,truth=NULL,rt_tolerance=.12,min_cosine=.70) {
  if (!inherits(x,"GCMSbenchmark")) stop("x must be a GCMSbenchmark result.",call.=FALSE)
  if (length(rt_tolerance)!=1L || !is.finite(rt_tolerance) || rt_tolerance<=0 ||
      length(min_cosine)!=1L || !is.finite(min_cosine) || min_cosine<=0 || min_cosine>1)
    stop("Use a positive RT tolerance and 0 < min_cosine <= 1.",call.=FALSE)
  if (!is.null(truth)) {
    fields <- c("id","time","H","P","area","apex","acid","mass")
    ids <- sample_id(x$files)
    if (!is.list(truth) || !setequal(names(truth),ids) || anyDuplicated(names(truth)))
      stop("Supply one named truth record per input sample.",call.=FALSE)
    for (z in truth) {
      k <- length(z$id)
      if (!all(fields%in%names(z)) || anyDuplicated(z$id) ||
          anyDuplicated(round(z$mass)) || !all(is.finite(z$mass)) ||
          !identical(dim(z$H),c(length(z$time),k)) ||
          !identical(dim(z$P),c(k,length(z$mass))) ||
          any(!is.finite(z$time)) || is.unsorted(z$time,strictly=TRUE) ||
          !all(lengths(z[c("area","apex","acid")])==k) ||
          anyNA(z$acid) || any(!is.finite(z$apex)) ||
          (!is.null(z$reference_apex) && (length(z$reference_apex)!=k ||
             any(!is.finite(z$reference_apex)))))
        stop("Malformed benchmark truth record.",call.=FALSE)
    }
  }
  x$matches <- empty_matches()
  rows <- lapply(names(x$methods),function(method) {
    fit <- x$methods[[method]]
    pairs <- if (is.null(truth) || !identical(fit$status,"ok")) empty_matches() else
      do.call(rbind,lapply(names(truth),function(id) match_sample(fit,
        which(fit$peaks$sample==id),truth[[id]],id,method,rt_tolerance,min_cosine)))
    x$matches <<- rbind(x$matches,pairs)
    score_row(method,fit,truth,pairs)
  })
  x$summary <- do.call(rbind,rows)
  x
}
