library(GCMSprime)
library(logconNMF)

benchmark_samples16 <- function(output = tempfile("GCMSprime_samples16_"),
                                cores = c(1L, 8L, 16L)) {
  blas <- RhpcBLASctl::blas_get_num_procs()
  on.exit(RhpcBLASctl::blas_set_num_threads(blas))
  RhpcBLASctl::blas_set_num_threads(1L)
  input <- file.path(output, "cdf")
  dir.create(input, recursive = TRUE)
  sim <- simulate_GCMS_samples(samples = 16L, n = 1200L, k = 48L,
    preset = "time_area", presence = .5, irregular = 0, drift = 2,
    family = "quasipoisson", dispersion = 2, seed = 9118)
  truth <- lapply(sim$samples, function(s) {
    write_GCMS_cdf(s$Y, file.path(input, paste0(s$sample_id, ".cdf")),
      time = s$time, acquired_at = s$acquired_at)
    list(time = s$time, apex = apply(s$H_true, 2L, which.max))
  })
  P <- sim$P
  present <- t(sim$present)
  rm(sim)
  fits <- lapply(cores, function(ncores) {
    GCMSdeco(input, output_dir = file.path(output, paste0("workers", ncores)),
      family = quasipoisson("identity"), dispersion = 2, correlation = "none",
      nwindows = 2L, baseline_df = 4L, precluster = 32L, ncores = ncores,
      psf_control = list(maxiter = 3L),
      batch = list(coupling = list(cycles = 2L, initial_cycles = 1L, rounds = 1L, q = 32L)))
  })
  common <- fits[[1L]]$coupled
  for (z in fits) {
    stopifnot(isTRUE(all.equal(z$coupled$P, common$P, tolerance = 1e-9)),
      isTRUE(all.equal(z$areas, common$areas, tolerance = 1e-9)),
      identical(z$coupled$batch_grouping, common$batch_grouping))
  }

  # Match by shared spectrum and aligned apex, with one match per compound.
  K <- nrow(P); G <- nrow(common$P)
  full <- matrix(0, K, ncol(common$P))
  full[, match(colnames(P), common$mass)] <- P
  cos <- tcrossprod(common$P, full) / sqrt(outer(rowSums(common$P^2), rowSums(full^2)))
  rt <- vapply(seq_len(K), function(k) median(vapply(seq_along(truth), function(s) {
    z <- truth[[s]]
    common$alignment$maps[[s]]$evaluate(z$time[z$apex[k]])
  }, 0.)), 0.)
  estimated <- vapply(seq_len(G), function(g) {
    median(common$compounds$apex[common$compounds$root == g])
  }, 0.)
  gap <- abs(outer(estimated, rt, "-"))
  edges <- which(cos >= .9 & gap <= .15, arr.ind = TRUE)
  edges <- edges[order(-cos[edges], gap[edges]), , drop = FALSE]
  map <- rep(NA_integer_, G); used <- rep(FALSE, K)
  for (i in seq_len(nrow(edges))) {
    g <- edges[i, 1L]; k <- edges[i, 2L]
    if (is.na(map[g]) && !used[k]) { map[g] <- k; used[k] <- TRUE }
  }
  missed <- which(is.na(map)); label <- map
  label[missed] <- K + seq_along(missed)
  observed <- cbind(present, matrix(FALSE, nrow(present), length(missed)))
  recovered <- matrix(FALSE, nrow(observed), ncol(observed))
  recovered[, label] <- !is.na(common$areas)
  TP <- sum(observed & recovered); FP <- sum(!observed & recovered)
  FN <- sum(observed & !recovered); TN <- sum(!observed & !recovered)
  f1 <- 2*TP/(2*TP+FP+FN)
  mcc <- (as.double(TP)*TN-as.double(FP)*FN)/sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  out <- do.call(rbind, lapply(seq_along(fits), function(i) {
    z <- fits[[i]]
    data.frame(ncores = cores[i], workers = z$parallel$workers,
      clustering = z$seconds[['clustering']], initial = z$seconds[['initial']],
      coupled = z$seconds[['coupled']], exports = z$seconds[['exports']],
      total = z$seconds[['total']], GIC = GIC(z), shared = G,
      occurrences = sum(!is.na(z$areas)), TP = TP, FP = FP, FN = FN, TN = TN,
      F1 = f1, MCC = mcc,
      median_spectral_cosine = median(cos[cbind(which(!is.na(map)), na.omit(map))]))
  }))
  write.csv(out, file.path(output, "benchmark.csv"), row.names = FALSE)
  out
}

benchmark_samples16()
