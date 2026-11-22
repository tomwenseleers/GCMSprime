# Refit the two experimental examples, regenerate every downloadable result and
# draw the pkgdown figures. Run from the repository root:
#
# Rscript GCMSprime/pkgdown/build_single_examples.R [timing repetitions]

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the 'pkgload' package to rebuild the example outputs.",
    call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args)) as.integer(args[[1L]]) else 3L
if (length(repetitions) != 1L || is.na(repetitions) || repetitions < 1L) {
  stop("The timing repetition count must be a positive integer.", call. = FALSE)
}

package <- if (file.exists("GCMSprime/DESCRIPTION")) "GCMSprime" else "."
figures <- file.path(package, "man", "figures")
downloads <- file.path(package, "pkgdown", "assets", "gcms")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(downloads, recursive = TRUE, showWarnings = FALSE)

pkgload::load_all(package, quiet = TRUE)
data(gcmsData1, package = "logconNMF")
data(gcmsData2, package = "logconNMF")
stopifnot(
  identical(dim(gcmsData1), c(119L, 670L)),
  identical(dim(gcmsData2), c(1764L, 518L))
)

work <- tempfile("GCMSprime-single-examples-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
cdf1 <- file.path(work, "gcms-data1.cdf")
cdf2 <- file.path(work, "gcms-data2.cdf")
writeGCMSCDF(gcmsData1, cdf1)
writeGCMSCDF(gcmsData2, cdf2)

kernel1 <- logconNMF::convKernel("gaussian", 65L, list(sigma = 9))
kernel2 <- logconNMF::convKernel("gaussian", 65L, list(sigma = 6))

fit_hydrocarbon <- function(index) {
  gcmsDeco(
    cdf1,
    if (index == 1L) file.path(downloads, "gcms-data1.msp") else
      file.path(work, sprintf("gcms-data1-%d.msp", index)),
    xlsx_file = if (index == 1L)
      file.path(downloads, "gcms-data1.xlsx") else NULL,
    nwindows = 1L,
    precluster = 7L,
    baseline_df = 1L,
    sample_name = "gcmsData1",
    correlation = "auto",
    ncores = 8L,
    control = list(initial_relaxed = FALSE, fusion_passes = 0L),
    psf_control = list(kernel_est = kernel1, maxiter = 0L)
  )
}

fit_beer <- function(index) {
  gcmsDeco(
    cdf2,
    if (index == 1L) file.path(downloads, "gcms-data2.msp") else
      file.path(work, sprintf("gcms-data2-%d.msp", index)),
    xlsx_file = if (index == 1L)
      file.path(downloads, "gcms-data2.xlsx") else NULL,
    nwindows = 1L,
    precluster = 32L,
    sample_name = "gcmsData2",
    correlation = "auto",
    ncores = 8L,
    control = list(initial_relaxed = TRUE, fusion_passes = 0L),
    psf_control = list(kernel_est = kernel2, maxiter = 0L)
  )
}

fit1_runs <- lapply(seq_len(repetitions), fit_hydrocarbon)
fit2_runs <- lapply(seq_len(repetitions), fit_beer)
fit1 <- fit1_runs[[1L]]
fit2 <- fit2_runs[[1L]]
for (runs in list(fit1_runs, fit2_runs)) {
  ranks <- vapply(runs, function(f) ncol(f$model$H), 0L)
  stopifnot(all(ranks > 0L), length(unique(ranks)) == 1L,
            all(vapply(runs, function(f) nrow(f$model$S) == f$model$rank, TRUE)))
}

weighted_fit <- function(Y, fit) {
  model <- fit$model
  mu <- as.matrix(model$H %*% model$S) +
    as.matrix(model$background %*% model$alpha)
  offset <- model$offset
  if (length(offset) == 1L) mu <- mu + offset else mu <- mu + offset
  w <- 1 / pmax(mu, 0.1)
  means <- colSums(w * Y) / colSums(w)
  rss <- sum(w * (Y - mu)^2)
  tss <- sum(w * sweep(Y, 2L, means, "-")^2)
  r2 <- 1 - rss / tss
  p <- sum(model$S > 0) + sum(model$alpha > 0)
  n <- length(Y)
  c(
    R2_weighted = r2,
    adj_R2_weighted = 1 - (1 - r2) * (n - 1) / max(1, n - p - 1)
  )
}

example_row <- function(name, Y, runs) {
  fit <- runs[[1L]]
  quality <- weighted_fit(Y, fit)
  seconds <- vapply(runs, function(z) unname(z$seconds[["total"]]), 0)
  out <- data.frame(
    dataset = name,
    n = nrow(Y),
    m = ncol(Y),
    compounds = ncol(fit$model$H),
    seconds = stats::median(seconds),
    GIC = unname(logconNMF::GIC(fit$model)),
    R2_weighted = unname(quality[["R2_weighted"]]),
    adj_R2_weighted = unname(quality[["adj_R2_weighted"]]),
    rho = fit$model$spec$correlation,
    kappa = fit$model$spec$information_factor
  )
  for (i in seq_along(seconds)) out[[paste0("seconds_run", i)]] <- seconds[[i]]
  out
}

results <- rbind(
  example_row("GCMS_data1", gcmsData1, fit1_runs),
  example_row("GCMS_data2", gcmsData2, fit2_runs)
)
utils::write.csv(results, file.path(downloads, "example-results.csv"),
  row.names = FALSE)

png_plot <- function(path, width, height, draw) {
  grDevices::png(path, width = width, height = height, res = 180,
    bg = "white", type = if (.Platform$OS.type == "windows") "windows" else "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
}

# PDF downloads contain all compounds and retain vector spectra.
plot(fit1, type = "lines", file = file.path(downloads, "gcms-data1-traces.pdf"),
  nrows = 1L, minutes_per_panel = diff(range(fit1$time)) + 0.01,
  labels = TRUE)
plot(fit1, type = "heatmap", file = file.path(downloads,
  "gcms-data1-heatmaps.pdf"))
plot(fit1, type = "profiles", file = file.path(downloads,
  "gcms-data1-profiles-spectra.pdf"), nrows = 3L)
plot(fit2, type = "lines", file = file.path(downloads, "gcms-data2-traces.pdf"),
  nrows = 1L, minutes_per_panel = diff(range(fit2$time)) + 0.01,
  labels = TRUE)
plot(fit2, type = "heatmap", file = file.path(downloads,
  "gcms-data2-heatmaps.pdf"))
plot(fit2, type = "profiles", file = file.path(downloads,
  "gcms-data2-profiles-spectra.pdf"), nrows = 3L)

# Raster previews use the same fitted objects as the downloadable reports.
png_plot(file.path(figures, "gcms-data1-traces.png"), 1800, 600,
  function() plot(fit1, type = "lines", nrows = 1L,
    minutes_per_panel = diff(range(fit1$time)) + 0.01, labels = TRUE))
png_plot(file.path(figures, "gcms-data1-heatmaps.png"), 1800, 2200,
  function() plot(fit1, type = "heatmap"))
png_plot(file.path(figures, "gcms-data1-profiles-spectra.png"), 1800, 2200,
  function() plot(fit1, type = "profiles", nrows = 6L))
png_plot(file.path(figures, "gcms-data2-traces.png"), 1800, 600,
  function() plot(fit2, type = "lines", nrows = 1L,
    minutes_per_panel = diff(range(fit2$time)) + 0.01, labels = TRUE))
png_plot(file.path(figures, "gcms-data2-heatmaps.png"), 1800, 2200,
  function() plot(fit2, type = "heatmap"))

H2 <- as.matrix(fit2$model$H)
profile_mean <- colSums(H2 * fit2$time)
profile_sd <- sqrt(colSums(H2 * sweep(matrix(fit2$time,
  nrow(H2), ncol(H2)), 2L, profile_mean, "-")^2))
largest <- head(order(-rowSums(as.matrix(fit2$model$S))), 6L)
broadest <- which.max(profile_sd)
if (!broadest %in% largest) largest[[6L]] <- broadest
representative <- sort(unique(largest))
stopifnot(length(representative) == 6L)
png_plot(file.path(figures, "gcms-data2-profiles-spectra.png"), 1800, 2200,
  function() plot(fit2, type = "profiles", which = representative,
    nrows = 6L))

# Compact hydrocarbon overview used by both package homepages.
png_plot(file.path(figures, "gcms-data1-decomposition.png"), 2100, 1400,
  function() {
    Y <- gcmsData1
    model <- fit1$model
    background <- as.matrix(model$background %*% model$alpha)
    corrected <- pmax(Y - background, 0)
    profiles <- sweep(as.matrix(model$H), 2L,
      apply(model$S, 1L, max), "*")
    spectra <- sweep(as.matrix(model$S), 1L,
      apply(model$S, 1L, max), "/") * 100
    apex <- model$apex_scan
    # Match plot.GCMSdeco() so every component keeps the same colour in the
    # compact overview and the detailed profile/spectrum panels below it.
    colours <- grDevices::hcl(
      h = (seq_len(model$rank) * 137.508) %% 360,
      c = 75,
      l = 48
    )
    # Use the matrix dimnames as the canonical axes. writeGCMSCDF() round-trips
    # these to numerical precision, but the original values also make this
    # compact rendering byte-for-byte reproducible with the logconNMF build.
    mass <- as.numeric(colnames(Y))
    time <- as.numeric(rownames(Y))
    panels <- matrix(0L, ceiling(model$rank/3), 3L)
    panels[t(row(panels) + (col(panels)-1L)*nrow(panels))[seq_len(model$rank)]] <- 1L + seq_len(model$rank)
    graphics::layout(rbind(rep(1L, 3L), panels), heights=c(1.85, rep(1, nrow(panels))))
    top_max <- max(corrected, profiles)
    raw_ticks <- pretty(c(0, top_max), n = 5L)
    raw_ticks <- raw_ticks[raw_ticks >= 0 & raw_ticks <= top_max]
    graphics::par(mar = c(5.2, 6.0, 3.4, 1.2),
      mgp = c(3.6, 0.9, 0), las = 1)
    graphics::matplot(time, sqrt(corrected), type = "l", lty = 1,
      lwd = 0.45, col = grDevices::adjustcolor("grey45", alpha.f = 0.16),
      xlab = "Retention time (min)", ylab = "Intensity",
      ylim = sqrt(c(0, top_max * 1.05)), yaxt = "n",
      main = sprintf("gcmsData1 | square-root display | %d fitted compounds", model$rank))
    graphics::axis(2, at = sqrt(raw_ticks),
      labels = format(raw_ticks, scientific = TRUE), las = 1)
    graphics::abline(h = 0, col = "grey70", lty = 3)
    for (k in seq_len(model$rank)) {
      graphics::lines(time, sqrt(profiles[, k]), col = colours[k], lwd = 2.3)
      graphics::text(time[apex[k]], sqrt(profiles[apex[k], k]), k,
        col = colours[k], font = 2, pos = 3, cex = 0.85, offset = 0.35)
    }
    for (k in seq_len(model$rank)) {
      spectrum <- spectra[k, ]
      keep <- which(spectrum >= 0.1)
      graphics::par(mar = c(3.8, 3.8, 3.1, 0.8),
        mgp = c(2.25, 0.72, 0), las = 1)
      graphics::plot(range(mass), c(0, 100), type = "n", xlab = "m/z",
        ylab = "Relative intensity (%)",
        main = sprintf("Component %d | apex %.3f min", k, time[apex[k]]),
        cex.main = 0.83, cex.axis = 0.78, cex.lab = 0.82, xaxs = "i")
      graphics::segments(mass[keep], 0, mass[keep], spectrum[keep],
        col = colours[k], lwd = 1.15)
    }
  })

print(results, digits = 10, row.names = FALSE)
cat("Hydrocarbon apices (min):",
  paste(format(fit1$compounds$retention_time, digits = 7L), collapse = ", "),
  "\n")
cat("Beer preview components:", paste(representative, collapse = ", "), "\n")
