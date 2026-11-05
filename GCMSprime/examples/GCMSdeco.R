library(L0nmf)
\dontrun{
library(L0glm)
data(GCMS_data1)
data(GCMS_data2)
out <- file.path(tempdir(), "L0glm-GCMS")
dir.create(out, showWarnings = FALSE)
write_GCMS_cdf(GCMS_data1, file.path(out, "GCMS_data1.cdf"))
write_GCMS_cdf(GCMS_data2, file.path(out, "GCMS_data2.cdf"))

fit1 <- GCMSdeco(file.path(out, "GCMS_data1.cdf"),
  file.path(out, "GCMS_data1.msp"), nwindows = 1L, precluster = 16L,
  baseline_df = 1L,
  correlation = "auto", ncores = 8L,
  psf_control = list(kernel_est = convkernel("gaussian", 65L,
    list(sigma = 8)), maxiter = 0L))
fit2 <- GCMSdeco(file.path(out, "GCMS_data2.cdf"),
  file.path(out, "GCMS_data2.msp"), nwindows = 1L, precluster = 32L,
  correlation = "auto", ncores = 8L,
  psf_control = list(kernel_est = convkernel("gaussian", 65L,
    list(sigma = 6)), maxiter = 0L))
fit1$compounds
GIC(fit2)
coef(fit2, normalize = "l1")

plot(fit1, type = "lines", nrows = 1L,
  minutes_per_panel = diff(range(fit1$time)) + .01, labels = TRUE)
plot(fit1, type = "heatmap")
plot(fit1, type = "profiles", nrows = 3L)
plot(fit2, type = "lines", file = file.path(out, "GCMS_data2_traces.pdf"),
  nrows = 1L, minutes_per_panel = diff(range(fit2$time)) + .01)
plot(fit2, type = "heatmap", file = file.path(out, "GCMS_data2_heatmaps.pdf"))
plot(fit2, type = "profiles", file = file.path(out, "GCMS_data2_profiles_and_spectra.pdf"))

# For a full chromatogram with a separate calibration standard:
beer <- GCMSdeco("beer.cdf", "beer.msp", alkanes_file = "alkanes.cdf",
  xlsx_file = "beer.xlsx", nwindows = 10L, precluster = 32L)

# Automatic PSF estimation from a bundled window and its named alkane TIC.
psf <- estimate_psf(GCMS_data2, nwindows = 3L)
dim(psf$X)

# A CDF contains the chromatogram; pass the R calibration attribute explicitly.
calibrated <- GCMSdeco(file.path(out, "GCMS_data2.cdf"),
  file.path(out, "GCMS_data2_calibrated.msp"), nwindows = 3L, precluster = 32L,
  psf_control = list(standard = attr(GCMS_data2, "alkanes_TIC")))
}

\dontrun{
# Sixteen related samples, retaining all mass channels. Compounds occur in
# half the samples on average, with retention drift and concentration-dependent
# widths and tails. Here the simulated quasipoisson dispersion is known.
library(L0nmf)
sim <- simulate_GCMS_samples(samples=16L, n=1200L, k=48L,
  preset="time_area", presence=0.5, irregular=0, drift=2,
  family="quasipoisson", dispersion=2, seed=9118)
input <- file.path(tempdir(), "GCMSprime_samples16")
dir.create(input, showWarnings=FALSE)
for(s in sim$samples) {
  write_GCMS_cdf(s$Y, file.path(input, paste0(s$sample_id, ".cdf")),
    time=s$time, acquired_at=s$acquired_at)
}
rm(sim, s)                         # observations are read from the CDF files
experiment <- GCMSdeco(input, output_dir=paste0(input, "_results"),
  family=quasipoisson("identity"), dispersion=2, correlation="none",
  nwindows=2L, baseline_df=4L, precluster=32L, ncores=16L,
  psf_control=list(maxiter=3L),
  batch=list(coupling=list(cycles=2L, initial_cycles=1L, rounds=1L, q=32L)))
experiment$seconds
experiment$parallel
experiment$samples[c("sample_name", "compounds", "GIC")]
experiment$areas
plot(experiment, type="areas")
plot(experiment, sample=1L, type="lines")
plot(experiment, sample=1L, type="profiles", which=1:6)
saveRDS(experiment$channel_groups, file.path(input, "channel_groups.rds"))
}
