# GC/MS deconvolution for multiple samples

[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
accepts an experiment directory and analyzes its chromatograms together.
The batch can share PSF calibration. Each sample has its own native scan
grid, smooth background and log-concave compound profiles. Simultaneous
factorisation shares mass spectra across samples and returns one
consensus compound library with a sample-by-compound area table.

## Organize and fit an experiment

``` text
chromatograms/
  standards/
    alkanes.cdf
  beer1/
    injection001.cdf
    injection002.cdf
  beer2/
    injection003.cdf
```

``` r

library(GCMSprime)

experiment <- gcmsDeco("chromatograms", output_dir="beer_results",
  family=quasipoisson("identity"), precluster=32L,
  batch=list(plots=TRUE))

experiment$samples
experiment$areas
coef(experiment, normalize="l1")
GIC(experiment)
plot(experiment, sample=1L, type="lines", file="sample1_traces.pdf")
```

The `standards` subdirectory is reserved for calibration chromatograms.
Analytical samples are sorted using their acquisition timestamps,
including time-zone information. Group labels come from their relative
directories. Each sample uses the nearest alkane standard in acquisition
time. A single standard can calibrate the whole experiment;
`alkanes_file=` explicitly selects a particular calibration.

For files without acquisition timestamps, supply their relative paths in
the intended run order:

``` r

experiment <- gcmsDeco("chromatograms", output_dir="beer_results",
  batch=list(sample_order=c("beer1/injection001.cdf",
                            "beer1/injection002.cdf",
                            "beer2/injection003.cdf")))
```

## Sixteen full chromatograms

The same 16 simulated chromatograms were supplied to every program. Each
contains **10,000 scans × 550 nominal mass channels** (m/z 30–579),
covering approximately 50 minutes. Together they contain **88 million
simulated ion intensities**. Each of the **1,000 distinct compounds**
occurs in eight randomly chosen samples. The design therefore contains
**8,000 positive sample–compound pairs**. Rank here counts the distinct
generating compounds across the batch, each with a time-local elution
profile in each positive sample. Of the compounds, 990 have
exponentially modified Gaussian (EMG) peaks and ten have a rising
shark-fin profile with an abrupt trailing edge. Higher concentrations
broaden both types. The acquisition clock allows a smooth retention
shift of up to 2.93 seconds between the first and last injections and
approximately 2% variation in scan spacing.

GCMSprime, using simultaneous logconNMF, **recovered the true rank of
1,000 exactly** in **35.95 seconds per sample**. It recovered **8,000 of
8,000** compound occurrences and made **no false-positive calls**. The
recovered peaks included **80 of the 80 acid occurrences**. MSHub
returned **1,414 deconvolved components**, 41.4% above the true rank.
The other five workflows produce peak lists or aligned feature groups.
Their outputs do not define a common factorisation rank.

| Workflow | Seconds/sample | Estimated rank (true 1,000) | TP | FP | FN | Acid TP / 80 |
|----|---:|---:|---:|---:|---:|---:|
| GCMSprime | 35.95 | **1,000 (exact)** | 8,000 | 0 | 0 | 80 / 80 |
| MS-DIAL | 16.02 | — | 7,728 | 11,685 | 272 | 73 / 80 |
| mzmine GC-EI | 6.92 | — | 4,202 | 185 | 3,798 | 29 / 80 |
| eRah | 25.94 | — | 7,172 | 5,615 | 828 | 70 / 80 |
| AMDIS | 28.99 | — | 7,878 | 9,122 | 122 | 80 / 80 |
| xcms + CAMERA | 4.22 | — | 7,261 | 890 | 739 | 31 / 80 |
| MSHub | 15.05 | 1,414 | 7,179 | 12,184 | 821 | 50 / 80 |

### Local search for the shared rank

GCMSprime calibrated the time-varying PSF on a separately generated
20-peak standard. It used a 16-thread budget with up to eight
simultaneous sample initialisations, 32 pooled mass-channel groups and
the final local rank search
(`batch = list(workers = 8L, coupling = list(penalty_search = TRUE))`).
The search compares changes to the original-ion GIC. It receives neither
the expected number of compounds nor their spectra. The other programs
used their local GC–EI workflows. The settings and software versions
accompany the results.

Each program ran once, in sequence, on the same 64 GB Windows laptop
with eight physical cores and 16 logical processors. Seconds per sample
means wall time divided by 16 for the whole batch. The comparison
includes any required input conversion or staging and the program’s
native analysis. GCMSprime’s time includes reading and fitting an
independent PSF standard, the coupled fit and native results export.
Simulation, common-format result extraction and accuracy scoring are
outside the timed comparison. One run per program gives a descriptive
timing for these files and this machine. The observed peak summed
working set of the GCMSprime process tree was 6.8 GiB, sampled every
three seconds. The final original-ion GIC was 32,610,830. [Stage times
per
sample](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gcmsprime-stages.csv).

Among GCMSprime’s 8,000 matched occurrences, the median absolute
relative area error was **0.095%** and ordinary adjusted R² on the log₁₀
scale was **0.999982**.

The examples are the EMG compound and the acid with the greatest total
true area in their respective classes. Both were selected from the
simulation before comparing any fits. Each occurs in eight samples.
Solid curves show inferred profiles, dashed curves show the true
profiles, and the eight absent true profiles coincide at zero. Colours
order the 16 samples by true area. Heights retain original total-ion
intensity per scan. The mean cosine in each example uses all eight
positive occurrences, with a missing recovered profile scoring zero. The
all-acid tally covers all ten acids.

[![Recovered and simulated profiles for the selected EMG
compound](../reference/figures/simulated-16-regular-profiles.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-regular-profiles.png)

[![Recovered and simulated profiles for the selected shark-fin
acid](../reference/figures/simulated-16-acid-profiles.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-acid-profiles.png)

[Seven-workflow comparison and reproduction
code](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.html).

### Optimality gap to ground truth-initialised solution

The blind fit descended from an original-ion GIC of **46,273,866** to
**32,610,830**. The fit initialised from the true elution profiles
finished at **32,682,079**, giving a signed empirical optimality gap of
**-0.2180%** (the blind fit scored below the reference). Both
trajectories showed smooth empirical descent. GIC never increased across
the 29 blind or 16 reference checkpoints.

Both fits used the same noisy counts, penalty and stopping rules. Only
the reference received the true unit-area elution shapes and their
cross-sample identities. Both fits then updated the elution profiles,
their amplitudes, the shared spectra and the background. The red line
gives the best GIC reached from the ground truth initialisation under
these stopping rules. The global minimum is unknown.

[![Blind and ground truth-initialised GIC trajectories on the full
16-sample
benchmark](../reference/figures/full-benchmark-gic-gap.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/full-benchmark-gic-gap.png)

The blind fit started with 1,088 candidate factors. It crossed below the
reference after the support exchange at sweep 6 and finished after 8
sweeps. The reference stopped after 4. Both ended at rank 1,000 with
8,000 occurrences. The reference had 143,787 less scaled deviance. The
blind fit used 3,961 fewer profile degrees of freedom and 20 fewer
spectral degrees of freedom, which lowered its final penalised criterion
by 71,249.

The signed gap is 100 × (blind GIC − reference GIC) / \|reference GIC\|.
The left panel shows every recorded original-ion GIC. The right panel
enlarges the final half of each trajectory on the percentage scale. The
red dotted line marks the reference endpoint. The published 35.95
seconds per sample measures the ordinary blind run. Collecting both
traces was a separate diagnostic run.

[Block-level GIC
trace](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gic-trace.csv)
· [Final GIC
values](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gic-gap.csv).

## Outputs and memory use

The output directory contains:

- `samples.csv`: acquisition order, group, calibration file and sample
  result paths.
- `areas.csv`: integrated compound areas; undetected occurrences have
  `NA`.
- `consensus.msp`: the common spectral library, using the shared 0.1%
  base-peak cutoff.
- `experiment.xlsx`: sample information, consensus compounds, packed
  spectra and areas.
- `experiment.rds`: the experiment object.
- `samples/<sample_id>/`: fitted sample RDS, NIST-compatible MSP and
  Excel results, plus trace and heatmap PDFs when `batch$plots=TRUE`.

Spectral text cells are unwrapped and list fragments by decreasing
intensity. Sample objects retain factor matrices and CDF paths. Grouped
PSF regressions divide the core budget between independent samples and
native threads. The scheduler estimates how many regressions fit in
available memory, including Windows commit headroom. After those workers
exit, a single process shares the original ion counts and fitted means.
OpenMP updates elution profiles across samples. It fits each ion’s
shared spectrum and sample-specific background coefficients
independently, keeping a fixed component order within each ion. The
general identity-link solver can stream original-channel blocks from
disk. Use `batch=list(workers=4L)` to cap file-import and PSF worker
counts. For this workflow, `experiment$parallel$workers` and
`experiment$parallel$threads` describe the grouped PSF regressions.
`ncores` sets the total thread budget for joint refinement.

The same ion groups are used across samples for matching. They are
learned from a pooled, time-binned chromatogram or supplied through
`precluster` as a saved [channel
map](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.html),
for example one trained on reference spectra. Each sample’s PSF and scan
grid remain separate.

``` r

sample1 <- readRDS(experiment$paths[1L])
sample1$compounds
plot(sample1, type="profiles", file="sample1_profiles_spectra.pdf")
fitted(experiment, sample=1L)
```

The [simultaneous NMF
tutorial](https://tomwenseleers.github.io/logconNMF/articles/simultaneous-nmf.html)
gives the mathematical model, alignment, spectral grouping, local
updates and the full 16-sample example, including the two-start GIC
diagnostic and reproduction code.

## Scan spacing and exposure

Observed matrices stay on their native scan grids. Monotone alignment
changes profile coordinates and bin widths; its Jacobian preserves
integrated peak areas. The model is evaluated directly at the original
scan times.

By default, GCMSprime models the vendor’s exported intensities with unit
exposure. Quadrupole acquisition distributes scan time across the mass
range. The effective time per mass also depends on scan step size, scan
averaging and acquisition settings; see the [Agilent method
documentation](https://openlab.help.agilent.com/en/mergedProjects/GCMSHelp/MEAcquiUI.htm).
ANDI files can contain arbitrary abundance units, so scan duration and
the number of exported mass peaks alone do not identify physical count
exposure.

If integrated counts and their durations are known, supply exposure
explicitly:

``` r

experiment <- gcmsDeco("chromatograms", output_dir="count_results",
  family=poisson("identity"),
  batch=list(exposure=function(metadata, time) metadata$scan_duration))
```

Here `scan_duration` is appropriate only when it represents the
integration duration of each reported count. A calibrated per-mass
duration can instead be returned by the function. Metadata vectors are
restricted to the fitted scans. With an identity link, exposure
multiplies the latent mean. It is incorporated into profile and
background design values during fitting.

## Calibrated beer example

An alkane standard supplies narrow peak shapes for both full
chromatograms. The conditional count workflow estimates a pooled
time-varying PSF and transfers it to each sample’s observed scan times.
In this example, `calibration_range=c(6, Inf)` excludes the solvent and
early contaminated region from PSF calibration; the analytical samples
retain their full time range and all 518 measured mass channels.

``` r

standard <- "chromatograms/standards/alkanes.cdf"
experiment <- gcmsDeco("chromatograms", output_dir="beer_results",
  alkanes_file=standard, calibration_range=c(6, Inf), precluster=32L,
  ncores=8L, batch=list(workers=2L))

# Use the calibrated narrow-peak width for onset-based RI coordinates.
reference <- readRDS(experiment$paths[
  grepl("54040014_a", experiment$paths)][[1L]])
cal <- calibrateRI(standard, carbon=7:23, psf=reference$psf_standard)
experiment <- annotateRI(experiment, cal, phase="DB5",
  msp_file="beer_results/consensus.msp",
  xlsx_file="beer_results/experiment.xlsx")
```

A 26-column spline basis with ion-specific nonnegative coefficients
allows the background spectrum to change through each chromatogram. The
fit used four threads per worker and retained **643 shared identities
and 1,076 occurrences**.

| Sample         | Occurrences | Dispersion | Sample quasi-GIC | Adjusted weighted R² |
|:---------------|------------:|-----------:|-----------------:|---------------------:|
| 502131900304_a |         549 |     631.78 |       14,310,415 |              0.98416 |
| 54040014_a     |         527 |     617.21 |       15,068,092 |              0.98345 |

The joint quasi-GIC was **28,467,211**, with each shared spectrum
charged once. The two beer CDFs took **339.26 s** from input through
export. Preparation and PSF calibration took 145.32 s, grouped
regression and joint refinement 186.38 s, and result construction and
export 6.94 s. RI calibration added 3.11 s and RI annotation and export
added 4.36 s, bringing the total to **346.73 s** before plotting.

Using the Agilent NIST records, a ±30 RI window and a minimum spectral
cosine of 0.85, the two beers yielded **212 and 203 distinct reciprocal
identifications**. There were 88 and 73 additional forward assignments
to identities already called in the respective samples. These are
tentative identifications in beers of unknown composition.

Diagnostic ions identify all 17 alkanes from C7 to C23. The retention
coordinate is the 1%-height onset plus half the local narrow-peak width
at that height. The natural cubic RI curve passes through all accepted
alkane anchors exactly. MSP files contain phase-labelled DB5 indices;
the consensus library uses their medians. Values outside the accepted
calibration range remain unassigned.

The full beer fits still contain residual fragmentation and occasional
overly broad profiles in difficult overlapping regions. Component counts
are estimates under the fitted model, rather than a verified inventory
of chemical compounds.

The following panels zoom into the 4.5–6 minute region of these
full-chromatogram fits. They are not separately fitted cropped windows.
Grey lines show observed ions after subtracting the jointly modelled
background; coloured lines are the fitted compounds.

![Local full-fit profiles for beer sample
502131900304_a](../reference/figures/beer-502131900304_a-local.png)

Local full-fit profiles for beer sample 502131900304_a

![Local full-fit profiles for beer sample
54040014_a](../reference/figures/beer-54040014_a-local.png)

Local full-fit profiles for beer sample 54040014_a

The fitted components can be inspected sample by sample, and common
spectra can be searched against a library with each sample’s
alkane-derived retention indices. Pass `library=` and `search=` to
[`gcmsDeco()`](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.md)
to include those identifications in per-sample MSP and Excel exports.

Set `coupled=FALSE` to export independent fits. Set
`batch=list(resume=TRUE)` to reuse independently fitted samples whose
source-file metadata and fitting settings match the cache. Simultaneous
refinement is then computed from those fits.

See
[GCMSbatch](https://tomwenseleers.github.io/GCMSprime/reference/GCMSbatch.md)
for options and [single-sample
deconvolution](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.md)
for diagnostic plots, profile fitting and spectral recovery.
