# Benchmarks against other software

## Sixteen full simulated chromatograms

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

The native GCMSprime factor IDs mapped one-to-one to the 1,000
generating compounds across all eight positive samples per compound. For
MSHub, the correctly matched sample occurrences belonged to 1,072 of its
1,414 native components. [Occurrence-to-factor
assignments](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/factor-matches.csv).

Each program ran once, in sequence, on the same 64 GB Windows laptop
with eight physical cores and 16 logical processors. Seconds per sample
means wall time divided by 16 for the whole batch. The comparison
includes any required input conversion or staging and the program’s
native analysis. GCMSprime’s time includes reading and fitting an
independent PSF standard, the coupled fit and native results export.
Simulation, common-format result extraction and accuracy scoring are
outside the timed comparison. One run per program gives a descriptive
timing for these files and this machine.

[![Time per sample, sensitivity and precision for all seven
workflows](../reference/figures/simulated-16-performance.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-performance.png)

### Timing boundaries

| Workflow | Native analysis (s/sample) | Input conversion or staging (s/sample) | Through result extraction (s/sample) |
|----|---:|---:|---:|
| GCMSprime | 35.95 | — | 36.40 |
| MS-DIAL | 12.26 | 3.76 | 17.40 |
| mzmine GC-EI | 3.51 | 3.41 | 6.93 |
| eRah | 25.94 | — | 26.46 |
| AMDIS | 25.75 | 3.24 | 29.57 |
| xcms + CAMERA | 4.22 | — | 4.99 |
| MSHub | 14.99 | 0.06 | 15.38 |

The headline time adds native analysis and separately timed input
conversion or staging. A dash means the workflow read its inputs
directly. The last column also includes extracting peak tables and
spectra into the common in-memory representation, before saving that
record or comparing it with truth.

### Simulation and scoring

Individual chromatograms contain 479–517 true peaks. The L0glm library
supplies spectra on m/z 41–298. The selected spectra have compound
signal in 245 of the 550 recorded channels. A smoothly varying
background fitted to the measured gcmsData2 window supplies the
background on m/z 41–298. The additional channels carry low background.
Counts have quasi-Poisson variance proportional to the mean, with
dispersion 783.18 estimated from gcmsData2. Measured on the recorded
scan grids, the median within-compound Spearman correlation between area
and half-height width was 0.86 for EMGs and 0.94 for acids. The nominal
compound centres follow a regular grid with one deliberately close pair.
Concentration-dependent broadening and the subset present in each
injection determine the actual overlap. The retention-index (RI) scale
is synthetic and time-calibrated for each injection. Reference spectra,
occurrence labels and true rank are used in scoring.

GCMSprime calibrated the time-varying PSF on a separately generated
20-peak standard. It used a 16-thread budget with up to eight
simultaneous sample initialisations, 32 pooled mass-channel groups and
the final local rank search
(`batch = list(workers = 8L, coupling = list(penalty_search = TRUE))`).
The search compares changes to the original-ion GIC. It receives neither
the expected number of compounds nor their spectra. The other programs
used their local GC–EI workflows. The settings and software versions
accompany the results.

eRah used a minimum peak width of 2 seconds. The xcms matched filter
used a 2.6-second FWHM and a signal-to-noise threshold of 1.5. CAMERA
grouped peaks by FWHM with `perfwhm = 0.8` and required four ions. AMDIS
used `PEAKWIDTH = 16`, `SENSIT = 60`, `RESOLUTION = 0` and
`PEAKSHAPE = 0`. MSHub used its standalone deconvolution workflow. The
simulation and CDF writer ran with logconNMF and GCMSprime 0.1.0.9012.
Their generating functions, reference spectra and empirical background
were checked against the versions used for fitting (0.1.0.9013).

Each recovered spectrum competes against all 1,000 references, including
compounds absent from that sample. A true positive requires a reciprocal
best match with spectral cosine ≥0.70, within ±30 simulated RI units, to
a compound present in that sample. The score allows one credit per
sample and reference. Extra calls, including duplicates, count as false
positives. Sensitivity uses all 8,000 occurrences as its denominator.
Acid recovery uses all 80.

## Quantification and spectral recovery

| Workflow | Area adjusted R² | Median spectral cosine | Mean ion-support F1 | Median profile cosine |
|----|---:|---:|---:|---:|
| GCMSprime | 0.999982 | 0.9999995 | 0.995 | 0.999994 |
| MS-DIAL | 0.700129 | 0.9925397 | 0.901 | 0.996864 |
| mzmine GC-EI | 0.966560 | 0.9975148 | 0.874 | — |
| eRah | 0.853212 | 0.9957392 | 0.794 | 0.001678 |
| AMDIS | 0.499773 | 0.9972175 | 0.894 | 0.970942 |
| xcms + CAMERA | 0.951953 | 0.9807221 | 0.604 | 0.971638 |
| MSHub | 0.305059 | 0.9857457 | 0.913 | — |

Among GCMSprime’s 8,000 matched occurrences, the median absolute
relative area error was **0.095%** and ordinary adjusted R² on the log₁₀
scale was **0.999982**. The area panels use ordinary least-squares
regression with pointwise 95% prediction intervals. Square axes span 7.5
to 10.5 on both log₁₀ scales. The fits include every positive matched
pair, including observations outside the displayed limits. The dashed
diagonal marks exact area recovery.

[![Estimated versus true total-ion areas across all seven
workflows](../reference/figures/simulated-16-area-predictions.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-area-predictions.png)

Areas use each program’s exported ion intensities and its acquisition
clock to convert to integrated total-ion equivalents. Accuracy is
conditional on a correct match. The first table gives the missed and
extra calls alongside it.

[![Reciprocal reference spectrum matches across the 16
chromatograms](../reference/figures/simulated-16-spectral-violin.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-spectral-violin.png)

[![Ion-support recovery across the seven
workflows](../reference/figures/simulated-16-spectral-support.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-spectral-support.png)

The first distribution includes all reciprocal library hits. Crosses
identify hits to references absent from the sample. Ion-support F1
compares the nonzero ions at or above 0.1% of each spectrum’s own base
peak and is averaged over correct matches.

## Correctness of inferred elution profiles

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

CAMERA supplies a representative-ion trace scaled to the grouped area.
The mzmine and MSHub exports provide spectra and areas without a
sample-specific profile in this workflow.

### Recovery of all ten shark-fin acids

Each entry counts correctly identified occurrences out of eight for one
acid. The downloadable per-compound table also gives the median
estimated-to-true area ratio and mean matched profile cosine where a
profile is available.

| Acid ID | GCMSprime | MS-DIAL | mzmine | eRah | AMDIS | xcms+CAMERA | MSHub |
|---------|----------:|--------:|-------:|-----:|------:|------------:|------:|
| 250     |         8 |       7 |      5 |    5 |     8 |           3 |     8 |
| 311     |         8 |       7 |      3 |    7 |     8 |           6 |     4 |
| 372     |         8 |       8 |      6 |    6 |     8 |           1 |     8 |
| 433     |         8 |       7 |      4 |    8 |     8 |           2 |     8 |
| 494     |         8 |       7 |      4 |    7 |     8 |           3 |     0 |
| 556     |         8 |       7 |      7 |    8 |     8 |           5 |     8 |
| 617     |         8 |       8 |      0 |    8 |     8 |           2 |     6 |
| 678     |         8 |       8 |      0 |    8 |     8 |           7 |     0 |
| 739     |         8 |       7 |      0 |    7 |     8 |           2 |     8 |
| 800     |         8 |       7 |      0 |    6 |     8 |           0 |     0 |

## Reproduce the benchmark

Install the programs and set their local executable paths in the named
list `control`. The included R recipe generates the common CDFs and
their scoring key, fits the programs sequentially, then scores the saved
fits. [Engine setup and control
list](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/INSTALL.md).

``` r

bench <- system.file("benchmarks", package = "GCMSbenchmark", mustWork = TRUE)
source(file.path(bench, "simulated-16.R"))
directory <- "full-chromatograms"
prepareFullBenchmark(directory)
methods <- c("GCMSprime", "MS-DIAL", "mzmine", "MSHub", "AMDIS", "xcms+CAMERA", "eRah")
for (m in methods) fitFullBenchmark(directory, m, control[[m]], ncores = 16L)
for (m in methods) scoreFullBenchmark(directory, m)
collectFullBenchmark(directory, methods, file.path(directory, "results"))
source(file.path(bench, "plot-simulated-16.R"))
plotBenchmark16(file.path(directory, "results"), file.path(directory, "figures"))
```

[CDF
checksums](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/inputs.csv)
· [Reference
spectra](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/reference-library.rds)
· [Occurrence and area
matrices](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/settings.rds)
· [Summary and timing
components](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/benchmark.csv)
· [Matched
peaks](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/matches.csv)
· [Per-acid
recovery](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/acid-by-compound.csv)
· [Recovery by abundance
quartile](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/sensitivity-by-area.csv)
· [Software
versions](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/versions.csv)
· [Serial run
log](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/run-order.csv).
The numerical input hashes cover scan times, masses and intensities. The
whole-file hashes also include the CDF export timestamp.

## Comparison with other software

GCMSprime combines a formally **estimated time-varying point-spread
function (PSF)** with **compound-specific peak shapes**, **grouped peak
selection across mass channels** and a **mathematically specified
information criterion** for automatic compound-count selection.

For its default quasi-Poisson model, the fitted mean and variance
satisfy \operatorname{Var}(Y\_{it})=\phi\mu\_{it}. The final quasi-GIC
combines the corresponding fit criterion with a penalty for the number
of compounds and their occurrences. Grouped L0 selection uses the mass
channels jointly, and hard log-concavity constrains each component
profile. This combination distinguishes GCMSprime from every evaluated
configuration. The other workflows use peak detection, representative
ion traces, correlations, local covariance or factorisation of peak
integrals. Their preprocessing noise estimates and selection rules are
summarized below. [Mathematical
basis](https://tomwenseleers.github.io/GCMSprime/articles/algorithm.md)

| Workflow | How it deconvolves and where that can limit recovery | Access |
|:---|:---|:---|
| GCMSprime | Estimates a time-varying PSF and uses grouped L0 selection across ions to call peaks. A log-concave factor model fits compound-specific profiles, full spectra and smooth background, with GIC or quasi-GIC selecting the compound count. Samples can share spectra while retaining their own profiles. | Free, GPL open source. [Repository](https://github.com/tomwenseleers/GCMSprime) |
| MS-DIAL | Extracts spectra from detected model peaks within each run, then aligns the features. Smoothing, detection thresholds and the deconvolution window affect which components are recovered. [GC–MS tutorial](https://github.com/systemsomicslab/mtbinfo.github.io/blob/master/MS-DIAL/tutorial.md) | Free, LGPL-3 source. Some vendor readers have separate terms. [MsdialWorkbench](https://github.com/systemsomicslab/MsdialWorkbench#source-code-license) |
| mzmine / ADAP-GC | The benchmark’s GC-EI module groups resolved ions by retention time and shape similarity. Pseudospectra depend on upstream detection, the grouping thresholds and the minimum number of ions. Batch alignment follows extraction. [GC-EI guide](https://mzmine.github.io/mzmine_documentation/latest/module_docs/featdet_spectraldeconvolutiongc/spectraldeconvolutiongc.html) | Free Community licence for eligible academic users. Commercial use requires a paid plan. [Public source](https://github.com/mzmine/mzmine#license) · [Plans](https://mzio.io/plans/) |
| xcms + CAMERA | Detects ion peaks first and then groups them into pseudospectra. Filter widths and grouping rules affect the coverage of weak and coeluting peaks. CAMERA was developed for LC–MS feature annotation. [Kuhl et al. (2012)](https://doi.org/10.1021/ac202450g) | Free, GPL open source. [xcms](https://bioconductor.org/packages/xcms) · [CAMERA](https://bioconductor.org/packages/CAMERA) |
| eRah | Uses local covariance and orthogonal signal deconvolution after filtering and baseline removal. Recovery depends on width and signal thresholds. Cross-sample alignment follows per-sample deconvolution. [eRah manual](https://cran.r-project.org/web/packages/erah/vignettes/erah-manual.html) | Free, GPL open source. [CRAN](https://cran.r-project.org/package=erah) |
| AMDIS | Uses representative ion traces as model peak shapes within each chromatogram. Component width, sensitivity and shape settings govern the calls and extracted spectra. [NIST guide](https://chemdata.nist.gov/mass-spc/amdis/docs/amdis.pdf) | Free standalone binary. [NIST distribution](https://chemdata.nist.gov/dokuwiki/doku.php?id=chemdata%3Aamdis) |
| GNPS/MSHub | Requires a sample series for its NMF workflow. It factorizes aligned peak integrals, so the model resolves spectra and abundances from variation across samples and drops the scan-level profile information. [Aksenov et al. (2021)](https://doi.org/10.1038/s41587-020-0700-3) · [Workflow guide](https://ccms-ucsd.github.io/GNPSDocumentation/gc-ms-deconvolution/) | Free hosted workflow with a GNPS account. Public source under the repository’s nonprofit-use licence. [Source](https://github.com/CCMS-UCSD/GNPS_Workflows/tree/master/mshub-gc) · [Licence](https://github.com/CCMS-UCSD/GNPS_Workflows/blob/733ed89372c3754a4c4d132319da5b4313c47288/LICENSE) |
