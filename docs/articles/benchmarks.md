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
identify hits to references absent from the sample. All cosines use
their measured values on the logit scale. Of the 8,000 GCMSprime
spectra, 5,656 scored above 0.999999 and the maximum was 0.9999999913.
Exact ion-support matches occupy a separate position at 1, with the
violin showing the remaining observations. Ion-support F1 compares the
nonzero ions at or above 0.1% of each spectrum’s own base peak and is
averaged over correct matches.

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

## Optimality gap to ground truth-initialised solution

The blind fit descended from **46,273,866** to **32,610,830**. The fit
initialised with the true elution profiles reached **32,672,206** after
the same eight sweeps. The blind fit finished **0.1879% below** that
reference. GIC never increased across the 29 blind or 25 reference
checkpoints. Both fits ended at rank **1,000** with **8,000
occurrences**.

Only the reference received the true unit-area profiles and compound
identities. Both fits estimated spectra, areas, profiles and background
from the same noisy counts. The open circle marks the fourth sweep, when
the reference met the ordinary stopping rule. The red line marks its
final GIC after four further full sweeps. These trajectories give an
empirical comparison. The global minimum is unknown.

[![Blind and ground truth-initialised GIC trajectories on the full
16-sample
benchmark](../reference/figures/full-benchmark-gic-gap.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/full-benchmark-gic-gap.png)

The curves cross twice. At the start, the reference already fitted the
ion counts better, with scaled deviance 19.56 million versus 30.90
million. Its initial spectral regression nevertheless used 542,801 free
ion coefficients, compared with 77,974 for the blind start. The larger
penalty put the reference above the blind fit. Its first spectral update
reduced that count to 73,871. By the end of sweep one the reference GIC
was 32,976,439 and the blind GIC was 36,473,383.

The reference met the 0.1% relative objective-change threshold at sweep
four. The blind fit needed six initial sweeps, accepted support changes
and made two further refits. Extending the reference to eight full
sweeps reduced its GIC from 32,682,079 to 32,672,206. The last sweep
still gained 819 GIC units. The blind fit overtook it again during the
profile update of sweep seven. At the endpoint the reference had 150,960
less scaled deviance. The blind fit used 3,915 fewer profile degrees of
freedom and 16 fewer spectral degrees of freedom, which gave it the
lower penalised GIC by 61,376.

The signed gap is 100 × (blind GIC − reference GIC) / \|reference GIC\|,
giving **-0.1879%**. The right panel enlarges the last four sweeps.
Individual support changes can lower GIC at the same sweep number. The
ordinary blind benchmark took 35.95 seconds per sample. The two-start
diagnostic was timed separately.

[Block-level GIC
trace](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gic-trace.csv)
· [Final GIC
values](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gic-gap.csv).

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
combines the dispersion-scaled fit criterion with charges for the free
nonzero ion proportions in each shared spectrum and the fitted dimension
of each log-concave profile. Grouped L0 selection uses the mass channels
jointly, and hard log-concavity constrains each component profile. This
combination distinguishes GCMSprime from every evaluated configuration.
The other workflows use peak detection, representative ion traces,
correlations, local covariance or factorisation of peak integrals. Their
preprocessing noise estimates and selection rules are summarized below.
[Mathematical
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

### Commercial GC–MS deconvolution software

Two decisions matter when an ion trace contains a shoulder. The software
must decide which ions belong to the same elution profile and whether a
second compound explains enough signal to justify another peak. The
commercial workflows below make some of these decisions visible through
peak-width, correlation or detection settings. Several also use
information from multiple samples.

| Software | Documented workflow and practical implications |
|:---|:---|
| **SpectralWorks AnalyzerPro** | Users set detection thresholds for peak area, height, width and signal-to-noise, along with scan windows and chromatographic resolution. Detected components can populate a target library across a sequence. A component missed during extraction contributes no spectrum from that file to the library. [AnalyzerPro metabolomics workflow](https://spectralworks.com/2016/07/08/using-analyzerpro-in-the-metabolomics-workflow/) |
| **Agilent MassHunter Unknowns Analysis** | This module supports untargeted deconvolution and retention-index library searches. Its retention-time window factor controls how extracted ion peaks are grouped. Agilent documents a resolution–runtime trade-off. Smaller settings separate closer peaks, return more components and take longer. Larger settings suit broader peaks, and several widths can be searched together. [Unknowns Analysis guide, p. 27](https://www.agilent.com/cs/library/usermanuals/public/G3335-90243_Unknowns_Analysis_Familiarization-en.pdf#page=27) |
| **Thermo Scientific TraceFinder** | In the documented unit-resolution plug-in (v1.3 for TraceFinder 4.1), a user-set minimum ion overlap determines which extracted peaks form a component. A threshold on the summed ion peak controls entry in the output list before optional library searching. Both settings therefore affect which candidates reach identification. [Thermo technical note 10624, pp. 4–5](https://documents.thermofisher.com/TFS-Assets/CMD/Technical-Notes/tn-10624-gc-ms-deconvolution-plugin-tracefinder-tn10624-en.pdf#page=4) |
| **Shimadzu LabSolutions Insight Explore GCMS** | Peak detection, deconvolution and library searching can run across a batch. The software extracts component spectra even for co-elutions represented by a single integrated peak. Shimadzu describes the batch workflow and its outputs. How it chooses a shared compound count across injections is unspecified. [Shimadzu features](https://www.shimadzu.com/an/products/software-informatics/software-option/labsolutions-insight-explore-gcms/features.html) |
| **Markes / SepSolve ChromSpace 1D** | Supports multiple input formats, automated background subtraction, deconvolution and library searching. Markes lists automatic alignment and between-sample comparison under the ChromCompare+ upgrade. ChromSpace extends processing to GC×GC. [ChromSpace 1D specifications](https://markes.com/separationanalysisdata/data-analysis/chromspace-1d) |
| **LECO ChromaTOF / ChromaTOF Sync** | Sync explicitly uses evidence across the entire sample set for peak discovery, spectral improvement and a single aligned peak table. This makes it a particularly relevant multisample comparator. LECO specifies Pegasus BT-series data for Sync and a separate Sync 2D workflow for GC×GC. The published description leaves its component-count criterion unspecified. [LECO ChromaTOF Sync](https://www.leco.com/instrument/chromatof-sync-data-alignment/) |
| **Wiley KnowItAll GC Expert / MS Expert** | Wiley’s worked examples deconvolve a chromatogram, search each extracted spectrum and review the proposed identities. A library hit supports the identity of an extracted component. Peak splitting and the allocation of total area require separate checks against known profiles. [GC Expert tutorial](https://sciencesolutions.wiley.com/wp-content/uploads/2025/09/07_Training-GC-MS_Process_And_Analysis.pdf) · [MS Expert tutorial](https://sciencesolutions.wiley.com/wp-content/uploads/2024/12/07_Training-GC-MS_Process_And_Analysis.pdf) |
| **MsMetrix GC-Analyzer / GCxGC-Analyzer** | GC-Analyzer groups ions by retention time and shape, with differential analysis for sample–control comparisons. GCxGC-Analyzer offers total-ion or individual-ion processing for two-dimensional separations. The choice of separation and detection mode matters when comparing the resulting peak counts. [GC-Analyzer](https://www.msmetrix.com/our-products/gc-analyzer/) · [GCxGC-Analyzer](https://www.msmetrix.com/our-products/gcxgc-analyzer/) |
| **ACD/Labs Spectrus (MS Structure ID Suite)** | IntelliXtract follows ion traces to separate co-eluting components. Spectrus also offers targeted screening, library and structure tools, and adjustable peak detection. The extraction and structure-assignment steps answer different questions, and the documented outputs can be assessed separately for peak recovery and identification. [ACD/Labs MS Structure ID Suite](https://www.acdlabs.com/products/spectrus-platform/ms-structure-id-suite/) |

An independent experiment illustrates why component width deserves
attention. In GC–TOF-MS metabolite mixtures, Lu et al. found that
AnalyzerPro returned 28 spectra with library match factors above 850 at
500 µM and five at 50 µM. They identified component width as the most
influential processing setting across the three programs they tested.
Those results describe the versions and mixtures studied in 2008. [Lu et
al. (2008)](https://doi.org/10.1016/j.trac.2007.11.004)

GCMSprime addresses peak grouping by fitting all original ion channels
with one spectrum per compound across samples. Each occurrence has its
own log-concave profile, so a resolved occurrence in one injection can
inform a co-eluting occurrence in another. The calibrated point spread
function supplies starting widths. The subsequent fit can change widths
and asymmetry, and its GIC charges for both spectral coefficients and
profile complexity. That makes the evidence for retaining a split
inspectable. The full benchmark recovered all 1,000 compounds and all
8,000 occurrences with zero false calls. Recovery, area error and
runtime for the commercial products on these same files remain
unmeasured.

### Matrix and tensor factorisation approaches

The models below describe the methodological alternatives to the seven
timed workflows.

| Approach | What is shared and how peaks are represented | How the component count enters |
|:---|:---|:---|
| **GCMSprime / logconNMF** | A common spectrum across samples and one log-concave elution profile per occurrence. Width, asymmetry, area and retention time can vary by sample. | Grouped PSF regression proposes peaks. Original-ion GIC scores changes to the shared factors and their support. |
| **MCR-ALS, conventional GC–MS NMF and unimodal NMF** | Alternating regressions estimate elution profiles and spectra. Nonnegativity and unimodality can be imposed. Poisson MCR also models signal-dependent noise. [pyMCR](https://doi.org/10.6028/jres.124.018) · [GC–MS NMF](https://doi.org/10.1016/j.talanta.2004.09.017) · [Unimodal NMF](https://doi.org/10.1109/ICASSP39728.2021.9414631) · [Poisson MCR](https://doi.org/10.1016/j.chemolab.2016.02.006) | The cited fits use a chosen component count. Diagnostics can guide that choice. |
| **Simultaneous NMF and global parametric GC–MS analysis** | Simultaneous NMF shares a factor across data matrices. Parametric GC–MS analysis fits shared spectra with an explicit elution model such as an exponentially modified Gaussian. [Shared-factor NMF](https://doi.org/10.1371/journal.pone.0048238) · [Global GC–MS analysis](https://doi.org/10.1016/j.chemolab.2008.10.004) | The analyst specifies a candidate model order. |
| **PARAFAC and shifted PARAFAC** | Trilinear components share their spectrum and elution shape up to amplitude. Shifted PARAFAC also allows translations of that shape. [Shifted factor models](https://three-mode.leidenuniv.nl/pdf/h/harshmanhonglundy2003i.pdf) | Candidate ranks can be compared using diagnostics such as core consistency. [Bro and Kiers](https://doi.org/10.1002/cem.801) |
| **Classical PARAFAC2 and PARADISe** | Elution profiles can vary across samples subject to a common matrix of cross-products. PARADISe applies PARAFAC2 to chromatographic intervals. [PARAFAC2 for retention shifts](https://research.rug.nl/en/publications/parafac2-part-ii-modeling-chromatographic-data-with-retention-tim/) · [PARADISe](https://doi.org/10.1016/j.chroma.2017.04.052) | Candidate ranks are fitted per interval. Automated model assessment is also available. [Johnsen et al.](https://doi.org/10.1002/cem.2575) |
| **Flexible and fully constrained PARAFAC2** | Flexible coupling relaxes the cross-product condition. AO-ADMM can retain the condition while imposing constraints on all modes, including the varying elution profiles. [Cohen and Bro](https://doi.org/10.1007/978-3-319-93764-9_9) · [Roald et al.](https://doi.org/10.1137/21M1450033) | Both formulations fit a specified number of components. |
| **Shift-invariant trilinearity and PARAFAC2×N** | Shift-invariant trilinearity uses Fourier-domain information to accommodate elution shifts. PARAFAC2×N extends drift modelling to multiple separation modes in multidimensional chromatography. [Shift-invariant trilinearity](https://doi.org/10.1002/cem.3501) · [PARAFAC2×N](https://doi.org/10.1016/j.aca.2023.340909) | Each model is fitted at a chosen component count. |

The 1,000 shared factors matched the 1,000 generating compounds
one-to-one across all eight occurrences of each compound. This is exact
automatic recovery of the chemical component count from 550 measured ion
channels. Scan order distinguishes compounds that elute at different
times, and each compound contributes a shared spectrum wherever it
occurs.

Across the 160,000 scans, the median number of simultaneously eluting
true peaks was one, the 95th percentile was three and the maximum was
five, counting a peak wherever it reached 1% of its own apex. Thus 1,000
measures the number of distinct compounds across the chromatograms. An
unrestricted nonnegative factorisation of the stacked scan-by-ion matrix
has nonnegative rank at most 550, since the 550 ion columns themselves
provide a nonnegative basis.

Published rank-selection experiments put the scale in perspective.
Distributed NMF (DnMFk) recovered seven planted factors in a
one-terabyte matrix. Eswar et al. compared NMF rank selectors at true
rank 10. Fu et al. tested nonnegative tensor rank selection at true
ranks 15–40. The 1,000 compounds here are **25 times the largest
generating rank in that tensor experiment**. These studies use different
data models and assess different kinds of rank. [Chennupati et
al.](https://doi.org/10.1007/s11227-020-03181-6) · [Eswar et
al.](https://www.osti.gov/servlets/purl/3002147) · [Fu et
al.](https://doi.org/10.3390/e21070632)
