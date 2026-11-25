# Benchmarks against other software

## Small-scale 16-sample software comparison

This page reports the **small-scale** seven-workflow comparison. It is
distinct from the [large-scale 16-sample resource
benchmark](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html#benchmark-scales)
with 10,000 scans per sample and 1,000 reference compounds across the
experiment. The linked overview reports both timings, settings and short
reproduction recipes. Here GCMSprime uses 32 common ion groups, two
blind-PSF windows, at most three PSF iterations, and up to 12 initial
conditional sweeps on the original ions.

The **16 chromatograms** share a library of **100 reference compounds**.
Each spans a 20-minute retention window, with 1,200 scans and all 258
mass channels. Compounds are absent from about half the samples, giving
**769 true occurrences**. Library retention indices determine their
positions. One compound has a broad shark-fin-shaped acid profile; the
remaining profiles are exponentially modified Gaussians. Higher
concentrations extend peaks towards later times from approximately fixed
leading edges, and a small positive drift shifts successive injections
to the right.

The smooth B-spline background is derived from the experimental
`gcmsData2` window and contributes 3.03% of the integrated compound
signal. The simulation retains its native abundance scale and uses
quasi-Poisson noise with **dispersion 783.18**, estimated from that
window. GCMSprime estimates dispersion separately in each sample, with a
**median estimate of 829.55**. The generating retention-index maps are
used for scoring; the deconvolution receives neither the known areas nor
the true profiles.

The optional search recovered **true rank 100 exactly (estimated rank
100)** and **769 of 769 occurrences with zero false positives**. The
acid was recovered in all seven samples where it was present.

| Method        | Seconds |  Estimated rank |  TP |   FP | Sensitivity | Peak F1 | Ion F1 |
|:--------------|--------:|----------------:|----:|-----:|------------:|--------:|-------:|
| GCMSprime     |  196.28 | **100 (exact)** | 769 |    0 |      100.0% |   1.000 |  0.995 |
| MS-DIAL       |   17.61 |               — | 675 |  641 |       87.8% |   0.647 |  0.858 |
| mzmine GC-EI  |   38.47 |               — | 622 |   70 |       80.9% |   0.851 |  0.860 |
| eRah          |   17.07 |               — | 645 | 1912 |       83.9% |   0.388 |  0.771 |
| AMDIS         |   92.58 |               — | 671 | 1153 |       87.3% |   0.518 |  0.854 |
| xcms + CAMERA |   24.06 |               — | 635 |  146 |       82.6% |   0.819 |  0.617 |
| MSHub         |   20.58 |             281 | 635 | 2797 |       82.6% |   0.302 |  0.844 |

The true rank is 100 distinct compounds present in at least one sample.
The rank column counts shared fitted factors. A dash identifies
workflows without a comparable cross-sample factorisation rank.

Times are wall-clock analysis times on an HP ZBook Fury G8 Mobile
Workstation with an Intel Core i9-11950H at 2.60 GHz, eight physical
cores and 16 logical threads, Windows and R 4.5.0. AMDIS runs serially
and xcms + CAMERA uses four workers. Input conversion is timed
separately. The six comparator rows use their saved outputs and timings
from the identical simulated inputs.

A shorter search recovered **true rank 100 exactly** in **104.48 s**,
with all 769 true occurrences and only 1 extra occurrence. The same
control recovered **true rank 1,000 exactly** on the independent large
example. [Compare local search at both
scales](https://tomwenseleers.github.io/GCMSprime/articles/multisample.html#local-rank-search).

The small experiment used 100 extra selective sweeps, a temporary
fourfold penalty and a second blind PSF. Preparation, PSF calibration
and ion grouping took 14.59 s, grouped regression and joint fitting
177.10 s, and export 4.33 s. The enclosing benchmark took **196.28 s**
from CDF input through results. The large resource benchmark used the
default support settings.

With the default stopping and support settings, the same CDF files took
**54.66 s** and gave **estimated rank 115** against **true rank 100**,
with 769 true positives and only 26 false positives. The
software-comparison figures use the 196.28 s extended fit. The
original-channel GIC of the selected rank-100 fit was **4,146,698**.

GCMSprime used four two-thread workers for grouped PSF regressions and
could use eight threads for the joint fit. Profile and background
updates run across samples. Shared spectra run across ions. Local merge
decisions and file exports also contribute to elapsed time.

## Sensitivity to the starting peak width

The same 16 chromatograms were also fitted from Gaussian PSFs of 3.64
and 14.56 scans at half height. The calibrated and both Gaussian starts
recovered **estimated rank 100, equal to true rank 100**, all 769
occurrences and zero false positives after the blind PSF exchange.
Across starts, the largest difference between two estimates of a matched
peak area was 0.86%. The acid area in sample 009 was 99.76–99.81% of
truth. The calibration, support proposals and GIC comparisons used only
observed counts.

![False positive and missed calls across three PSF starts and three
searches](../reference/figures/psf-start-recovery.png)

False positive and missed calls across three PSF starts and three
searches

[Rank, original-price GIC, recovery and area measurements for all nine
fits](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/psf-starts.csv).

## Reciprocal spectral matches and integrated areas

The reference library contains all 100 compounds, including those absent
from an individual sample. Candidate spectra must fall within **30 RI
units** and have cosine similarity of at least **0.70**. A reciprocal
pair selects each other as the best eligible spectral hit. Ties are
resolved by RI difference and then ID. A reciprocal hit to a compound
actually present is a true positive; every other returned component is a
false positive. Additional calls to the same reference compound are
included in those counts. Their areas are not added to the selected
match.

[![Ordinary log–log peak-area fits and 95% prediction
intervals](../reference/figures/simulated-16-area-predictions.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-area-predictions.png)

All panels use **total-ion count equivalents**. Spectra are normalized
to unit L1 norm, with their magnitude transferred to the corresponding
profile or area. Native model-ion areas are corrected for their ion
fraction and integration-time units. The conversions use each engine’s
exported factors and acquisition metadata, without calibration against
the known answers. The native values and conversion factors remain in
the result objects.

For a component with profile h and spectrum s, let c=\sum_j s_j, p=s/c
and \widetilde h=ch. Then hs^\top=\widetilde h p^\top and \sum_j p_j=1.
An area measured for a particular ion q becomes
A\_{\mathrm{TIC}}=A_q/p_q, after converting its integration-time units.

The solid line fits an ordinary least-squares regression of log₁₀
estimated area on log₁₀ true area to all positive matched pairs,
including the outliers. The shaded region is the pointwise 95%
prediction interval for an individual peak at a given true area. The
dashed line shows exact agreement. Darker points indicate overlap
(opacity 0.5).

GCMSprime recovered **all 769 occurrences**, with **adjusted R²
0.999993** and a **median absolute area error of 0.079%**. Even the
lowest estimate was 98.3% of its true area. Area agreement is
conditional on recovered true matches, so it is shown together with
sensitivity and false-positive counts.

[![Reciprocal spectral cosine similarity with individual
points](../reference/figures/simulated-16-spectral-violin.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-spectral-violin.png)

Each point is a reciprocal spectral hit. Crosses identify hits to
compounds absent from that sample; the cosine axis uses logit spacing.
Display values at or above 0.999999 are shown at 1, while the statistics
retain their original values. Compact horizontal tick labels keep the
scale readable.

## Correctness of spectral ion support

[![Within-spectrum ion-support F1 with the shared 0.1%
cutoff](../reference/figures/simulated-16-spectral-support.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-spectral-support.png)

Ion support uses a shared **0.1% of base-peak intensity cutoff**,
applied independently to reference and recovered spectra for every
method: I_j\>0 and I_j\geq0.001\max_l I_l. Within-spectrum
F_1=2TP/(2TP+FP+FN) is averaged over correctly matched true spectra.
GCMSprime achieves **mean ion F1 0.995**, compared with **0.860** for
the strongest comparator. Individual points accompany the violins; the
logit-spaced axis retains labels in the original F1 units.

This is also GCMSprime’s default spectral reporting cutoff for
coefficient summaries, library queries and MSP/Excel output.
`ion_threshold=0` retains all positive ions. Reporting preserves fitted
total ion counts, peak areas, profiles and information criteria.

## Observed signal and recovered peak areas

[![Observed and true signal with square-root intensity
colours](../reference/figures/simulated-16-sample-heatmap.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-sample-heatmap.png)

The observed counts and true mean, including background, use the plotMat
palette with **square-root intensity colours**. Red lines mark true
apices.

[![Method-by-compound peak areas with yellow false
negatives](../reference/figures/simulated-16-area-heatmap.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-area-heatmap.png)

The ground-truth row and method rows show one sample’s true or inferred
total-ion areas. **Bright yellow marks false negatives**; grey marks no
signal or call. Compounds follow increasing RI. A [pooled
view](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-area-heatmap-pooled.png)
summarizes all 16 samples; the [sample-by-sample
PDF](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-area-heatmaps-all-samples.pdf)
shows every acquisition.

## Correctness of inferred elution profiles

[![Solid inferred and dashed true regular-compound profiles, blue to red
by true
area](../reference/figures/simulated-16-regular-profiles.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-regular-profiles.png)

[![Solid inferred and dashed true acid-compound profiles, blue to red by
true
area](../reference/figures/simulated-16-acid-profiles.png)](https://tomwenseleers.github.io/GCMSprime/reference/figures/simulated-16-acid-profiles.png)

The first panel set shows the regular compound with the largest total
true area across samples (present in **11 of 16 samples**); the second
shows the acid compound (present in **7 of 16**). These are the same
saved profiles shown on the logconNMF homepage. All 16 samples are
included and ordered by true abundance, from **blue for the smallest
area to red for the largest**. **Dashed curves are ground truth; solid
curves are returned profiles.** Intensities retain their original
total-ion scale, without unit-area or peak-height normalization. Absent
compounds have zero true signal. Compatible fragments remain separate
curves so splitting stays visible. The generating retention-index
calibration supplies the common horizontal axis; profile amplitudes
remain on their native total-ion scale.

GCMSprime recovers all occurrences of both the regular EMG and the
shark-fin-shaped acid. Mean profile cosine, counting missed occurrences
as zero, is **0.999996 for the regular compound** and **0.999983 for the
acid**; median cosine among matched occurrences is **0.999998** and
**0.999990**, respectively. Other configurations often truncate
profiles, split their contributions or miss the acid. MS-DIAL also
resolves the regular compound (mean cosine 0.998). It recovers five of
the seven acid occurrences and has acid mean cosine 0.449. CAMERA
supplies a representative-ion trace scaled to the grouped total-ion
area. The mzmine export here contains spectra and areas without
component scan profiles. MSHub factorizes aligned ion-peak integrals
across samples. Its model has no common compound elution profile
constrained to be shared across mass channels, so component scan
profiles are absent from the model itself. See [Aksenov et
al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7971188/).

## Detection settings and reproducibility

The xcms + CAMERA configuration uses an **8.6-second filter FWHM**,
signal-to-noise threshold **1.5**, a minimum of **four ions**, and
CAMERA FWHM grouping with `perfwhm=0.8`. Settings were selected on
odd-numbered samples and checked on even-numbered samples. This
configuration recovered 635 of 769 occurrences across all samples (82.6%
sensitivity).

### Recovery across true peak sizes

Sensitivity is also calculated separately in quartiles of true
integrated peak area, including every true occurrence. This makes the
preferential recovery of large peaks visible alongside the conditional
area regressions. GCMSprime recovered **all 193 peaks in the smallest
quarter**.

| Method        | Smallest quarter | Lower-middle | Upper-middle | Largest quarter |
|:--------------|-----------------:|-------------:|-------------:|----------------:|
| GCMSprime     |           100.0% |       100.0% |       100.0% |          100.0% |
| MS-DIAL       |            79.3% |        84.4% |        91.1% |           96.4% |
| mzmine GC-EI  |            70.5% |        79.2% |        82.3% |           91.7% |
| eRah          |            72.0% |        78.6% |        89.1% |           95.8% |
| AMDIS         |            78.2% |        83.3% |        89.6% |           97.9% |
| xcms + CAMERA |            66.3% |        81.8% |        84.9% |           97.4% |
| MSHub         |            73.6% |        82.8% |        86.5% |           87.5% |

GCMSbenchmark includes the seeded simulation, empirical-background
template, engine settings, scoring tables and plotting code. Its
[installation
guide](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/INSTALL.md)
covers the R packages and external engines.

``` r

library(GCMSbenchmark)
bench <- system.file("benchmarks", package="GCMSbenchmark", mustWork=TRUE)
source(file.path(bench, "simulated-16.R"), encoding="UTF-8")
source(file.path(bench, "plot-simulated-16.R"), encoding="UTF-8")

# Recreate the plots from the included results.
plotBenchmark16(file.path(bench, "results"), "benchmark-figures")

# Generate the same samples and rerun selected installed workflows.
fit <- runBenchmark16("benchmark-run",
  method=c("GCMSprime", "xcms+CAMERA"), ncores=16L, support_search=TRUE)
```

The simulation uses seed 9118. The tested versions are MS-DIAL
5.5.261105, mzmine 4.11.0, eRah 2.2.0, AMDIS 2.71, xcms 4.6.4, CAMERA
1.64.0 and standalone MSHub revision `59fe5d5d7ad5`. Native engine
exports are retained. The MS-DIAL method and mzmine batch are included.
AMDIS uses component width 16 scans, sensitivity 60, resolution 0 and
peak-shape setting 0.

Download the [summary
table](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/benchmark.csv),
[matched
pairs](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/matches.csv),
[all component
identifications](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/identifications.csv),
[stage timings and coupled
GIC](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/native-timings.csv),
and [xcms tuning
results](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/xcms-training-results.csv).

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
ion traces, correlations, local covariance, factorization of peak
integrals or PARAFAC2 least squares. Their preprocessing noise estimates
and selection rules are summarized below. [Mathematical
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
| PARADISe | Fits PARAFAC2 within retention intervals. Intervals and ranks can require review, and the constant cross-product constraint restricts how peak overlap can change across samples. [Johnsen et al. (2017)](https://doi.org/10.1016/j.chroma.2017.04.052) · [Protocol](https://doi.org/10.21203/rs.3.pex-2143/v1) | Free Community and paid Professional editions. [Kromath](https://www.kromath.com/) |
