\name{GCMSbenchmark}
\alias{GCMSbenchmark}
\alias{scoreGCMSbenchmark}
\title{Run and score GC/MS deconvolution workflows}
\usage{
GCMSbenchmark(cdf_file, output_dir, method = "GCMSprime", truth = NULL,
  ncores = 1L, control = list(), timeout = 7200, rt_tolerance = 0.12,
  min_cosine = 0.70, matching = c("assignment", "reciprocal"),
  ri_tolerance = 30, reference = NULL, ion_threshold = 0.001, ...)

scoreGCMSbenchmark(x, truth = NULL, rt_tolerance = 0.12,
  min_cosine = 0.70, matching = c("assignment", "reciprocal"),
  ri_tolerance = 30, reference = NULL, ion_threshold = 0.001)
}
\arguments{
\item{cdf_file}{Local ANDI-MS CDF paths with unique basenames, or a directory of CDF files.}
\item{output_dir}{A fresh directory for native exports, standardized RDS files and tables.}
\item{method}{One or more of \code{"GCMSprime"}, \code{"MS-DIAL"}, \code{"mzmine"},
\code{"eRah"}, \code{"AMDIS"}, \code{"xcms+CAMERA"} and \code{"MSHub"}. Use \code{"all"} to run the available
adapters for the input structure and operating system. MSHub requires at least
ten samples. AMDIS execution requires Windows. Existing ELU exports can be imported
on any platform.}
\item{truth}{Optional list keyed by input basename without extension. Each record contains
\code{id}, \code{time} in minutes, \code{mass}, \code{H} (scan by compound),
\code{P} (compound by mass), \code{apex} in minutes, \code{area} in integrated
counts and a logical \code{acid} vector. Supply \code{reference_apex} if known for
comparison with methods that return aligned retention times. For RI scoring,
add \code{ri_calibration}, an increasing data frame with \code{time} (minutes)
and \code{ri}. Supply \code{reference_ri_calibration} for aligned RT outputs.}
\item{ncores}{Requested maximum number of local threads or workers. AMDIS is serial.}
\item{control}{Named list of backend-specific argument lists. For external software,
set the executable path in \code{exe} if it is absent from PATH. The
README lists the other settings and software requirements. For
\code{"xcms+CAMERA"}, \code{detection$fwhm} sets the xcms filter width in
seconds, \code{fwhm_grouping} sets CAMERA's retention-time grouping
arguments, and \code{grouping = FALSE} uses FWHM grouping alone.}
\item{timeout}{Timeout in seconds for each external command.}
\item{rt_tolerance}{Maximum absolute retention-time discrepancy in minutes for a match.}
\item{min_cosine}{Minimum cosine similarity between nominal-mass spectra.}
\item{matching}{\code{"assignment"} maximises qualifying one-to-one matches within
the RT window. \code{"reciprocal"} requires mutually best spectral matches within
the RI window using the supplied reference library.}
\item{ri_tolerance}{Maximum absolute RI difference for reciprocal matching.}
\item{reference}{For reciprocal matching, a list with unique compound \code{id},
compound-by-mass matrix \code{P}, nominal \code{mass} values and compound
retention indices \code{ri}. Include library compounds absent from particular
samples. If omitted when rescoring, the object's stored reference is reused.}
\item{ion_threshold}{Ion-support cutoff as a fraction of each spectrum's base
peak, default 0.001 (0.1\%). The same rule is applied to reference and estimated
spectra for every method. Zero counts every positive ion as present.
This cutoff also sets the GCMSprime reporting option.}
\item{...}{Extra arguments passed to \code{GCMSprime::gcmsDeco()} only.}
\item{x}{A \code{GCMSbenchmark} result.}
}
\description{
Run local deconvolution backends on common input files, retain their native
exports, and compare their extracted components with optional simulated truth.
}
\details{
Peak matching within each sample first maximises the number of qualifying
one-to-one pairs, then maximises the sum of their spectral cosines.
With \code{matching="reciprocal"}, spectra are restricted by RI and the cosine
cutoff before choosing the best recovered spectrum for each library entry and
the best library entry for each recovered spectrum. A pair is retained when
both choices agree. Equal cosines are ordered by RI difference and then ID.
The matching step sees the full reference library. A reciprocal hit to a
compound absent from that sample remains a false positive.

\code{identifications} records each called component's best eligible library
entry and whether it is reciprocal and present. Extra calls to the same entry
are counted in \code{duplicate_calls}. All unmatched calls count as false
positives. Areas use the recorded area of the selected component.
Profile cosine uses linearly interpolated intensities on the original
time grid. Area errors are calculated among matched peaks with positive
areas. Native exports are retained, while standardized areas use total-ion count
equivalents and spectra have unit L1 norm. Model-ion areas are divided by the
reported model-ion fraction; mzmine minute integrals are converted to seconds,
and time integrals are converted using the acquisition scan interval. eRah's
exported factor scale is restored. CAMERA profiles remain explicitly labelled
TIC-scaled representative-ion traces. These conversions use native exports and
acquisition metadata only.

\code{area_adjusted_R2} describes the robust log--log regression among matched
positive pairs using \code{robustbase::lmrob(setting="KS2014")} and its corrected
adjusted R-squared. \code{area_OLS_adjusted_R2} retains ordinary least-squares
R-squared. Prediction bands use the robust residual scale and coefficient
covariance, with approximately Gaussian errors for the main log-area
distribution. Every observed point remains in the plot.
\code{area_pairs} records the
number of pairs used. These conditional area metrics should be read alongside
TP, FP, FN and the number of matched pairs.

\code{spectral_F1} is the arithmetic mean of within-spectrum ion-support F1
across correctly matched true occurrences. An ion is present when its positive
intensity is at least \code{ion_threshold} times its spectrum's base peak.
Per-spectrum TP, FP and FN count
ions on the union of observed and reference mass grids, with
\eqn{F_1=2TP/(2TP+FP+FN)}. Matching and spectral cosine use the reported
intensities without this additional support cutoff.
\code{dispersion_estimate} records the median fitted
quasi-Poisson dispersion when the engine supplies it.

The adapters retain native areas, conversion factors, alignment groups and
source-output labels alongside standardized results. Data conversion,
analysis and adapter elapsed times are recorded separately where available.
}
\value{
A \code{GCMSbenchmark} object containing \code{methods}, \code{files},
\code{directory}, \code{summary} and \code{matches}. A successful method
contains a peak table with retention times in minutes and lists of spectra
and profiles in peak-table order. A failed backend records its status and
message. \code{scoreGCMSbenchmark()} adds or updates the scoring tables and,
for reciprocal matching, the identification table.
}
\examples{
\dontrun{
b <- GCMSbenchmark("samples", "comparison",
  method = c("GCMSprime", "eRah"), ncores = 8L)
print(b)
}
}
