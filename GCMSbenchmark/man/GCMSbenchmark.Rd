\name{GCMSbenchmark}
\alias{GCMSbenchmark}
\alias{scoreGCMSbenchmark}
\title{Run and score GC/MS deconvolution workflows}
\usage{
GCMSbenchmark(cdf_file, output_dir, method = "GCMSprime", truth = NULL,
  ncores = 1L, control = list(), timeout = 7200, rt_tolerance = 0.12,
  min_cosine = 0.70, ...)

scoreGCMSbenchmark(x, truth = NULL, rt_tolerance = 0.12,
  min_cosine = 0.70)
}
\arguments{
\item{cdf_file}{Local ANDI-MS CDF paths with unique basenames, or a directory of CDF files.}
\item{output_dir}{A fresh directory for native exports, standardized RDS files and tables.}
\item{method}{One or more of \code{"GCMSprime"}, \code{"MS-DIAL"}, \code{"mzmine"},
\code{"eRah"}, \code{"AMDIS"} and \code{"MSHub"}. Use \code{"all"} to run the available
adapters for the input structure. MSHub requires at least ten samples.}
\item{truth}{Optional list keyed by input basename without extension. Each record contains
\code{id}, \code{time} in minutes, \code{mass}, \code{H} (scan by compound),
\code{P} (compound by mass), \code{apex} in minutes, \code{area} in integrated
counts and a logical \code{acid} vector. Supply \code{reference_apex} if known for
comparison with methods that return aligned retention times.}
\item{ncores}{Requested maximum number of local threads or workers. AMDIS is serial.}
\item{control}{Named list of backend-specific argument lists. For external software,
set the executable path in \code{exe} if it is absent from PATH. The
README lists the other settings and software requirements.}
\item{timeout}{Timeout in seconds for each external command.}
\item{rt_tolerance}{Maximum absolute retention-time discrepancy in minutes for a match.}
\item{min_cosine}{Minimum cosine similarity between nominal-mass spectra.}
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
Profile cosine uses linearly interpolated intensities on the original
time grid. Area errors are calculated among matched peaks with positive
areas. Direct count-scale errors are available when the backend area has
the same units as truth. Scale-adjusted errors centre the log area ratios.

The adapters preserve differences in native output, including model-ion
and component TIC profiles, alignment groups and area units. Data conversion,
analysis and adapter elapsed times are recorded separately where available.
}
\value{
A \code{GCMSbenchmark} object containing \code{methods}, \code{files},
\code{directory}, \code{summary} and \code{matches}. A successful method
contains a peak table with retention times in minutes and lists of spectra
and profiles in peak-table order. A failed backend records its status and
message. \code{scoreGCMSbenchmark()} adds or updates the two scoring tables.
}
\examples{
\dontrun{
b <- GCMSbenchmark("samples", "comparison",
  method = c("GCMSprime", "eRah"), ncores = 8L)
print(b)
}
}
