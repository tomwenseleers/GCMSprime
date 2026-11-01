# Plot GC/MS signals, compound profiles and recovered spectra

Draw four signal heatmaps, background-corrected observed mass traces
with fitted compound profiles, or recovered mass spectra. Observations
are read from the source CDF for the requested interval and released
after plotting.

## Usage

``` r
# S3 method for class 'GCMSdeco'
plot(
  x,
  type = c("lines", "heatmap", "spectra", "profiles"),
  file = NULL,
  time_range = NULL,
  mass = NULL,
  which = NULL,
  nrows = 4L,
  ncols = 2L,
  minutes_per_panel = 2.5,
  transform = c("sqrt", "identity", "log1p"),
  peak_scale = c("basepeak", "total"),
  labels = FALSE,
  max_time = 1600L,
  clip = 0.999,
  gamma = 0.13,
  residual_limit = 3,
  order = c("rt", "area", "height", "index"),
  ...
)
```

## Arguments

- x:

  A fitted
  [`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md)
  object.

- type:

  `"lines"`, `"heatmap"`, `"spectra"` or `"profiles"`. The last places
  each compound's total-ion profile beside its spectrum, with one
  compound per row.

- file:

  Optional PDF output path. `NULL` uses the current device. PDF reports
  can contain multiple pages.

- time_range:

  Retention-time bounds in minutes; `NULL` uses the full fitted
  interval.

- mass:

  Optional vector of recorded nominal m/z values.

- which:

  Compound indices, as returned in `x$compounds$compound`. `NULL` shows
  all compounds.

- nrows, ncols:

  Rows per PDF page, and columns for spectral panels.
  Profile-and-spectrum panels always use two columns.

- minutes_per_panel:

  Retention-time span of each line panel, default 2.5 minutes.
  Consecutive panels cover the whole selected interval.

- transform:

  Line-display transform: signed square root (default), identity or
  signed `log1p`. Negative background-corrected traces are retained.
  Tick labels use the original intensity scale.

- peak_scale:

  `"basepeak"` scales each profile by its strongest ion in the displayed
  mass range, matching the individual grey mass traces. `"total"` uses
  the compound's summed ion intensity and shows the observed total ion
  trace in grey.

- labels:

  Label line-profile apices by compound index.

- max_time:

  Maximum heatmap time bins. Data and fitted signals are averaged over
  the same bins, retaining all selected mass channels.

- clip:

  Upper quantile used for heatmap colour saturation, default .999. Set
  to one to display the full transformed range.

- gamma:

  Positive heatmap contrast exponent, default .13 as in
  [`plotmat()`](https://tomwenseleers.github.io/L0glm/reference/plotmat.html).

- residual_limit:

  Symmetric Pearson-residual colour limit, default 3.

- order:

  Spectral panel order: `"rt"` (default) uses increasing apex retention
  time; `"area"` or `"height"` uses decreasing compound area or
  total-ion peak height; `"index"` preserves the selected indices.

- ...:

  Graphical arguments passed to the line or spectral panel.

## Value

Invisibly, `x`. Plotting does not attach observed data to it.

## Details

Heatmaps show the observed signal, observed signal minus modeled
background, fitted chromatogram and Pearson residuals. Spline background
and any known offset are subtracted in the second heatmap and grey line
traces. Signal heatmaps use the
[`plotmat()`](https://tomwenseleers.github.io/L0glm/reference/plotmat.html)
black-blue-green-yellow-orange-red palette and truncate negative
corrected intensities at zero. Pearson residuals use blue-black-red for
negative-zero-positive deviations, without an intensity transform.
Display time bins contain their average residual. Heatmaps are
interpolated onto a regular display grid and drawn as rasters; red
vertical lines mark the selected compound apices. Spectral panels use
base-peak percentages. The Excel compound sheet places a compact
spectrum beside the integrated area; `coef(x,normalize="l1")` returns
unit-L1 spectra.

## See also

[`GCMSdeco()`](https://tomwenseleers.github.io/GCMSprime/reference/GCMSdeco.md),
[`L0nmf::splogconNMF()`](https://tomwenseleers.github.io/L0nmf/reference/splogconNMF.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- GCMSdeco("beer.cdf", "beer.msp", xlsx_file="beer.xlsx")
plot(fit, type="heatmap", file="beer_heatmaps.pdf")
plot(fit, type="lines", file="beer_traces.pdf", minutes_per_panel=2.5)
plot(fit, type="spectra", file="beer_spectra.pdf",
     which=head(fit$compounds$compound[order(-fit$compounds$area)],24))
} # }
```
