#' Add alkane-calibrated retention indices without a library search
#'
#' Annotate a single-sample or multisample deconvolution with retention indices
#' and optionally export annotated MSP and Excel files. Profiles and spectra
#' are retained unchanged.
#' @param x A [GCMSdeco()] or [GCMSbatch] result.
#' @param calibration A [calibrate_ri()] result. `NULL` uses each sample's
#'   `alkanes_file` and calibrated PSF.
#' @param phase Column class, default `"DB5"`.
#' @param extrapolate Linearly extend the RI curve beyond its alkane anchors.
#'   The default returns `NA` outside the anchor range. Extended values are
#'   flagged in the compound table and MSP comments.
#' @param msp_file,xlsx_file Optional output paths. For multisample results these
#'   name the consensus outputs; the corresponding per-sample files are also
#'   refreshed. Updated sample RDS objects remain in the batch output directory.
#' @return The annotated object, including the calibration, retention coordinates,
#'   indices, column classes and extrapolation flags. Multisample consensus RIs
#'   are medians over detected occurrences.
#' @details Peak coordinates use the 1%-height onset plus half the calibrated
#'   narrow-peak width, measured at 1% height in the local window. A natural
#'   cubic spline interpolates the corresponding alkane coordinates and indices.
#'   MSP exports use NIST's phase-labelled `Retention_index`
#'   syntax: `SemiStdNP` for DB5, `StdNP` for DB1 and `StdPolar` for DBWAX.
#' @references NIST MS Search 2.3 manual, pp. 53--54, RI import/export:
#'   <https://chemdata.nist.gov/dokuwiki/lib/exe/fetch.php?media=chemdata%3Anist17%3Anistms_ver23man.pdf>
#' @examples
#' \dontrun{
#' cal <- calibrate_ri("alkanes.cdf", psf=standard_psf)
#' fit <- annotate_ri(fit, cal, phase="DB5",
#'   msp_file="sample_RI.msp", xlsx_file="sample_RI.xlsx")
#' fit$compounds[c("retention_time", "retention_index", "ri_extrapolated")]
#' }
#' @seealso [calibrate_ri()], [identify_compounds()]
#' @md
#' @export
annotate_ri <- function(x,calibration=NULL,phase=c("DB5","DB1","DBWAX"),
                        extrapolate=FALSE,msp_file=NULL,xlsx_file=NULL) {
  phase<-match.arg(phase)
  if(!is.logical(extrapolate)||length(extrapolate)!=1L||is.na(extrapolate))stop("extrapolate must be TRUE or FALSE.")
  begun<-proc.time()[['elapsed']]
  if(inherits(x,"GCMSbatch")) {
    tables<-calibrations<-vector('list',length(x$paths))
    for(s in seq_along(x$paths)) {
      fit<-readRDS(x$paths[s])
      fit<-annotate_ri(fit,calibration,phase,extrapolate,
        if(is.null(msp_file))NULL else fit$msp_file,
        if(is.null(xlsx_file))NULL else fit$xlsx_file)
      saveRDS(fit,x$paths[s]);tables[[s]]<-fit$compounds;calibrations[[s]]<-fit$ri_calibration
    }
    x$ri_calibration<-if(is.null(calibration))calibrations else calibration
    x$samples$RI_available<-vapply(tables,function(t)sum(is.finite(t$retention_index)),0L)
    if(!is.null(x$coupled)) {
      all<-do.call(rbind,tables);G<-nrow(x$coupled$P)
      value<-vapply(seq_len(G),function(g) {
        z<-all$retention_index[all$coupled_compound==g];z<-z[is.finite(z)]
        if(length(z))median(z)else NA_real_
      },0.)
      extra<-vapply(seq_len(G),function(g)any(all$ri_extrapolated[all$coupled_compound==g]),FALSE)
      x$retention_indices<-data.frame(compound=seq_len(G),retention_index=value,
        ri_phase=phase,ri_extrapolated=extra)
      if(!is.null(msp_file)||!is.null(xlsx_file))GCMS_batch_consensus(x,msp_file,xlsx_file)
    }
    utils::write.csv(x$samples,file.path(x$output_dir,'samples.csv'),row.names=FALSE)
    x$seconds['retention_indices']<-proc.time()[['elapsed']]-begun
    x$seconds['total']<-x$seconds['total']+x$seconds['retention_indices']
    saveRDS(x,file.path(x$output_dir,'experiment.rds'));return(x)
  }
  if(!inherits(x,'GCMSdeco'))stop('x must be a GCMSdeco or GCMSbatch result.')
  query<-search_queries(x,calibration,NULL,NULL,extrapolate)
  q<-query$queries;calibration<-query$calibration;c<-x$compounds;at<-match(c$compound,q$compound)
  c$peak_onset<-q$peak_onset[at];c$narrow_width<-q$narrow_width[at]
  c$retention_coordinate<-q$retention_coordinate[at]
  c$retention_index<-q$query_RI[at];c$ri_phase<-rep(phase,nrow(c))
  bounds<-range(calibration$anchors$retention_time)
  c$ri_extrapolated<-is.finite(c$retention_index)&
    (c$retention_coordinate<bounds[1L]|c$retention_coordinate>bounds[2L])
  x$compounds<-c;x$ri_calibration<-calibration
  if(!is.null(msp_file)) {
    GCMS_write_msp(x$model,x$time,x$mass,msp_file,c)
    x$msp_file<-normalizePath(msp_file,winslash='/')
  }
  if(!is.null(xlsx_file)) {
    GCMS_write_xlsx(x$model,x$time,x$mass,c,xlsx_file,x$cdf_file,x$alkanes_file,x$seconds,
      matches=if(is.null(x$identification))NULL else x$identification$matches,calibration=calibration)
    x$xlsx_file<-normalizePath(xlsx_file,winslash='/')
  }
  x$seconds['retention_indices']<-proc.time()[['elapsed']]-begun
  x$seconds['total']<-x$seconds['total']+x$seconds['retention_indices'];x
}
