#' Annotate recovered GC/MS compounds and export library matches
#'
#' Add retention indices and ranked library matches to an existing deconvolution,
#' without refitting its profiles or spectra.
#' @param x A GCMSdeco object.
#' @param library A GCMSlibrary.
#' @param ... Search options passed to library_search(), including calibration,
#'   phase, ri_tolerance, method and n_matches.
#' @param msp_file Optional new annotated MSP output path.
#' @param xlsx_file Optional Excel output path with compound, match and calibration sheets.
#' @return The GCMSdeco object with identification (a GCMSmatches object), retention
#'   indices and best-match columns in compounds. The first requested method supplies
#'   the best-match columns; all requested methods remain in identification$matches.
#' @details The original compound names, factorization and fitted ion totals are
#'   preserved. Library names are tentative annotations. MSP files include the
#'   RI, column class, original fitted ion total and best-match description.
#'   Excel includes every returned match and the accepted alkane calibration anchors.
#'   Explicit assignments of unlabelled reference RI phases are recorded in the
#'   compound table, match table, Run sheet and MSP annotation.
#'   Supply output paths to export companion files from a saved fit.
#' @seealso [calibrate_ri()], [library_search()], [GCMSdeco()]
#' @examples
#' # Given a saved deconvolution and an imported spectral library:
#' \dontrun{
#' cal <- calibrate_ri("alkanes.cdf", carbon=7:23)
#' identified <- identify_compounds(fit, lib, calibration=cal, phase="DB5",
#'   ri_tolerance=30, n_matches=5, method=c("stein_scott","hellinger"),
#'   msp_file="sample_identified.msp", xlsx_file="sample_identified.xlsx")
#' identified$identification$matches
#' }
#' @export
identify_compounds <- function(x,library,...,msp_file=NULL,xlsx_file=NULL) {
  if(!inherits(x,"GCMSdeco")) stop("x must be a GCMSdeco fit.")
  begun<-proc.time()[['elapsed']]
  found<-library_search(x,library,...)
  q<-found$queries;c<-x$compounds;at<-match(c$compound,q$compound)
  c$retention_coordinate<-q$retention_coordinate[at];c$retention_index<-q$query_RI[at]
  c$ri_phase<-rep(found$settings$phase,nrow(c));c$identification_status<-q$status[at]
  hits<-found$matches
  best<-hits[hits$method==found$settings$method[1L]&hits$rank==1,,drop=FALSE]
  at<-match(c$compound,best$compound)
  for(field in c('name','cas','library_id','score','delta_RI','method','source','origin','ri_origin','ri_phase','ri_assumed')) {
    value<-best[[field]]
    c[[paste0('match_',field)]]<-if(is.null(value)) rep(NA,nrow(c)) else value[at]
  }
  x$compounds<-c;x$identification<-found
  x$seconds['library_search']<-found$seconds
  if(!is.null(msp_file)) {
    GCMS_write_msp(x$model,x$time,x$mass,msp_file,c)
    x$msp_file<-normalizePath(msp_file,winslash='/')
  }
  if(!is.null(xlsx_file)) {
    GCMS_write_xlsx(x$model,x$time,x$mass,c,xlsx_file,x$cdf_file,x$alkanes_file,x$seconds,
      matches=found$matches,calibration=found$calibration,search_settings=found$settings)
    x$xlsx_file<-normalizePath(xlsx_file,winslash='/')
  }
  x$seconds['total']<-x$seconds['total']+proc.time()[['elapsed']]-begun
  x
}
