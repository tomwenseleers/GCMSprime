GCMS_report_spectra <- function(S,ion_threshold=.001,normalize=FALSE) {
  if(!is.numeric(ion_threshold)||length(ion_threshold)!=1L||
      !is.finite(ion_threshold)||ion_threshold<0||ion_threshold>1)
    stop("ion_threshold must be a base-peak fraction between 0 and 1.",call.=FALSE)
  if(!nrow(S))return(S)
  total<-rowSums(S)
  S[S < ion_threshold*apply(S,1L,max)]<-0
  retained<-rowSums(S)
  S*ifelse(retained>0,if(normalize)1/retained else total/retained,0)
}

GCMS_ion_threshold <- function(object) {
  if(is.null(object$ion_threshold)) .001 else object$ion_threshold
}

GCMS_xlsx_path <- function(path) {
  if(!is.character(path)||length(path)!=1L||is.na(path)||!nzchar(path)||
      !grepl("[.]xlsx$",path,ignore.case=TRUE)) stop("xlsx_file must end in .xlsx.",call.=FALSE)
  if(!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx for Excel export.",call.=FALSE)
  invisible(path)
}

GCMS_spectrum_text <- function(S,mass,ion_threshold=.001) {
  S<-GCMS_report_spectra(S,ion_threshold)
  vapply(seq_len(nrow(S)),function(k) {
    s<-S[k,];intensity<-if(max(s)>0)as.integer(round(9999*s/max(s)))else integer(length(s))
    keep<-which(intensity>0);keep<-keep[order(-intensity[keep],mass[keep])]
    paste0(mass[keep]," (",intensity[keep],")",collapse=", ")
  },"")
}

GCMS_write_xlsx <- function(model,time,mass,compounds,path,cdf_file,alkanes_file,seconds,matches=NULL,calibration=NULL,search_settings=NULL,ion_threshold=.001) {
  GCMS_xlsx_path(path)
  wb<-openxlsx::createWorkbook(creator="GCMSprime")
  header<-openxlsx::createStyle(fgFill="#102d4d",fontColour="#FFFFFF",textDecoration="bold",
    fontSize=11,wrapText=TRUE,valign="center")
  compounds$Spectrum<-GCMS_spectrum_text(model$S[compounds$compound,,drop=FALSE],mass,ion_threshold)
  if(any(nchar(compounds$Spectrum)>32767L))
    stop("A spectrum exceeds Excel's cell text limit; use the MSP or R spectra.",call.=FALSE)
  first<-c("compound","retention_time","area","Spectrum")
  compounds<-compounds[,c(first,setdiff(names(compounds),first)),drop=FALSE]
  metadata<-data.frame(field=c("sample","source CDF","calibration CDF","noise family",
    "dispersion","theta","scans","mass channels","compounds","retention time unit",
      "area unit","area_min unit","spectral display","ion threshold (base-peak fraction)","fitting seconds"),
    value=c(if(nrow(compounds)) compounds$sample[1L] else tools::file_path_sans_ext(basename(cdf_file)),
      cdf_file,if(is.null(alkanes_file)) "" else alkanes_file,model$spec$name,
      format(model$spec$dispersion),format(model$spec$theta),length(time),length(mass),nrow(compounds),
      "minutes","intensity-scans","intensity-minutes","m/z (intensity), decreasing intensity; base peak 9999; shared relative-ion cutoff followed by display rounding.",format(ion_threshold),
      if("total"%in%names(seconds)) seconds[["total"]] else sum(seconds)),
    stringsAsFactors=FALSE)
  if(!is.null(search_settings)) metadata<-rbind(metadata,data.frame(
    field=c("RI phase class","RI tolerance","unlabelled reference RI phase","spectral methods",
      "reference count convention","RI timing convention","small-peak width source"),
    value=c(search_settings$phase,search_settings$ri_tolerance,
      if(is.null(search_settings$unknown_phase)) "excluded from constrained search" else
        paste("explicitly assumed",search_settings$unknown_phase),
      paste(search_settings$method,collapse=", "),search_settings$reference_count,
      if(is.null(calibration)) "explicit query RI" else calibration$timing,
      if(is.null(calibration$width_source)) "not supplied" else calibration$width_source)))
  tables<-list(Compounds=compounds,Run=metadata)
  if(!is.null(matches)) tables$Matches<-matches
  if(inherits(calibration,"RIcalibration")) tables$RI_calibration<-calibration$anchors
  for(name in names(tables)) {
    openxlsx::addWorksheet(wb,name,gridLines=FALSE)
    openxlsx::writeData(wb,name,tables[[name]],headerStyle=header,withFilter=name!="Run")
    openxlsx::freezePane(wb,name,firstActiveRow=2L,firstActiveCol=if(name=="Compounds") 4L else 1L)
    openxlsx::setRowHeights(wb,name,1L,32)
    openxlsx::setColWidths(wb,name,seq_len(ncol(tables[[name]])),if(name=="Run") c(28,90) else 16)
  }
  openxlsx::setColWidths(wb,"Compounds",4L,160)
  if(nrow(compounds)) {
    rows<-seq_len(nrow(compounds))+1L
    openxlsx::addStyle(wb,"Compounds",openxlsx::createStyle(numFmt="0.00000"),rows=rows,cols=2L,gridExpand=TRUE)
    openxlsx::addStyle(wb,"Compounds",openxlsx::createStyle(numFmt="0.000E+00"),rows=rows,
      cols=match(c("area","area_min"),names(compounds)),gridExpand=TRUE)
    openxlsx::addStyle(wb,"Compounds",openxlsx::createStyle(wrapText=FALSE,valign="center"),
      rows=rows,cols=4L,gridExpand=TRUE)
    openxlsx::setRowHeights(wb,"Compounds",rows,22)
  }
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  openxlsx::saveWorkbook(wb,path,overwrite=TRUE)
  invisible(path)
}
