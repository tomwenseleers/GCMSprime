#' Import a GCMS-ID predicted spectral library
#'
#' @param file Extracted GCMS-ID CSV file from its public downloads page.
#' @param backend,db_file,batch_size Storage options as in [readMSP()].
#' @param max_spectra Maximum number of CSV records to import.
#' @return A GCMSlibrary object with predicted spectra and predicted RI provenance.
#' @details The GCMS-ID columns compound_name, inchikey, compound_formula,
#'   retention_index_value, exact_mass and peak_list are read in batches. The
#'   StdNP, SemiStdNP and StdPolar name suffixes identify DB1, DB5 and DBWAX
#'   phase classes. Records with different RI predictions are retained separately.
#'   Unknown suffixes retain an unknown phase. All source columns are preserved.
#' @seealso <https://gcms-id.ca/downloads>
#' @export
readGCMSID <- function(file,backend=c("auto","memory","sqlite"),db_file=NULL,
                         batch_size=1000L,max_spectra=Inf) {
  backend<-match.arg(backend)
  if(!file.exists(file)) stop("CSV file does not exist.")
  if(backend=="auto") {
    backend<-if(file.info(file)$size>100e6) "sqlite" else "memory"
    if(backend=="sqlite"&&is.null(db_file)) db_file<-paste0(file,".sqlite")
  }
  input<-file(file,"rt",encoding="UTF-8");on.exit(close(input))
  header<-names(utils::read.csv(text=readLines(input,n=1L),check.names=FALSE))
  needed<-c("compound_name","inchikey","compound_formula","retention_index_value","exact_mass","peak_list")
  if(!all(needed%in%header)) stop("The GCMS-ID CSV schema is not recognized.")
  if(length(batch_size)!=1L||!is.finite(batch_size)||batch_size<1||batch_size!=floor(batch_size)) stop("Use a positive integer batch_size.")
  if(length(max_spectra)!=1L||is.na(max_spectra)||max_spectra<1||max_spectra!=floor(max_spectra)) stop("Use a positive max_spectra.")
  tmp<-tempfile(fileext=".msp");on.exit(unlink(tmp),add=TRUE)
  out<-file(tmp,"wt",encoding="UTF-8");on.exit(try(close(out),silent=TRUE),add=TRUE)
  n<-0L
  repeat {
    z<-tryCatch(utils::read.csv(input,header=FALSE,col.names=header,nrows=min(batch_size,max_spectra-n),
      check.names=FALSE,stringsAsFactors=FALSE,colClasses="character",na.strings=character()),
      error=function(e) if(grepl("no lines available",conditionMessage(e),fixed=TRUE)) data.frame() else stop(e))
    if(!nrow(z)) break
    for(i in seq_len(nrow(z))) {
      phase<-if(grepl("_SemiStdNP$",z$compound_name[i])) "DB5" else
        if(grepl("_StdNP$",z$compound_name[i])) "DB1" else
        if(grepl("_StdPolar$",z$compound_name[i])) "DBWAX" else "unknown"
      peaks<-strsplit(z$peak_list[i],";",fixed=TRUE)[[1]];peaks<-trimws(peaks);peaks<-peaks[nzchar(peaks)]
      fields<-c(paste0("Name: ",z$compound_name[i]),paste0("InChIKey: ",z$inchikey[i]),
        paste0("Formula: ",z$compound_formula[i]),paste0("RI_",phase,": ",z$retention_index_value[i]),
        paste0("Exact Mass: ",z$exact_mass[i]),paste0("Source ID: ",n+i),
        paste0("GCMS-ID_",header,": ",as.character(z[i,])))
      writeLines(c(fields,paste0("Num Peaks: ",length(peaks)),paste0(peaks,";"),""),out)
    }
    n<-n+nrow(z);if(n>=max_spectra) break
  }
  close(out)
  readMSP(tmp,backend=backend,db_file=db_file,batch_size=batch_size,
    source=paste0("GCMS-ID: ",basename(file)),origin="predicted")
}
