library_require <- function() {
  for(p in c("DBI","RSQLite","jsonlite")) if(!requireNamespace(p,quietly=TRUE))
    stop("Install ",p," for spectral-library storage.",call.=FALSE)
}
library_connect <- function(x) {
  library_require()
  con<-DBI::dbConnect(RSQLite::SQLite(),x$path,flags=RSQLite::SQLITE_RO)
  version<-DBI::dbGetQuery(con,"SELECT version FROM format")$version
  if(!identical(version,1L)) {DBI::dbDisconnect(con);stop("Unsupported spectral-library format.")}
  con
}
library_decode <- function(x) readBin(x,double(),n=length(x)%/%8L,size=8L,endian="little")
library_encode <- function(x) writeBin(as.double(x),raw(),size=8L,endian="little")
library_text <- function(x,encoding) {
  value<-iconv(x,from=encoding,to="UTF-8")
  missing<-is.na(value)&!is.na(x)
  if(any(missing)) value[missing]<-iconv(x[missing],from="windows-1252",to="UTF-8",sub="byte")
  value
}
library_record <- function(x,id,source,encoding,origin) {
  tags<-library_text(x$tags,encoding);values<-library_text(x$values,encoding)
  raw_fields<-data.frame(tag=tags,value=values)
  comments<-values[toupper(tags)%in%c("COMMENTS","COMMENT")]
  pairs<-unlist(lapply(comments,function(value)
    regmatches(value,gregexpr('"[^"=]+=[^"]*"',value,perl=TRUE))[[1]]),use.names=FALSE)
  if(length(pairs)) {
    pairs<-substring(pairs,2L,nchar(pairs)-1L)
    at<-regexpr("=",pairs,fixed=TRUE)
    tags<-c(tags,substring(pairs,1L,at-1L));values<-c(values,substring(pairs,at+1L))
  }
  # OpenChrom preserves NIST's column-specific RI summary in comments.
  ri_pairs<-unlist(lapply(comments,function(value) regmatches(value,
    gregexpr('(?i)\\b(?:SemiStdNP|StdNP|StdPolar)\\s*=\\s*[0-9][^[:space:];"]*',value,perl=TRUE))[[1]]),use.names=FALSE)
  if(length(ri_pairs)) {
    at<-regexpr("=",ri_pairs,fixed=TRUE)
    tags<-c(tags,trimws(substring(ri_pairs,1L,at-1L)))
    values<-c(values,trimws(substring(ri_pairs,at+1L)))
  }
  keys<-toupper(gsub("[^[:alnum:]]","",tags))
  field<-function(names) {i<-which(keys%in%names);if(length(i)) values[i[1]] else NA_character_}
  metadata<-data.frame(id=as.integer(id),name=field("NAME"),formula=field("FORMULA"),
    cas=field(c("CASNO","CAS","CASNUMBER")),inchi=field(c("INCHI","STDINCHI")),
    inchikey=field(c("INCHIKEY","STDINCHIKEY")),smiles=field(c("SMILES","CANONICALSMILES")),
    license=field("LICENSE"),source=source,source_id=field(c("DB","ID","REFID","DBNUMBER","ACCESSION","SOURCEID")),
    origin=origin,source_line=x$line,npeaks=length(x$mass),stringsAsFactors=FALSE)
  count<-suppressWarnings(as.numeric(field(c("PSEUDOCOUNT","IONTOTAL","TOTALIONCOUNT"))))
  metadata$pseudocount<-if(length(count)==1L&&is.finite(count)&&count>0) count else NA_real_
  ri<-which(keys%in%c("RI","RETENTIONINDEX","KOVATSRI","KOVATSINDEX","RIUNKNOWN","SEMISTDNP","STDNP","STDPOLAR")|
    grepl("^(RETENTIONINDEX|KOVATS|RI(DB|STANDARD|SEMI))",keys))
  ri<-ri[!duplicated(paste(keys[ri],trimws(values[ri]),sep="\r"))]
  indices<-lapply(ri,function(i) {
    number<-suppressWarnings(as.numeric(sub("^.*?([0-9]+(?:\\.[0-9]+)?).*$","\\1",values[i],perl=TRUE)))
    if(!is.finite(number)||number<=0) return(NULL)
    column<-field(c("COLUMN","COLUMNNAME","COLUMNPHASE","STATIONARYPHASE","GCCOLUMN"))
    label<-paste(tags[i],values[i]);compact<-toupper(gsub("[^[:alnum:]]","",label))
    phase<-if(grepl("DBWAX|STANDARDPOLAR|STDPOLAR|PEG|CARBOWAX",compact)) "DBWAX" else
      if(grepl("DB5|SEMISTANDARDNONPOLAR|SEMISTDNP|5MS",compact)) "DB5" else
      if(grepl("DB1|STANDARDNONPOLAR|STDNP",compact)) "DB1" else "unknown"
    if(phase=="unknown"&&!is.na(column)) {
      compact<-toupper(gsub("[^[:alnum:]]","",column))
      phase<-if(grepl("DBWAX|STANDARDPOLAR|STDPOLAR|PEG|CARBOWAX",compact)) "DBWAX" else
        if(grepl("DB5|SEMISTANDARDNONPOLAR|SEMISTDNP|5MS",compact)) "DB5" else
        if(grepl("DB1|STANDARDNONPOLAR|STDNP",compact)) "DB1" else "unknown"
    }
    data.frame(id=as.integer(id),value=number,phase_class=phase,
      column=if(is.na(column)) tags[i] else paste(tags[i],column,sep="; "),
      convention=if(grepl("KOVATS",compact)) "Kovats" else "unknown",
      origin=if(grepl("PREDICT",compact)||origin=="predicted") "predicted" else "unknown",
      source=source,raw=values[i],stringsAsFactors=FALSE)
  })
  list(metadata=metadata,fields=raw_fields,
    indices=do.call(rbind,indices),mass=x$mass,intensity=x$intensity)
}

#' Import a mass spectral library
#'
#' Stream NIST-compatible MSP files into sparse in-memory or disk-backed storage.
#' OpenChrom Agilent and Shimadzu exports and MoNA/Golm MSP fields are accepted.
#'
#' @param file Local MSP file. Compressed files should first be extracted.
#' @param backend Storage: auto uses SQLite for files larger than 100 MB.
#' @param db_file SQLite output path; required for the sqlite backend. Auto storage
#'   uses a new file beside the MSP when this is omitted.
#' @param batch_size Number of records buffered per transaction.
#' @param source Library provenance, defaulting to the input filename.
#' @param encoding Text encoding; UTF-8 input with Windows-1252 fallback is supported.
#' @param origin Whether spectra are experimental, predicted or of unknown origin.
#' @param strict Check the declared number of peaks against the parsed pairs.
#' @param max_spectra Optional import limit, useful for previews.
#' @return A GCMSlibrary object. Memory objects contain sparse spectra;
#'   disk objects retain a path and open read-only connections when queried.
#' @details MSP tags are case-insensitive. Repeated tags, synonyms and unknown
#'   fields remain available through library_metadata(raw=TRUE). Raw masses and
#'   intensities are retained as double precision arrays. Sparse matrix extraction
#'   sums duplicate masses; nominal binning and normalization are explicit options.
#'   Explicit Pseudocount/Ion_total tags retain a supplied full-spectrum total;
#'   library_search uses it for count-sensitive comparisons without changing raw peaks.
#'   In-memory libraries store a dgCMatrix plus compact metadata. Exceptional
#'   unsorted, duplicate or zero-intensity peak records retain exact raw arrays.
#'
#'   SQLite stores one pair of little-endian float64 blobs per spectrum, indexed
#'   metadata and a long retention-index table. Import is transactional and uses
#'   bounded batches. Library files can be moved and reopened with open_library().
#'   Existing databases are never overwritten. No live connection is serialized.
#'
#'   Retention indices retain their raw field, phase class and provenance.
#'   NIST/OpenChrom comment summaries SemiStdNP, StdNP and StdPolar map to
#'   DB5, DB1 and DBWAX; slash-delimited uncertainty fields remain in the raw value.
#'   DB1, DB5 and DBWAX denote nonpolar, semi-standard nonpolar and polar classes;
#'   a separate column/stationary-phase field also supplies the RI phase class.
#'   An unspecified column remains unknown. Zero RI placeholders are omitted from
#'   the index table. Supplied SMILES are preserved without chemical canonicalization.
#' @examples
#' f <- tempfile(fileext=".msp")
#' writeLines(c("Name: example", "Formula: C2H6O", "RI_DB5: 560",
#'              "Num Peaks: 3", "31 100; 45 60; 46 20;"), f)
#' lib <- read_msp(f)
#' lib
#' library_metadata(lib)
#' retention_indices(lib)
#' spectra(lib)
#' @export
read_msp <- function(file,backend=c("auto","memory","sqlite"),db_file=NULL,
                     batch_size=1000L,source=basename(file),encoding="UTF-8",
                     origin=c("unknown","experimental","predicted"),strict=TRUE,
                     max_spectra=Inf) {
  backend<-match.arg(backend);origin<-match.arg(origin)
  if(!file.exists(file)||dir.exists(file)) stop("MSP file does not exist.")
  if(length(batch_size)!=1L||!is.finite(batch_size)||batch_size<1||batch_size!=floor(batch_size)) stop("Use a positive integer batch_size.")
  if(length(max_spectra)!=1L||is.na(max_spectra)||max_spectra<1||max_spectra!=floor(max_spectra)) stop("Use a positive max_spectra.")
  if(backend=="auto") {
    backend<-if(file.info(file)$size>100e6) "sqlite" else "memory"
    if(backend=="sqlite"&&is.null(db_file)) db_file<-paste0(file,".sqlite")
  }
  handle<-msp_open_cpp(normalizePath(file,winslash="/",mustWork=TRUE));on.exit(msp_close_cpp(handle),add=TRUE);n<-0L;records<-list()
  if(backend=="sqlite") {
    library_require()
    if(is.null(db_file)||length(db_file)!=1L||file.exists(db_file)) stop("Supply a new db_file path.")
    dir.create(dirname(db_file),recursive=TRUE,showWarnings=FALSE)
    temporary<-tempfile(tmpdir=dirname(db_file),fileext=".sqlite")
    con<-DBI::dbConnect(RSQLite::SQLite(),temporary)
    on.exit({DBI::dbDisconnect(con);if(file.exists(temporary)) unlink(temporary)},add=TRUE)
    DBI::dbWriteTable(con,"format",data.frame(version=1L,source=source))
    DBI::dbExecute(con,"PRAGMA journal_mode=DELETE")
  }
  repeat {
    batch<-msp_read_cpp(handle,as.integer(min(batch_size,max_spectra-n)),isTRUE(strict))
    if(!length(batch)) break
    parsed<-lapply(seq_along(batch),function(i) library_record(batch[[i]],n+i,source,encoding,origin))
    if(backend=="memory") records<-c(records,parsed) else {
      meta<-do.call(rbind,lapply(parsed,`[[`,"metadata"))
      meta$fields<-vapply(parsed,function(x) jsonlite::toJSON(x$fields,auto_unbox=TRUE),"")
      meta$mass<-lapply(parsed,function(x) library_encode(x$mass))
      meta$intensity<-lapply(parsed,function(x) library_encode(x$intensity))
      ri<-do.call(rbind,lapply(parsed,`[[`,"indices"))
      DBI::dbWithTransaction(con,{
        if(DBI::dbExistsTable(con,"spectra")) DBI::dbAppendTable(con,"spectra",meta) else DBI::dbWriteTable(con,"spectra",meta)
        if(!is.null(ri)) {if(DBI::dbExistsTable(con,"retention_indices")) DBI::dbAppendTable(con,"retention_indices",ri) else DBI::dbWriteTable(con,"retention_indices",ri)}
      })
    }
    n<-n+length(parsed);if(n>=max_spectra) break
  }
  if(!n) stop("No spectra imported.")
  if(backend=="sqlite") {
    library_require()
    DBI::dbExecute(con,"CREATE UNIQUE INDEX spectrum_id ON spectra(id)")
    for(field in c("name","cas","inchikey","source_id")) DBI::dbExecute(con,paste0("CREATE INDEX spectrum_",field," ON spectra(",field,")"))
    if(DBI::dbExistsTable(con,"retention_indices")) {
      DBI::dbExecute(con,"CREATE INDEX retention_id ON retention_indices(id)")
      DBI::dbExecute(con,"CREATE INDEX retention_phase_value ON retention_indices(phase_class,value,id)")
    }
    # Close before renaming on Windows; the on.exit handler is then replaced.
    DBI::dbDisconnect(con);on.exit({msp_close_cpp(handle);if(file.exists(temporary)) unlink(temporary)},add=FALSE)
    if(!file.rename(temporary,db_file)) stop("Cannot finalize library database.")
    return(open_library(db_file))
  }
  metadata<-do.call(rbind,lapply(records,`[[`,"metadata"));rownames(metadata)<-NULL
  fields<-lapply(records,function(x) stats::setNames(x$fields$value,x$fields$tag))
  mass<-lapply(records,`[[`,"mass");intensity<-lapply(records,`[[`,"intensity")
  grid<-sort(unique(unlist(mass)))
  matrix<-Matrix::drop0(Matrix::sparseMatrix(i=match(unlist(mass),grid),
    j=rep(seq_len(n),lengths(mass)),x=unlist(intensity),dims=c(length(grid),n),
    dimnames=list(as.character(grid),as.character(seq_len(n)))))
  exceptions<-lapply(seq_len(n),function(i) {
    if(is.unsorted(mass[[i]],strictly=TRUE)||any(intensity[[i]]==0))
      cbind(mass=mass[[i]],intensity=intensity[[i]]) else NULL
  })
  structure(list(backend="memory",n=n,matrix=matrix,mass=grid,metadata=metadata,
    fields=fields,exceptions=exceptions,indices=do.call(rbind,lapply(records,`[[`,"indices")),
    source=source),class="GCMSlibrary")
}

#' @rdname read_msp
#' @param path Existing spectral-library SQLite file.
#' @export
open_library <- function(path) {
  if(!file.exists(path)) stop("Library database does not exist.")
  x<-structure(list(backend="sqlite",path=normalizePath(path,winslash="/")),class="GCMSlibrary")
  con<-library_connect(x);on.exit(DBI::dbDisconnect(con))
  x$n<-DBI::dbGetQuery(con,"SELECT COUNT(*) AS n FROM spectra")$n
  x$source<-DBI::dbGetQuery(con,"SELECT source FROM format")$source;x
}
library_ids <- function(x,which) {
  if(is.null(which)) return(seq_len(x$n))
  if(!is.numeric(which)||anyNA(which)||any(which<1|which>x$n|which!=floor(which))) stop("Invalid spectrum indices.")
  as.integer(which)
}
library_fetch <- function(x,which,peaks=FALSE) {
  ids<-library_ids(x,which)
  if(x$backend=="memory") return(lapply(ids,function(i) {
    raw<-NULL
    if(peaks) {
      raw<-x$exceptions[[i]]
      if(is.null(raw)) {
        column<-x$matrix[,i,drop=FALSE]
        raw<-cbind(mass=x$mass[column@i+1L],intensity=column@x)
      }
    }
    list(metadata=x$metadata[i,,drop=FALSE],
      fields=data.frame(tag=names(x$fields[[i]]),value=unname(x$fields[[i]])),
      mass=if(peaks) unname(raw[,1]) else NULL,intensity=if(peaks) unname(raw[,2]) else NULL)
  }))
  if(!length(ids)) return(list())
  con<-library_connect(x);on.exit(DBI::dbDisconnect(con))
  columns<-if(peaks) "*" else paste(DBI::dbQuoteIdentifier(con,
    setdiff(DBI::dbListFields(con,"spectra"),c("mass","intensity"))),collapse=",")
  z<-DBI::dbGetQuery(con,paste0("SELECT ",columns," FROM spectra WHERE id IN (",paste(unique(ids),collapse=","),")"))
  z<-z[match(ids,z$id),,drop=FALSE]
  lapply(seq_along(ids),function(i) list(metadata=z[i,setdiff(names(z),c("mass","intensity","fields")),drop=FALSE],
    fields=jsonlite::fromJSON(z$fields[i]),mass=if(peaks) library_decode(z$mass[[i]]) else NULL,
    intensity=if(peaks) library_decode(z$intensity[[i]]) else NULL))
}

#' Inspect spectral-library metadata
#' @param x A GCMSlibrary object.
#' @param which Spectrum indices; NULL selects every record.
#' @param raw Return a list of original tag/value tables, preserving repeated tags.
#' @return A metadata data frame, or a list of raw field tables.
#' @export
library_metadata <- function(x,which=NULL,raw=FALSE) {
  z<-library_fetch(x,which)
  if(raw) return(lapply(z,`[[`,"fields"))
  do.call(rbind,lapply(z,`[[`,"metadata"))
}

#' Extract spectral-library peaks or a sparse mass matrix
#' @param x A GCMSlibrary object.
#' @param which Spectrum indices; select batches for large databases.
#' @param format Sparse matrix (mass rows, spectrum columns) or raw peak lists.
#' @param nominal Bin masses to the nearest integer; FALSE preserves exact masses.
#' @param normalize Optional unit-sum (l1) or unit-Euclidean-norm (l2) spectra.
#' @return A dgCMatrix or a list of two-column mass/intensity matrices.
#' @export
spectra <- function(x,which=NULL,format=c("matrix","list"),nominal=FALSE,
                    normalize=c("none","l1","l2")) {
  format<-match.arg(format);normalize<-match.arg(normalize);ids<-library_ids(x,which)
  z<-library_fetch(x,ids,TRUE)
  raw<-lapply(z,function(x)cbind(mass=x$mass,intensity=x$intensity))
  if(format=="list") return(raw)
  masses<-lapply(raw,function(x) if(nominal) floor(x[,1]+.5) else x[,1]);grid<-sort(unique(unlist(masses)))
  ans<-Matrix::sparseMatrix(i=match(unlist(masses),grid),j=rep(seq_along(z),lengths(masses)),
    x=unlist(lapply(raw,function(x)x[,2])),dims=c(length(grid),length(z)),
    dimnames=list(as.character(grid),as.character(ids)))
  ans<-Matrix::drop0(ans)
  if(normalize!="none") {
    d<-if(normalize=="l1") Matrix::colSums(ans) else sqrt(Matrix::colSums(ans^2))
    ans<-ans%*%Matrix::Diagonal(x=ifelse(d>0,1/d,0));colnames(ans)<-as.character(ids)
  }
  ans
}

#' Extract retention indices and their provenance
#' @param x A GCMSlibrary object.
#' @param which Spectrum indices.
#' @param phase Optional DB1, DB5, DBWAX or unknown phase filter.
#' @return A long data frame with one row per recorded index.
#' @export
retention_indices <- function(x,which=NULL,phase=NULL) {
  ids<-library_ids(x,which)
  if(x$backend=="memory") {
    z<-x$indices
    if(!is.null(z)) z<-z[z$id%in%ids,,drop=FALSE]
  } else {
    con<-library_connect(x);on.exit(DBI::dbDisconnect(con))
    z<-if(DBI::dbExistsTable(con,"retention_indices")&&length(ids)) DBI::dbGetQuery(con,
      paste0("SELECT * FROM retention_indices WHERE id IN (",paste(unique(ids),collapse=","),")")) else NULL
  }
  if(is.null(z)) return(data.frame(id=integer(),value=numeric(),phase_class=character(),
    column=character(),convention=character(),origin=character(),source=character(),raw=character()))
  if(!is.null(phase)) z<-z[z$phase_class%in%phase,,drop=FALSE]
  z[order(match(z$id,ids)),,drop=FALSE]
}

#' @rdname read_msp
#' @param x A GCMSlibrary object.
#' @param ... Additional printing arguments.
#' @method print GCMSlibrary
#' @export
print.GCMSlibrary <- function(x,...) {
  cat(format(x$n,big.mark=",")," spectra; ",x$backend," library\n",sep="")
  if(x$backend=="sqlite") cat(x$path,"\n")
  invisible(x)
}
