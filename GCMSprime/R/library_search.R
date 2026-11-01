search_methods <- c("stein_scott","cosine","hellinger","jensen_shannon","multinomial")

search_queries <- function(x,calibration,ri,pseudocount,extrapolate,require_ri=TRUE) {
  fit<-inherits(x,"GCMSdeco")
  S<-if(fit) as.matrix(x$model$S[x$compounds$compound,,drop=FALSE]) else as.matrix(x)
  mass<-if(fit) x$mass else suppressWarnings(as.numeric(colnames(S)))
  if(!is.numeric(S)||!ncol(S)||any(!is.finite(S)|S<0)||length(mass)!=ncol(S)||
      any(!is.finite(mass)|mass<=0)) stop("Use a nonnegative compound-by-mass matrix with numeric mass names.")
  # Nominal bins match ANDI-MS import; retain all acquired bins, including zeros.
  mass<-round(mass);colnames(S)<-as.character(mass)
  if(anyDuplicated(mass)) S<-t(rowsum(t(S),mass,reorder=TRUE))
  mass<-sort(unique(mass));S<-S[,order(as.numeric(colnames(S))),drop=FALSE]
  colnames(S)<-as.character(mass)
  total<-rowSums(S)
  if(any(!is.finite(total))) stop("Spectrum totals must be finite.")
  if(!is.null(pseudocount)) {
    if(!is.numeric(pseudocount)||!length(pseudocount)%in%c(1L,nrow(S))||
        any(!is.finite(pseudocount)|pseudocount<=0)||any(total<=0)) stop("pseudocount must be positive for every nonempty spectrum.")
    S<-S*rep(pseudocount,length.out=nrow(S))/total;total<-rowSums(S)
  }
  q<-data.frame(query=seq_len(nrow(S)),compound=if(fit) x$compounds$compound else seq_len(nrow(S)),
    name=if(fit) x$compounds$name else if(is.null(rownames(S))) paste0("spectrum_",seq_len(nrow(S))) else rownames(S),
    query_total=as.numeric(total),retention_coordinate=rep(NA_real_,nrow(S)),query_RI=rep(NA_real_,nrow(S)))
  if(!nrow(S)) return(list(S=S,mass=mass,queries=q,calibration=calibration))
  if(!is.null(ri)) {
    if(!is.numeric(ri)||length(ri)!=nrow(S)||any(!is.na(ri)&!is.finite(ri))) stop("ri must contain one RI per query spectrum.")
    q$query_RI<-ri
  } else if(fit) {
    if(is.null(calibration)&&!is.null(x$alkanes_file)) calibration<-calibrate_ri(x$alkanes_file)
    if(!require_ri&&(!inherits(calibration,"RIcalibration")||is.null(calibration$width_curve)))
      return(list(S=S,mass=mass,queries=q,calibration=calibration))
    if(!inherits(calibration,"RIcalibration")) stop("Supply an RI calibration, alkane file in the fitted object, or explicit ri values.")
    if(is.null(calibration$width_curve)) stop("Calibration requires narrow_width values for timing profiles; alternatively supply query ri.")
    q$peak_onset<-vapply(seq_len(ncol(x$model$H)),function(j)
      ri_crossings(x$time,x$model$H[,x$compounds$compound[j]],which.max(x$model$H[,x$compounds$compound[j]]),calibration$fraction)[['start']],0.)
    q$narrow_width<-ri_predict_width(calibration,q$peak_onset)
    q$retention_coordinate<-q$peak_onset+q$narrow_width/2
    q$query_RI<-as.numeric(predict(calibration,q$retention_coordinate,extrapolate=extrapolate))
  }
  list(S=S,mass=mass,queries=q,calibration=calibration)
}

search_ri <- function(library,con,phase,unknown_phase,range) {
  if(is.null(con)) r<-library$indices else if(DBI::dbExistsTable(con,"retention_indices")) {
    phases<-c(phase,if(identical(unknown_phase,phase)) "unknown")
    sql<-paste0("SELECT * FROM retention_indices WHERE phase_class IN (",
      paste(DBI::dbQuoteString(con,phases),collapse=","),")")
    if(all(is.finite(range))) sql<-paste0(sql," AND value BETWEEN ",range[1]," AND ",range[2])
    r<-DBI::dbGetQuery(con,sql)
  } else r<-NULL
  if(is.null(r)||!nrow(r)) return(data.frame(id=integer(),value=numeric(),phase_class=character()))
  r<-r[r$phase_class==phase|r$phase_class=="unknown"&identical(unknown_phase,phase),,drop=FALSE]
  r[r$value>=range[1]&r$value<=range[2],,drop=FALSE]
}

search_matrix <- function(library,con,ids,mass) {
  if(is.null(con)) {
    a<-library$matrix[,ids,drop=FALSE];raw_total<-as.numeric(Matrix::colSums(a))
    count<-library$metadata$pseudocount[match(ids,library$metadata$id)]
    i<-match(round(library$mass),mass);j<-rep.int(seq_along(ids),diff(a@p))
    valid<-!is.na(i[a@i+1L]);ii<-i[a@i[valid]+1L];jj<-j[valid];xx<-a@x[valid]
  } else {
    fields<-if("pseudocount"%in%DBI::dbListFields(con,"spectra")) ",pseudocount" else ""
    a<-DBI::dbGetQuery(con,paste0("SELECT id,mass,intensity",fields," FROM spectra WHERE id IN (",paste(ids,collapse=","),") ORDER BY id"))
    a<-a[match(ids,a$id),,drop=FALSE];count<-a$pseudocount
    mz<-lapply(a$mass,library_decode);intensity<-lapply(a$intensity,library_decode)
    raw_total<-vapply(intensity,sum,0.);i<-match(round(unlist(mz,use.names=FALSE)),mass)
    j<-rep.int(seq_along(ids),lengths(mz));x<-unlist(intensity,use.names=FALSE)
    valid<-!is.na(i);ii<-i[valid];jj<-j[valid];xx<-x[valid]
  }
  if(length(count)) {
    scale<-ifelse(is.finite(count)&count>0&raw_total>0,count/pmax(raw_total,.Machine$double.xmin),1)
    xx<-xx*scale[jj];raw_total<-raw_total*scale
  }
  list(matrix=Matrix::sparseMatrix(i=ii,j=jj,x=xx,dims=c(length(mass),length(ids))),raw_total=raw_total)
}

#' Identify deconvoluted spectra using retention indices and spectral libraries
#'
#' Search acquired nominal-mass bins with fast native spectral comparisons and
#' column-specific retention-index filtering. Fitted ion totals are preserved.
#'
#' @param x A GCMSdeco fit, or a compound-by-mass matrix with numeric mass names.
#' @param library A GCMSlibrary returned by read_msp(), read_gcms_id() or open_library().
#' @param calibration RIcalibration from calibrate_ri(). For GCMSdeco fits, the
#'   stored alkane standard is used automatically when available.
#' @param ri Optional RI vector, one value per query; overrides profile-based timing.
#' @param phase RI column class: DB1 (nonpolar), DB5 (semi-standard nonpolar),
#'   or DBWAX (polar).
#' @param ri_tolerance Maximum absolute RI difference, default 30. Inf requests
#'   an unconstrained spectral comparison, also allowing references without RI.
#' @param n_matches Maximum number of matches per spectrum and scoring method.
#' @param method One or more of stein_scott, cosine, hellinger, jensen_shannon,
#'   multinomial, or all. The default is the Stein--Scott weighted dot product.
#' @param pseudocount Optional total for each query (or one common total), useful
#'   when passing L1-normalized spectra. GCMSdeco already retains fitted ion totals.
#' @param extrapolate Allow linear RI extrapolation beyond the calibrated anchors.
#' @param unknown_phase Optional explicit column-class assignment for library
#'   RI records whose phase is unknown. NULL excludes these from constrained searches.
#' @param ncores Native search threads; NULL uses the physical core count.
#' @param batch_size Maximum library spectra loaded per batch; bounds SQLite memory.
#' @return A GCMSmatches list with matches (ranked long table), queries (including
#'   unassigned spectra), calibration, search settings and elapsed seconds.
#'   as.data.frame() extracts the match table. Matches are tentative annotations.
#' @details RI filtering precedes spectral scoring. Only the chosen phase class
#'   is accepted within ri_tolerance; missing query RI produces no constrained
#'   matches. Multiple library RI measurements are compared using the closest.
#'   A compound with fewer eligible references returns fewer than n_matches.
#'   The same onset-plus-half-narrow-width coordinate is used for standards and
#'   fitted profiles. Unknown library column classes require an explicit assignment.
#'
#'   Let a and b be the acquired-bin intensities, N=sum(a), M=sum(b),
#'   p=a/N and q=b/M. stein_scott is the squared cosine between
#'   a^0.6 * mass^3 and b^0.6 * mass^3; cosine uses unweighted intensities.
#'   hellinger uses sqrt(1-sum(sqrt(p*q))). jensen_shannon uses
#'   sqrt((KL(p,r)+KL(q,r))/(2*log(2))) with r=(p+q)/2.
#'   multinomial returns the symmetric two-sample deviance
#'   2*(N*KL(p,u)+M*KL(q,u)), u=(a+b)/(N+M), retaining both totals.
#'   Zero contributions use 0*log(0)=0. These distances are finite with disjoint support.
#'
#'   Scores run from 0 to 1000, with 1000 denoting identical spectral shape:
#'   squared cosine, Bhattacharyya coefficient, 1-normalized-JS, and 1/(1+deviance),
#'   respectively. Compare scores within a method. The returned distance is the
#'   quantity minimized in ranking. Relative library intensities are pseudocounts;
#'   their arbitrary scale affects multinomial deviance. An explicit library
#'   Pseudocount/Ion_total tag supplies the original total when available. Hellinger and equal-weight
#'   Jensen--Shannon compare shapes independently of that scale. query_total,
#'   library_total and library_mass_coverage make the count and acquisition conventions explicit.
#'
#'   References are restricted to mass bins actually acquired for the query.
#'   Spectra are streamed from SQLite in bounded batches. Native OpenMP workers
#'   compare queries independently with deterministic ranking and no nested fitting.
#' @references Stein, S. E. and Scott, D. R. (1994). J. Am. Soc. Mass Spectrom.
#'   5, 859--866. doi:10.1016/1044-0305(94)87009-8.
#'   Kim, S. et al. (2012). Bioinformatics 28, 1158--1163.
#'   doi:10.1093/bioinformatics/bts083.
#' @examples
#' f <- tempfile(fileext=".msp")
#' writeLines(c("Name: compound A", "RI_DB5: 1000", "Num Peaks: 3",
#'              "43 100; 57 80; 71 20;", "Name: compound B", "RI_DB5: 1200",
#'              "Num Peaks: 2", "43 20; 91 100;"), f)
#' lib <- read_msp(f)
#' q <- matrix(c(500,400,100,0),1,dimnames=list(NULL,c(43,57,71,91)))
#' hit <- library_search(q,lib,ri=1005,phase="DB5",method="all",ncores=1)
#' as.data.frame(hit)
#' @export
library_search <- function(x,library,calibration=NULL,ri=NULL,phase=c("DB5","DB1","DBWAX"),
                           ri_tolerance=30,n_matches=5L,method="stein_scott",pseudocount=NULL,
                           extrapolate=FALSE,unknown_phase=NULL,ncores=NULL,batch_size=4096L) {
  begun<-proc.time()[['elapsed']];phase<-match.arg(phase)
  if(!inherits(library,"GCMSlibrary")) stop("Supply a GCMSlibrary.")
  if(identical(method,"all")) method<-search_methods
  if(!is.character(method)||!length(method)||any(!method%in%search_methods)) stop("Unknown spectral scoring method.")
  method<-unique(method)
  integer_option<-function(z) is.numeric(z)&&length(z)==1L&&is.finite(z)&&z>=1&&z<=.Machine$integer.max&&z==floor(z)
  if(!integer_option(n_matches)||!integer_option(batch_size)) stop("Use positive integer n_matches and batch_size.")
  if(!is.numeric(ri_tolerance)||length(ri_tolerance)!=1L||is.na(ri_tolerance)||ri_tolerance<0) stop("ri_tolerance must be nonnegative.")
  if(!is.null(unknown_phase)&&!identical(unknown_phase,phase)) stop("unknown_phase must be NULL or the chosen phase.")
  if(is.null(ncores)) ncores<-parallel::detectCores(logical=FALSE)
  if(is.na(ncores)) ncores<-1L
  if(!integer_option(ncores)) stop("ncores must be a positive integer.")
  query<-search_queries(x,calibration,ri,pseudocount,extrapolate,is.finite(ri_tolerance))
  if(nrow(query$queries)&&is.finite(ri_tolerance)&&all(is.na(query$queries$query_RI)))
    warning("No query RI is available; constrained matches remain unassigned.",call.=FALSE)
  con<-if(library$backend=="sqlite") library_connect(library) else NULL
  if(!is.null(con)) on.exit(DBI::dbDisconnect(con),add=TRUE)
  valid<-query$queries$query_RI[is.finite(query$queries$query_RI)]
  range<-if(is.finite(ri_tolerance)&&length(valid)) range(valid)+c(-1,1)*ri_tolerance else c(-Inf,Inf)
  indices<-search_ri(library,con,phase,unknown_phase,range)
  ids<-if(is.finite(ri_tolerance)) sort(unique(indices$id)) else seq_len(library$n)
  if(is.finite(ri_tolerance)&&!length(valid)) ids<-integer()
  hits<-list();by_id<-split(indices$value,indices$id)
  for(at in seq.int(1L,max(1L,length(ids)),by=as.integer(batch_size))) {
    batch<-ids[seq.int(at,min(length(ids),at+batch_size-1L))]
    if(!length(ids)) break
    prepared<-search_matrix(library,con,batch,query$mass)
    ris<-lapply(as.character(batch),function(id) {a<-by_id[[id]];if(is.null(a)) numeric() else as.numeric(a)})
    out<-library_search_cpp(t(query$S),prepared$matrix,query$mass,ris,query$queries$query_RI,
      ri_tolerance,as.integer(n_matches),match(method,search_methods),as.integer(ncores))
    if(nrow(out)) {
      out$library_mass_coverage<-out$library_total/prepared$raw_total[out$library_id]
      out$library_id<-batch[out$library_id];hits[[length(hits)+1L]]<-out
    }
  }
  hits<-if(length(hits)) do.call(rbind,hits) else data.frame(query=integer(),library_id=integer(),
    method=integer(),distance=numeric(),library_RI=numeric(),library_total=numeric(),library_mass_coverage=numeric())
  if(nrow(hits)) {
    hits$query_RI<-query$queries$query_RI[hits$query];hits$delta_RI<-hits$library_RI-hits$query_RI
    hits<-hits[order(hits$query,hits$method,hits$distance,abs(hits$delta_RI),hits$library_id,na.last=TRUE),,drop=FALSE]
    hits$rank<-ave(hits$query,paste(hits$query,hits$method),FUN=seq_along)
    hits<-hits[hits$rank<=n_matches,,drop=FALSE]
    kept<-unique(hits$library_id)
    meta<-if(is.null(con)) library$metadata[match(kept,library$metadata$id),,drop=FALSE] else
      DBI::dbGetQuery(con,paste0("SELECT id,name,cas,inchikey,source,source_id,origin FROM spectra WHERE id IN (",paste(kept,collapse=","),")"))
    a<-match(hits$library_id,meta$id)
    for(nm in intersect(c("name","cas","inchikey","source","source_id","origin"),names(meta))) hits[[nm]]<-meta[[nm]][a]
    nearest<-match(paste(hits$library_id,hits$library_RI),paste(indices$id,indices$value))
    hits$ri_phase<-indices$phase_class[nearest]
    hits$ri_assumed<-!is.na(hits$ri_phase)&hits$ri_phase=="unknown"&!is.null(unknown_phase)
    hits$ri_origin<-if('origin'%in%names(indices)) indices$origin[nearest] else NA_character_
    hits$ri_source<-if('source'%in%names(indices)) indices$source[nearest] else NA_character_
    hits$query_total<-query$queries$query_total[hits$query]
    hits$compound<-query$queries$compound[hits$query]
    hits$score<-1000*ifelse(hits$method%in%c(1L,2L),1-hits$distance,
      ifelse(hits$method==3L,1-hits$distance^2,ifelse(hits$method==4L,1-hits$distance^2,1/(1+hits$distance))))
    hits$method<-search_methods[hits$method];rownames(hits)<-NULL
  } else hits$method<-character()
  query$queries$status<-ifelse(query$queries$query_total<=0,"empty_spectrum",
    ifelse(is.finite(ri_tolerance)&is.na(query$queries$query_RI),"uncalibrated",
      ifelse(query$queries$query%in%hits$query,"matched","no_eligible_reference")))
  structure(list(matches=hits,queries=query$queries,calibration=query$calibration,
    settings=list(phase=phase,ri_tolerance=ri_tolerance,method=method,unknown_phase=unknown_phase,
      ncores=ncores,n_matches=n_matches,mass=query$mass),
    seconds=proc.time()[['elapsed']]-begun),class="GCMSmatches")
}

#' @rdname library_search
#' @param row.names,optional Additional data-frame arguments.
#' @param ... Additional method arguments.
#' @method as.data.frame GCMSmatches
#' @export
as.data.frame.GCMSmatches <- function(x,row.names=NULL,optional=FALSE,...) x$matches

#' @rdname library_search
#' @method print GCMSmatches
#' @export
print.GCMSmatches <- function(x,...) {
  cat(sum(x$queries$status=="matched"),"/",nrow(x$queries)," spectra annotated in ",
    format(x$seconds,digits=3)," s; ",x$settings$phase," RI tolerance ",x$settings$ri_tolerance,"\n",sep="")
  if(nrow(x$matches)) print(x$matches[x$matches$rank==1,
    c("compound","method","name","score","query_RI","library_RI","delta_RI")],row.names=FALSE)
  invisible(x)
}
