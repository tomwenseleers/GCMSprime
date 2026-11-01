# Ion rarity weights, count balancing and angular metrics on separate cohorts.
library(GCMSprime)
library_file <- "codex_sandbox/library_import/mona_import.rds"
out <- "codex_sandbox/library_search_weights"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
lib <- readRDS(library_file)
stopifnot(lib$backend=="memory")
mass<-33:550;a<-GCMSprime:::search_matrix(lib,NULL,seq_len(lib$n),mass)$matrix
key<-ifelse(!is.na(lib$metadata$inchikey)&nzchar(lib$metadata$inchikey),lib$metadata$inchikey,
  ifelse(!is.na(lib$metadata$cas)&nzchar(lib$metadata$cas),paste0('CAS:',lib$metadata$cas),NA))
groups<-split(seq_len(lib$n),key);groups<-groups[lengths(groups)>1]
set.seed(20261101);groups<-groups[sample(length(groups))]
query_id<-integer();excluded<-integer()
for(g in groups) {
  total<-Matrix::colSums(a[,g,drop=FALSE]);g<-g[total>0]
  if(length(g)<2) next
  x<-as.matrix(a[,g,drop=FALSE]);pp<-sweep(x,2,colSums(x),'/')
  if(sum(pp[,1]>0)<8||max(colSums(abs(pp-pp[,1])))<.05) next
  excluded<-c(excluded,g[colSums(abs(pp-pp[,1]))<=.05])
  query_id<-c(query_id,g[1]);if(length(query_id)==500) break
}
q<-t(as.matrix(a[,query_id,drop=FALSE]));colnames(q)<-mass
excluded<-unique(excluded);a[,excluded]<-0;a<-Matrix::drop0(a)
label<-ifelse(seq_along(query_id)<=250,'development','validation')
lib$matrix<-a;lib$mass<-mass
d<-list(q=q,key=key,truth=key[query_id],label=label,mass=mass)
cases<-list(none=NULL,idf=library_ion_weights(lib,mass=d$mass,power=1,ncores=8),
  idf_squared=library_ion_weights(lib,mass=d$mass,power=2,ncores=8))
cat('Weight preparation seconds:',attr(cases$idf,'seconds'),'; spectra:',attr(cases$idf,'spectra'),'\n')
summary<-list();matches<-list()
for(name in names(cases)) {
  z<-library_search(d$q,lib,ri_tolerance=Inf,method='all',n_matches=5,ncores=8,
    ion_weights=cases[[name]],reference_count='query')
  m<-z$matches;m$correct<-!is.na(d$key[m$library_id])&d$key[m$library_id]==d$truth[m$query]
  for(pair in list(c('cosine','angular'),c('stein_scott','stein_scott_angular'),c('jensen_shannon','multinomial'))) {
    a<-m[m$method==pair[1],c('query','library_id','rank')];b<-m[m$method==pair[2],c('query','library_id','rank')]
    rownames(a)<-rownames(b)<-NULL
    stopifnot(isTRUE(all.equal(a,b)))
  }
  matches[[name]]<-m
  for(cohort in unique(d$label)) for(method in z$settings$method) {
    h<-m[m$method==method&d$label[m$query]==cohort,]
    summary[[length(summary)+1L]]<-data.frame(weight=name,cohort=cohort,method=method,
      queries=sum(d$label==cohort),top1=mean(h$correct[h$rank==1]),
      top5=mean(tapply(h$correct,h$query,any)),seconds=z$seconds)
  }
}
raw<-library_search(d$q,lib,ri_tolerance=Inf,method='multinomial',reference_count='library',n_matches=5,ncores=8)
h<-raw$matches;h$correct<-!is.na(d$key[h$library_id])&d$key[h$library_id]==d$truth[h$query]
for(cohort in unique(d$label)) {
  hh<-h[d$label[h$query]==cohort,]
  summary[[length(summary)+1L]]<-data.frame(weight='none',cohort=cohort,method='multinomial_library',
    queries=sum(d$label==cohort),top1=mean(hh$correct[hh$rank==1]),
    top5=mean(tapply(hh$correct,hh$query,any)),seconds=raw$seconds)
}
summary<-do.call(rbind,summary)
write.csv(summary,file.path(out,'production_accuracy.csv'),row.names=FALSE)
saveRDS(matches,file.path(out,'production_matches.rds'));saveRDS(cases,file.path(out,'production_weights.rds'))
times<-do.call(rbind,lapply(c('stein_scott','cosine','angular','stein_scott_angular','multinomial','jensen_shannon','hellinger'),function(m) {
  z<-library_search(d$q,lib,ri_tolerance=Inf,method=m,n_matches=5,ncores=8,ion_weights=cases$idf)
  data.frame(method=m,queries=nrow(d$q),seconds=z$seconds)
}))
write.csv(times,file.path(out,'production_timing.csv'),row.names=FALSE)
print(summary[summary$cohort=='validation',],row.names=FALSE);print(times)
