# Compare independent records of the same identity, with different spectral shapes.
# Import a local MSP once, then saveRDS(lib, library_file). Nothing is downloaded.
library(GCMSprime)
library_file <- "codex_sandbox/library_import/mona_import.rds"
out <- "codex_sandbox/librarySearch"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
lib <- readRDS(library_file)
# Form independent-record comparisons: different nominal shape, same full InChIKey/CAS.
mass<-33:550
a<-GCMSprime:::search_matrix(lib,NULL,seq_len(lib$n),mass)$matrix
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
  query_id<-c(query_id,g[1]);if(length(query_id)==250) break
}
q<-t(as.matrix(a[,query_id,drop=FALSE]));colnames(q)<-mass
truth<-key[query_id]
excluded<-unique(excluded)
lib$matrix[,excluded]<-0 # exclude query records and near-duplicate same-identity shapes
all<-librarySearch(q,lib,method='all',ri_tolerance=Inf,n_matches=5,ncores=8)
saveRDS(all,file.path(out,'replicate_matches.rds'))
tab<-do.call(rbind,lapply(all$settings$method,function(m) {
  z<-all$matches[all$matches$method==m,];correct<-key[z$library_id]==truth[z$query]
  data.frame(method=m,queries=nrow(q),top1=mean(correct[z$rank==1]),
    top5=mean(vapply(seq_len(nrow(q)),function(i) any(correct[z$query==i]),FALSE)))
}))
timing<-do.call(rbind,lapply(c(1,8),function(threads) do.call(rbind,lapply(all$settings$method,function(m) {
  z<-librarySearch(q,lib,method=m,ri_tolerance=Inf,n_matches=5,ncores=threads)
  data.frame(method=m,threads=threads,seconds=z$seconds)
}))))
write.csv(tab,file.path(out,'replicate_accuracy.csv'),row.names=FALSE)
write.csv(timing,file.path(out,'replicate_timing.csv'),row.names=FALSE)
cat('Independent records with different shapes:',nrow(q),'; reference catalogue:',lib$n-length(excluded),'\n')
print(tab);print(timing)
