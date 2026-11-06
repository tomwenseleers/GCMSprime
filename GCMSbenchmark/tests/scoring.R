library(GCMSbenchmark)
internal <- asNamespace("GCMSbenchmark")
time <- c(.9,1,1.05,1.2)
truth <- list(sample=list(id=c("a","b"),time=time,mass=c(44,55),
  H=cbind(c(0,3,0,0),c(0,0,4,0)),P=diag(2),area=c(3,4),
  apex=c(1,1.05),acid=c(FALSE,TRUE)))
peaks <- internal$peak_table("sample",1:3,c(1,1.05,1.06),c(3,4,4))
f <- internal$records(peaks,list(rbind(c(44,1),c(44.01,1)),
  rbind(c(55,3)),rbind(c(55,2))),
  list(cbind(time,c(0,3,0,0)),cbind(time,c(0,0,4,0)),NULL),"integrated counts")
f$status <- "ok"
f$message <- ""
x <- structure(list(methods=list(test=f),files="sample.cdf"),class="GCMSbenchmark")
s <- scoreGCMSbenchmark(x,truth)
stopifnot(s$summary$TP==2L,s$summary$FP==1L,s$summary$FN==0L,
  s$summary$acid_TP==1L,s$summary$area_median_error_factor==0,
  all(s$matches$spectral_cosine==1),all(s$matches$profile_cosine==1))
stopifnot(is.na(scoreGCMSbenchmark(x)$summary$F1))

area <- internal$area_regression(data.frame(area_true=c(1,10,100,1,NA),
  area_estimated=c(5,50,500,0,10)))
stopifnot(nrow(area$data)==3L,isTRUE(all.equal(area$adjusted_R2,1)),
  isTRUE(all.equal(unname(stats::coef(area$fit)[2L]),1)))
stopifnot(is.na(internal$area_regression(data.frame(area_true=rep(1,3),
  area_estimated=1:3))$adjusted_R2))

f <- tempfile(fileext=".elu")
writeLines(c("NAME: |SC2|CN1|AM3|FR1-3|RT1|", "RE","0 3 0",
  "NUM PEAKS: 2", "(44,100 )(55,9 B0.1)"),f)
e <- internal$read_elu(f,"sample",time)
stopifnot(e$peaks$area==3,e$peaks$amdis_base_ion_signal==3,
  all(e$spectra[[1L]][,1L]==c(44,55)),all(e$profiles[[1L]][,1L]==time[1:3]))
