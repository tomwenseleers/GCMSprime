library(GCMSbenchmark)
internal <- asNamespace("GCMSbenchmark")
time <- c(0,1,2)
truth <- list(sample=list(id="a",time=time,mass=c(44,55),H=matrix(c(0,1,0),3,1),
  P=matrix(c(1,0),1,2),area=10,apex=1,reference_apex=1,acid=FALSE,
  ri_calibration=data.frame(time=c(0,2),ri=c(900,1100)),
  reference_ri_calibration=data.frame(time=c(0,2),ri=c(900,1100))))
reference <- list(id=c("a","b"),P=diag(2),mass=c(44,55),ri=c(1000,1005))
peaks <- internal$peak_table("sample",1:3,c(1,1.02,1.05),c(10,2,4))
fit <- internal$records(peaks,list(matrix(c(44,1),1,2),
  matrix(c(44,1),1,2),matrix(c(55,1),1,2)),area_kind="integrated counts")
fit$status <- "ok"
fit$message <- ""
x <- structure(list(methods=list(test=fit),files="sample.cdf"),class="GCMSbenchmark")
y <- scoreGCMSbenchmark(x,truth,matching="reciprocal",reference=reference,ri_tolerance=10)
stopifnot(y$summary$TP==1L,y$summary$FP==2L,y$summary$FN==0L,
  y$summary$duplicate_calls==1L,y$summary$false_library_hits==1L,
  identical(y$matches$estimate_id,"1"),y$matches$area_estimated==10,
  sum(y$identifications$reciprocal)==2L)
# The absent library entry remains eligible, and its hit is a false call.
stopifnot(y$identifications$reference_id[3L]=="b",!y$identifications$present[3L])
# A strict RI window rejects a perfect spectral match at the wrong RI.
shifted <- x
shifted$methods$test$peaks$rt <- shifted$methods$test$peaks$rt+.4
z <- scoreGCMSbenchmark(shifted,truth,matching="reciprocal",reference=reference,ri_tolerance=10)
stopifnot(z$summary$TP==0L,z$summary$FP==3L,z$summary$FN==1L)
# Existing ELU files remain readable across platforms, while execution is Windows-only.
amdis <- internal$bench_amdis
environment(amdis) <- list2env(list(.Platform=list(OS.type="unix")),parent=environment(amdis))
err <- tryCatch(amdis("sample.cdf",tempdir(),1L,list(),1,list()),error=identity)
stopifnot(inherits(err,"error"),grepl("requires Windows",conditionMessage(err),fixed=TRUE))
