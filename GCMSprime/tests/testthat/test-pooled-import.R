test_that("parallel pooled imports keep scan order and noise estimates", {
  skip_if_not_installed("ncdf4")
  root <- tempfile()
  dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  set.seed(1438)
  inputs <- file.path(root, paste0("sample",1:2,".cdf"))
  for (i in seq_along(inputs)) {
    y <- matrix(rpois(120L*6,30+i),120,6)
    writeGCMSCDF(y,inputs[i],time=seq_len(nrow(y))/60,mass=41:46)
  }
  read_group <- function(workers) {
    out <- file.path(root,paste0("workers",workers))
    dir.create(file.path(out,"samples"),recursive=TRUE)
    jobs <- lapply(seq_along(inputs),function(i) {
      list(id=paste0("sample",i),output=out,file=inputs[i],time_range=NULL,
           family=quasipoisson("identity"),dispersion=NULL,theta=NULL,
           correlation="none",threads=1L)
    })
    GCMS_batch_apply(jobs,workers,1L,fun=GCMS_read_pooled_sample)
  }
  one <- read_group(1L)
  two <- read_group(2L)
  for (i in seq_along(inputs)) {
    expect_identical(two[[i]]$time,one[[i]]$time)
    expect_identical(two[[i]]$total,one[[i]]$total)
    expect_equal(two[[i]]$spec$dispersion,one[[i]]$spec$dispersion)
    expect_identical(readRDS(two[[i]]$y),readRDS(one[[i]]$y))
  }
})
