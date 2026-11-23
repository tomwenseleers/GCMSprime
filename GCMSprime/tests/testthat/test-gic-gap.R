test_that("CDF and directory fits return correctly scored diagnostics", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("openxlsx")
  d <- tempfile()
  dir.create(d)
  paths <- file.path(d,c("A.cdf","B.cdf"))
  outputs <- paste0(d,c("_joint","_independent"))
  on.exit(unlink(c(d,outputs),recursive=TRUE))
  kernel <- convKernel("gaussian",15L,list(sigma=2))
  n <- 80L
  H <- as.matrix(convMatrix(kernel,n)[,c(20L,57L)])
  colnames(H) <- c("one","two")
  P <- rbind(c(600,200,0),c(0,200,450))
  for(i in seq_along(paths))
    writeGCMSCDF(round(4+i*H%*%P),paths[i],time=seq_len(n)/60,mass=40:42,
                 acquired_at=as.POSIXct("2026-11-01",tz="UTC")+i*3600)
  args <- list(ncores=1L,nwindows=1L,baseline_df=1L,correlation="none",
               psf_control=list(kernel_est=kernel,maxiter=1L))
  f <- do.call(gcmsDeco,c(list(cdf_file=paths[1],msp_file=file.path(d,"A.msp"),
                              truth.profiles=H),args))
  expect_equal(f$gic_gap$fitted_gic,as.numeric(GIC(f)))
  expect_equal(f$model$gic_gap,f$gic_gap)
  truth <- list(B=H[,2:1],A=H)
  for (joint in c(TRUE,FALSE)) {
    batch <- do.call(gcmsDeco,c(list(cdf_file=d,output_dir=outputs[if(joint)1 else 2],
                         coupled=joint,truth.profiles=truth,batch=list(workers=1L)),args))
    if (joint) {
      expect_equal(batch$gic_gap$fitted_gic,as.numeric(GIC(batch)))
      expect_equal(batch$gic_gap,readRDS(file.path(batch$output_dir,"experiment.rds"))$gic_gap)
    } else {
      expect_named(batch$gic_gap,batch$samples$sample_id)
      expect_identical(gicGap(batch),batch$gic_gap)
      expect_equal(unname(vapply(batch$gic_gap,"[[",0,"fitted_gic")),batch$samples$GIC)
    }
  }
})
