test_that("directory discovery uses acquisition times and separates standards", {
  skip_if_not_installed("ncdf4")
  d<-tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE))
  for(g in c("standards","beer1","beer2"))dir.create(file.path(d,g))
  t<-as.POSIXct("2026-11-01 08:00:00",tz="UTC");y<-matrix(1,12,3)
  write_GCMS_cdf(y,file.path(d,"beer1","z.cdf"),time=(1:12)/60,mass=40:42,acquired_at=t)
  write_GCMS_cdf(y,file.path(d,"beer2","a.cdf"),time=(1:12)/60,mass=40:42,acquired_at=t+3600)
  write_GCMS_cdf(y,file.path(d,"standards","alkanes.cdf"),time=(1:12)/60,mass=40:42,acquired_at=t-3600)
  z<-GCMS_batch_files(d,paste0(d,"_results"))
  expect_identical(z$manifest$sample_name,c("z","a"))
  expect_identical(z$manifest$group,c("beer1","beer2"))
  expect_true(all(basename(z$manifest$alkanes_file)=="alkanes.cdf"))
  expect_equal(as.numeric(z$metadata[[1]]$acquired_at),as.numeric(t))
  nc<-ncdf4::nc_open(file.path(d,"beer1","z.cdf"),write=TRUE)
  ncdf4::ncatt_put(nc,0,"experiment_date_time_stamp","");ncdf4::nc_close(nc)
  expect_error(GCMS_batch_files(d,paste0(d,"_results")),"Acquisition dates missing")
  z<-GCMS_batch_files(d,paste0(d,"_results"),sample_order=c("beer1/z.cdf","beer2/a.cdf"))
  expect_identical(z$manifest$sample_name,c("z","a"))
})

test_that("coupled directory workflows export original-channel sample and consensus results", {
  skip_if_not_installed("ncdf4");skip_if_not_installed("openxlsx")
  set.seed(54);n<-90L;m<-5L;kern<-convkernel("gaussian",13L,list(sigma=1.5));X<-convmatrix(kern,n)
  S<-matrix(rexp(2*m,1/3000),2,m)
  d<-tempfile();out<-paste0(d,"_results");dir.create(d)
  on.exit(unlink(c(d,out),recursive=TRUE))
  for(s in 1:2) {
    mu<-as.matrix(X[,c(25,65)]%*%(S*s))+2;y<-matrix(rpois(n*m,mu),n,m)
    write_GCMS_cdf(y,file.path(d,paste0(s,".cdf")),time=seq_len(n)/60,mass=40L+seq_len(m),
      acquired_at=as.POSIXct("2026-11-01",tz="UTC")+s*3600)
  }
  f<-GCMSdeco(d,output_dir=out,ncores=1L,nwindows=1L,dispersion=1,correlation="none",
    psf_control=list(kernel_est=kern,maxiter=1L),
    control=list(profile_cycles=0L,original_cleanup=FALSE),
    batch=list(coupling=list(cycles=1L,initial_cycles=1L,rounds=1L)))
  expect_s3_class(f,"GCMSbatch");expect_equal(nrow(f$samples),2L)
  expect_output(print(f),"GC/MS experiment")
  expect_equal(dim(fitted(f,sample=1L)),c(n,m));expect_true(is.finite(GIC(f)))
  expect_true(all(file.exists(file.path(out,c("consensus.msp","experiment.xlsx","samples.csv","areas.csv")))))
  expect_identical(openxlsx::getSheetNames(file.path(out,"experiment.xlsx")),c("samples","compounds","areas"))
  z<-readRDS(f$paths[1L]);expect_null(z$model$y);expect_null(z$model$mu)
  expect_equal(unname(fitted(z)),unname(fitted(f$coupled,sample=1L)),tolerance=1e-12)
  expect_equal(ncol(coef(f)),m)
  again<-GCMSdeco(d,output_dir=out,ncores=1L,nwindows=1L,dispersion=1,correlation="none",
    psf_control=list(kernel_est=kern,maxiter=1L),
    control=list(profile_cycles=0L,original_cleanup=FALSE),
    batch=list(resume=TRUE,coupling=list(cycles=1L,initial_cycles=1L,rounds=1L)))
  expect_true(all(again$samples$cache_hit))
  expect_equal(f$areas,again$areas)
})
