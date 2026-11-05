test_that("cache signatures retain standard provenance and compact cached geometry", {
  sample<-tempfile();standard<-tempfile()
  writeLines('sample',sample);writeLines('standard',standard)
  on.exit(unlink(c(sample,standard)))
  geometry<-list(psf_kernels=matrix(c(1,2,1)/4,1),psf_offsets=-1:1,peak_timepoints=2)
  args<-list(family=quasipoisson('identity'),ncores=16L,
    psf_control=list(method='standard',standard=list(calibration=geometry,spec=list(noise='fixed'))))
  signature<-GCMS_batch_signature(sample,standard,args)
  expect_equal(signature$standard$size,unname(file.info(standard)$size))
  expect_identical(signature$args$psf_control$standard,geometry)
  args$ncores<-1L
  expect_identical(signature,GCMS_batch_signature(sample,standard,args))
})

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

test_that("sample workers preserve acquisition order, shared maps and full spectra", {
  skip_if_not_installed("ncdf4");skip_if_not_installed("openxlsx")
  set.seed(661);n<-100L;m<-8L;kernel<-dnorm(-7:7,sd=1.6)
  X<-L0nmf::convmatrix(kernel,n);S<-matrix(runif(2*m,1500,4500),2,m)
  d<-tempfile();dir.create(d);outputs<-c(tempfile(),tempfile())
  on.exit(unlink(c(d,outputs),recursive=TRUE))
  for(s in 1:2)write_GCMS_cdf(matrix(rpois(n*m,as.matrix(X[,c(28,72)]%*%S)+2),n,m),
    file.path(d,paste0(3-s,".cdf")),time=seq_len(n)/60,mass=40+seq_len(m),
    acquired_at=as.POSIXct("2026-11-01",tz="UTC")+s*3600)
  args<-list(cdf_file=d,nwindows=1L,dispersion=1,correlation="none",precluster=3L,
    psf_control=list(kernel_est=kernel,maxiter=1L),
    control=list(profile_cycles=0L,original_cleanup=FALSE),
    batch=list(coupling=list(cycles=1L,initial_cycles=1L,rounds=1L)))
  one<-do.call(GCMSdeco,c(args,list(output_dir=outputs[1],ncores=1L)))
  two<-do.call(GCMSdeco,c(args,list(output_dir=outputs[2],ncores=2L)))
  expect_identical(two$samples$sample_name,c("2","1"))
  expect_equal(two$parallel,list(ncores=2L,workers=2L,threads=1L))
  expect_identical(two$channel_groups,one$channel_groups)
  expect_equal(coef(two),coef(one),tolerance=1e-9)
  expect_equal(two$areas,one$areas,tolerance=1e-9)
  expect_equal(GIC(two),GIC(one),tolerance=1e-9)
  expect_equal(dim(fitted(two,sample=1L)),c(n,m))
  expect_null(two$coupled$store$cache)
  expect_true(all(file.exists(two$paths)))
  args$batch$resume<-TRUE
  resumed<-do.call(GCMSdeco,c(args,list(output_dir=outputs[2],ncores=1L)))
  expect_true(all(resumed$samples$cache_hit))
  expect_equal(resumed$areas,two$areas,tolerance=1e-9)
  cal<-calibrate_ri(anchors=data.frame(carbon=8:12,retention_time=c(.2,.5,.8,1.1,1.4),narrow_width=.08))
  annotated<-annotate_ri(two,cal,msp_file=file.path(outputs[2],'consensus.msp'),
    xlsx_file=file.path(outputs[2],'experiment.xlsx'))
  expect_equal(annotated$areas,two$areas)
  expect_equal(coef(annotated),coef(two))
  expect_true(all(annotated$samples$RI_available>0))
  expect_true(all(is.finite(annotated$retention_indices$retention_index)))
  expect_true('alkane_calibration'%in%openxlsx::getSheetNames(file.path(outputs[2],'experiment.xlsx')))
  expect_match(paste(readLines(file.path(outputs[2],'consensus.msp')),collapse='\n'),
    'Retention_index: SemiStdNP=')
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
