library(GCMSbenchmark)
internal<-asNamespace("GCMSbenchmark")
if(requireNamespace("GCMSprime",quietly=TRUE)&&requireNamespace("ncdf4",quietly=TRUE)) {
  file<-file.path(tempdir(),"unit-test.cdf")
  GCMSprime::writeGCMSCDF(matrix(1,5,2),file,time=(0:4)*2/60,mass=41:42)
  id<-"unit-test"
  make<-function(area,height=10,spectrum=c(10,30),profile=c(0,5,10,5,0)) {
    p<-internal$peak_table(id,"a",2/60,area,height)
    internal$records(p,list(cbind(41:42,spectrum)),list(cbind((0:4)*2/60,profile)))
  }
  f<-make(40);f$peaks$quant_mass<-41
  mz<-internal$tic_units(f,"mzmine",file)
  stopifnot(mz$peaks$area==40*60/2*4,
    abs(sum(mz$spectra[[1]][,2])-1)<1e-15,
    mz$peaks$native_area==40)
  ms<-internal$tic_units(make(40),"MS-DIAL",file)
  stopifnot(ms$peaks$area==80,ms$peaks$height==40,
    max(ms$profiles[[1]][,2])==40)
  er<-internal$tic_units(make(40,spectrum=c(1000,3000)),"eRah",file)
  stopifnot(er$peaks$area==160,max(er$profiles[[1]][,2])==40)
  for(method in c("GCMSprime","AMDIS")) {
    x<-internal$tic_units(make(40),method,file)
    stopifnot(x$peaks$area==40,identical(internal$tic_units(x,method,file),x))
  }
  unlink(file)
}

set.seed(100)
x<-seq(1,3,length.out=100)
y<-2+.9*x+rnorm(100,sd=.03)
y[1:4]<-y[1:4]+5
seed<-.Random.seed
z<-internal$area_regression(data.frame(area_true=10^x,area_estimated=10^y))
stopifnot(inherits(z$fit,"lmrob"),abs(coef(z$fit)[2]-.9)<.03,
  identical(seed,.Random.seed),z$adjusted_R2>z$OLS_adjusted_R2)
band<-predict(z$fit,data.frame(log_true=x),interval="prediction")
stopifnot(all(is.finite(band)),all(band[,"lwr"]<band[,"upr"]))
