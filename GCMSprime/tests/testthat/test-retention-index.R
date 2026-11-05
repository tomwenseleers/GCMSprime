test_that("retention coordinates use onset and narrow duration", {
  t<-seq(0,8,length.out=1601);mass<-33:200
  Y<-matrix(0,length(t),length(mass),dimnames=list(as.character(t),as.character(mass)))
  centers<-1:5;carbon<-7:11
  for(k in seq_along(centers)) {
    h<-1e4*exp(-.5*((t-centers[k])/.035)^2)
    for(j in c(43,57,71,85,99)) Y[,as.character(j)]<-Y[,as.character(j)]+h
    Y[,as.character(14*carbon[k]+2)]<-Y[,as.character(14*carbon[k]+2)]+.1*h
  }
  z<-calibrate_ri(Y,carbon=carbon)
  expect_equal(z$anchors$carbon,carbon)
  expect_true(all(z$anchors$accepted))
  expect_equal(z$anchors$retention_time,centers,tolerance=.005)
  expect_true(all(is.finite(z$anchors$start)))
  expect_true(all(z$anchors$narrow_width>0))
  expect_true(all(diff(predict(z,seq(1,5,length.out=100)))>=0,na.rm=TRUE))
  expect_match(capture.output(print(z)),"alkane anchors")
})

test_that("manual anchors, missing times and extrapolation behave consistently", {
  a<-data.frame(carbon=7:11,retention_time=1:5,narrow_width=.1)
  z<-calibrate_ri(anchors=a)
  expect_equal(predict(z,1:5),100*(7:11),tolerance=1e-5)
  expect_true(all(is.na(predict(z,c(0,6,NA)))))
  expect_equal(predict(z,c(0,6,NA),extrapolate=TRUE),c(600,1200,NA),tolerance=1e-5)
  expect_error(calibrate_ri(anchors=a[1:3,]),"four")
  a$carbon[4]<-a$carbon[3]
  expect_error(calibrate_ri(anchors=a),"increasing")
  expect_error(calibrate_ri(anchors=data.frame(x=1)),"carbon")
  expect_error(calibrate_ri(control=list(fraction=0)),"Invalid")
  expect_error(calibrate_ri(control=list(typo=0)),"Unknown")
  expect_error(calibrate_ri(matrix(0,10,5)),"numeric mass")
})

test_that("RI calibration is exact natural cubic interpolation", {
  a<-data.frame(carbon=c(7,8,10,11,14),retention_time=c(1,2.2,4.3,5.7,9),narrow_width=.1)
  z<-calibrate_ri(anchors=a)
  f<-stats::splinefun(a$retention_time,100*a$carbon,method='natural')
  t<-seq(.5,10,length.out=101)
  expect_identical(z$interpolation,'natural')
  expect_equal(predict(z,t,extrapolate=TRUE),f(t),tolerance=1e-12)
  expect_equal(z$anchors$fitted_RI,100*a$carbon,tolerance=1e-12)
})

test_that("irregular anchors keep their exact indices when the natural cubic reverses", {
  a<-data.frame(carbon=7:10,retention_time=c(1,2,2.01,3))
  z<-calibrate_ri(anchors=a)
  expect_identical(z$interpolation,'hyman')
  expect_equal(predict(z,a$retention_time),100*a$carbon,tolerance=1e-12)
  expect_true(all(diff(predict(z,seq(1,3,length.out=500)))>=0))
})

test_that("calibrated PSF supplies the width at the selected relative height", {
  offset<-seq(-.2,.2,length.out=1001);g<-exp(-.5*(offset/.02)^2)
  psf<-list(psf_kernels=matrix(g,nrow=1),psf_offsets=offset)
  a<-data.frame(carbon=7:11,retention_time=1:5)
  z<-calibrate_ri(anchors=a,psf=psf)
  expected<-2*.02*sqrt(2*log(100))
  # Linear crossing interpolation is accurate to less than 1e-5 minutes.
  expect_lt(max(abs(GCMSprime:::ri_predict_width(z,c(0,2,7))-expected)),1e-5)
  expect_identical(z$width_source,"calibrated_psf")
  psf$psf_kernels<-rbind(g,g);psf$peak_timepoints<-c(1,5)
  z<-calibrate_ri(anchors=a,psf=psf)
  expect_lt(abs(GCMSprime:::ri_predict_width(z,3)-expected),1e-5)
  psf$psf_offsets<-1:3
  expect_error(calibrate_ri(anchors=a,psf=psf),"psf_kernels")
})

test_that("time-varying PSF widths retain their local scale", {
  t<-c(1,3,5,7);offset<-seq(-.5,.5,length.out=2001)
  sigma<-.015*exp(.1*t)
  psf<-list(psf_kernels=t(vapply(sigma,function(s) exp(-.5*(offset/s)^2),numeric(length(offset)))),
    psf_offsets=offset,psf_timepoints=t)
  z<-calibrate_ri(anchors=data.frame(carbon=7:11,retention_time=1:5),psf=psf)
  expect_lt(max(abs(GCMSprime:::ri_predict_width(z,t)-2*sigma*sqrt(2*log(100)))),1e-5)
  expect_gt(GCMSprime:::ri_predict_width(z,6),GCMSprime:::ri_predict_width(z,2))
})
