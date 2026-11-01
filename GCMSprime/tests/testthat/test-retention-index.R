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
