test_that("balanced multinomial deviance preserves query counts and agrees with JS", {
  a<-search_fixture()
  z<-librarySearch(a$q,a$lib,ri_tolerance=Inf,method=c("jensen_shannon","multinomial"),n_matches=4)
  j<-z$matches[z$matches$method=="jensen_shannon",]
  d<-z$matches[z$matches$method=="multinomial",]
  expect_equal(j$library_id,d$library_id)
  expect_equal(d$distance,4*d$query_total*log(2)*j$distance^2,tolerance=1e-7)
  expect_equal(d$reference_total,d$query_total)
  expect_true(all(d$reference_count=="query"))
  expect_equal(d$library_total[d$library_id==1],rep(170,2))
  b<-a$lib;b$metadata$pseudocount<-c(5000,NA,NA,NA)
  known<-librarySearch(a$q,b,ri_tolerance=Inf,method="multinomial",n_matches=4)
  expect_equal(known$matches$reference_total[known$matches$library_id==1],rep(5000,2))
  expect_true(all(known$matches$reference_count[known$matches$library_id==1]=="known"))
  force<-librarySearch(a$q,b,ri_tolerance=Inf,method="multinomial",n_matches=4,reference_count="query")
  expect_equal(force$matches$reference_total,force$matches$query_total)
})

test_that("weighted distances agree with direct Jensen gaps and angular formulas", {
  a<-search_fixture();x<-a$q[2,,drop=FALSE];w<-c('43'=.5,'57'=1.2,'71'=3,'91'=2)
  p<-as.numeric(x)/sum(x);q<-c(100,50,20,0)/170;mass<-as.numeric(names(w))
  wc<-function(x,y) sum(w*x*y)/sqrt(sum(w*x^2)*sum(w*y^2))
  kl<-function(p,q) sum(w[p>0]*p[p>0]*log(p[p>0]/q[p>0]))
  js<-(kl(p,(p+q)/2)+kl(q,(p+q)/2))/2
  for(mode in c("query","library")) {
    N<-sum(x);M<-if(mode=="query") N else 170;r<-(N*p+M*q)/(N+M)
    truth<-c(1-wc(p^.6*mass^3,q^.6*mass^3)^2,1-wc(p,q)^2,
      sqrt(sum(w*(sqrt(p)-sqrt(q))^2)/(2*max(w))),sqrt(js/(max(w)*log(2))),
      2*(N*kl(p,r)+M*kl(q,r)),acos(wc(p,q)),acos(wc(p^.6*mass^3,q^.6*mass^3)))
    z<-librarySearch(x,a$lib,ri_tolerance=Inf,method="all",n_matches=4,reference_count=mode,ion_weights=w)
    h<-z$matches[z$matches$library_id==1,]
    expect_equal(h$distance,truth,tolerance=1e-7)
    expect_true(all(z$matches$distance>=0))
    z2<-librarySearch(x,a$lib,ri_tolerance=Inf,method="all",n_matches=4,reference_count=mode,
      ion_weights=w,ncores=8,batch_size=1)
    expect_equal(z$matches,z2$matches)
  }
  z<-librarySearch(x,a$lib,ri_tolerance=Inf,method="all",n_matches=4)
  one<-librarySearch(x,a$lib,ri_tolerance=Inf,method="all",n_matches=4,ion_weights=rep(1,4))
  expect_equal(z$matches,one$matches)
})

test_that("angular distances are metric and retain cosine rankings", {
  a<-search_fixture();w<-c(.5,1.2,3,2)
  z<-librarySearch(a$q,a$lib,ri_tolerance=Inf,method="all",n_matches=4,ion_weights=w)
  for(pair in list(c("cosine","angular"),c("stein_scott","stein_scott_angular"))) {
    x<-z$matches[z$matches$method==pair[1],];y<-z$matches[z$matches$method==pair[2],]
    expect_equal(x$library_id,y$library_id)
    expect_equal(y$distance,acos(sqrt(pmax(0,1-x$distance))),tolerance=1e-7)
  }
  set.seed(34);Q<-matrix(rexp(3*4),3,4,dimnames=list(NULL,c(43,57,71,91)))
  lib<-a$lib;lib$matrix<-as(Matrix::Matrix(t(Q),sparse=TRUE),"dgCMatrix");lib$mass<-as.numeric(colnames(Q))
  lib$n<-3L;lib$metadata<-lib$metadata[1:3,];lib$metadata$id<-1:3
  z<-librarySearch(Q,lib,ri_tolerance=Inf,method="angular",n_matches=3,ion_weights=w)
  D<-matrix(0,3,3);D[cbind(z$matches$query,z$matches$library_id)]<-z$matches$distance
  expect_equal(D,t(D),tolerance=1e-7);expect_equal(diag(D),rep(0,3),tolerance=1e-7)
  for(i in 1:3) for(j in 1:3) for(k in 1:3) expect_lte(D[i,k],D[i,j]+D[j,k]+1e-7)
})

test_that("ion frequencies are reusable and independent of storage and threads", {
  a<-search_fixture();b<-search_fixture("sqlite")
  w<-libraryIonWeights(a$lib,mass=c(43,57,71,91),ncores=1)
  df<-c(3,2,1,1);truth<-1+log(5/(df+1));truth<-truth/mean(truth)
  expect_equal(as.numeric(w),truth)
  expect_equal(as.numeric(attr(w,"frequency")),df);expect_equal(attr(w,"spectra"),4L)
  ww<-libraryIonWeights(b$lib,mass=c(43,57,71,91),ncores=8,batch_size=1)
  expect_equal(as.numeric(w),as.numeric(ww))
  expect_equal(attr(w,"frequency"),attr(ww,"frequency"))
  expect_equal(as.numeric(libraryIonWeights(a$lib,mass=c(43,57,71,91),power=0)),rep(1,4))
  expect_equal(GCMSprime:::search_weights(c('57'=2),c(43,57,71)),c(1,2,1))
  expect_error(librarySearch(a$q,a$lib,ion_weights=rep(0,4)),"positive weight")
  expect_error(librarySearch(a$q,a$lib,ion_weights=c(1,2)),"one per")
  expect_error(librarySearch(a$q,a$lib,ion_weights=c('bad'=1)),"numeric mass")
  q<-matrix(c(1,0,0,0),1,dimnames=list(NULL,c(43,57,71,91)))
  z<-librarySearch(q,a$lib,ri_tolerance=Inf,method="all",ion_weights=c(0,1,1,1))
  expect_equal(nrow(z$matches),0L);expect_equal(z$queries$status,"empty_weighted_spectrum")
})

test_that("large sample totals do not amplify a cancellation floor", {
  a<-search_fixture();q<-matrix(c(100,50,20,0),1,dimnames=list(NULL,c(43,57,71,91)))
  z<-librarySearch(q,a$lib,ri_tolerance=Inf,method="multinomial",pseudocount=1e12,n_matches=4)
  expect_lt(z$matches$distance[z$matches$library_id==1],1e-12)
  expect_equal(z$matches$reference_total,rep(1e12,4))
  expect_equal(z$matches$query_total,rep(1e12,4))
})
