search_fixture <- function(backend="memory") {
  f<-tempfile(fileext=".msp")
  writeLines(c("Name: A","CAS: 1","RI_DB5: 1000","RI_DB1: 900","RI_DB5: 1010","Num Peaks: 3",
    "43 100;57 50;71 20;","Name: B","RI_DB5: 1020","Num Peaks: 2","57 70;91 50;",
    "Name: C","RI: 1000","Num Peaks: 1","43 10;","Name: D","RI_DBWAX: 1000","Num Peaks: 1","43 10;"),f)
  lib<-read_msp(f,backend=backend,db_file=if(backend=="sqlite") tempfile(fileext=".sqlite") else NULL)
  list(lib=lib,q=matrix(c(20,10,4,0,2,15,0,10),2,byrow=TRUE,dimnames=list(NULL,c(43,57,71,91))))
}
test_that("native spectral distances match direct calculations and preserve totals", {
  a<-search_fixture();z<-library_search(a$q,a$lib,ri=c(1002,1018),n_matches=2,method="all",ncores=1)
  expect_equal(z$queries$query_total,c(34,27))
  expect_equal(nrow(z$matches),20L)
  b<-as.matrix(spectra(a$lib))[,1];b<-b[match(colnames(a$q),rownames(spectra(a$lib)))];b[is.na(b)]<-0
  x<-a$q[1,];p<-x/sum(x);q<-b/sum(b);m<-as.numeric(colnames(a$q))
  kl<-function(x,y) sum(x[x>0]*log(x[x>0]/y[x>0]))
  w<-function(x) x^.6*m^3
  truth<-c(1-sum(w(x)*w(b))^2/sum(w(x)^2)/sum(w(b)^2),
    1-sum(x*b)^2/sum(x^2)/sum(b^2),sqrt(1-sum(sqrt(p*q))),
    sqrt((kl(p,(p+q)/2)+kl(q,(p+q)/2))/(2*log(2))),
    2*(sum(x)*kl(p,(x+b)/sum(x+b))+sum(b)*kl(q,(x+b)/sum(x+b))))
  one<-z$matches[z$matches$query==1&z$matches$library_id==1,]
  expect_equal(one$distance,truth,tolerance=1e-7)
  expect_equal(one$library_RI,rep(1000,5))
  expect_equal(one$query_total,rep(34,5))
  expect_equal(one$library_total,rep(170,5))
  expect_equal(one$library_mass_coverage,rep(1,5))
  z2<-library_search(a$q,a$lib,ri=c(1002,1018),method="all",n_matches=2,ncores=8,batch_size=1)
  expect_equal(z$matches,z2$matches)
  z3<-library_search(a$q/rowSums(a$q),a$lib,ri=c(1002,1018),method="all",n_matches=2,pseudocount=c(34,27))
  expect_equal(z$matches,z3$matches,tolerance=1e-7)
  expect_equal(as.data.frame(z),z$matches)
})
test_that("RI phases, unknown indices, acquired masses and empty queries are handled", {
  a<-search_fixture()
  z<-library_search(a$q,a$lib,ri=c(1000,NA),ri_tolerance=0)
  expect_equal(z$matches$library_id,1L)
  expect_equal(z$queries$status,c("matched","uncalibrated"))
  z<-library_search(a$q,a$lib,ri=c(1000,1000),ri_tolerance=0,unknown_phase="DB5")
  expect_true(3L%in%z$matches$library_id)
  expect_false(4L%in%z$matches$library_id)
  z<-library_search(a$q,a$lib,ri=c(900,900),phase="DB1",ri_tolerance=0)
  expect_equal(unique(z$matches$library_id),1L)
  z<-library_search(a$q,a$lib,ri_tolerance=Inf,n_matches=4)
  expect_equal(sort(unique(z$matches$library_id)),1:4)
  z<-library_search(a$q[,c(1,2),drop=FALSE],a$lib,ri=c(1000,1000))
  expect_equal(z$matches$library_mass_coverage[z$matches$library_id==1],rep(150/170,2))
  expect_warning(z<-library_search(a$q,a$lib),"No query RI")
  expect_equal(nrow(z$matches),0L)
  expect_equal(z$queries$status,rep("uncalibrated",2))
  z<-library_search(a$q*0,a$lib,ri=c(1000,1000))
  expect_equal(z$queries$status,rep("empty_spectrum",2))
  expect_error(library_search(a$q,a$lib,n_matches=0),"positive integer")
  expect_error(library_search(a$q,a$lib,method="bogus"),"Unknown")
})
test_that("SQLite batching agrees with memory search", {
  skip_if_not_installed("RSQLite")
  a<-search_fixture();b<-search_fixture("sqlite")
  z<-library_search(a$q,a$lib,ri=c(1000,1020),method="all",n_matches=2)
  zz<-library_search(b$q,b$lib,ri=c(1000,1020),method="all",n_matches=2,batch_size=1)
  # Independently created temporary files have different source names.
  z$matches$source<-zz$matches$source<-NULL;z$matches$ri_source<-zz$matches$ri_source<-NULL
  expect_equal(z$matches,zz$matches)
})
test_that("count deviance is symmetric and shape metrics ignore intensity scales", {
  a<-search_fixture();x<-a$q[2,,drop=FALSE]
  z<-library_search(x,a$lib,ri_tolerance=Inf,method="all",n_matches=4)
  zz<-library_search(x*100,a$lib,ri_tolerance=Inf,method="all",n_matches=4)
  k<-z$matches$method!="multinomial"
  expect_equal(z$matches$distance[k],zz$matches$distance[k],tolerance=1e-7)
  expect_false(isTRUE(all.equal(z$matches$distance[!k],zz$matches$distance[!k])))
  f<-tempfile(fileext=".msp")
  writeLines(c("Name: query","Num Peaks: 3","43 2;57 15;91 10;"),f)
  reverse<-library_search(matrix(c(100,50,20,0),1,dimnames=list(NULL,c(43,57,71,91))),
    read_msp(f),ri_tolerance=Inf,method="multinomial")
  expect_equal(reverse$matches$distance,z$matches$distance[z$matches$method=="multinomial"&z$matches$library_id==1],tolerance=1e-7)
})
