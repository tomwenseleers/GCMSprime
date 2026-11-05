test_that("native spectral distances match direct calculations and preserve totals", {
  a<-search_fixture();z<-library_search(a$q,a$lib,ri=c(1002,1018),n_matches=2,method="all",ncores=1,reference_count="library")
  expect_equal(z$queries$query_total,c(34,27))
  expect_equal(nrow(z$matches),28L)
  b<-as.matrix(spectra(a$lib))[,1];b<-b[match(colnames(a$q),rownames(spectra(a$lib)))];b[is.na(b)]<-0
  x<-a$q[1,];p<-x/sum(x);q<-b/sum(b);m<-as.numeric(colnames(a$q))
  kl<-function(x,y) sum(x[x>0]*log(x[x>0]/y[x>0]))
  w<-function(x) x^.6*m^3
  truth<-c(1-sum(w(x)*w(b))^2/sum(w(x)^2)/sum(w(b)^2),
    1-sum(x*b)^2/sum(x^2)/sum(b^2),sqrt(1-sum(sqrt(p*q))),
    sqrt((kl(p,(p+q)/2)+kl(q,(p+q)/2))/(2*log(2))),
    2*(sum(x)*kl(p,(x+b)/sum(x+b))+sum(b)*kl(q,(x+b)/sum(x+b))),
    acos(min(1,sum(x*b)/sqrt(sum(x^2)*sum(b^2)))),
    acos(min(1,sum(w(x)*w(b))/sqrt(sum(w(x)^2)*sum(w(b)^2)))))
  one<-z$matches[z$matches$query==1&z$matches$library_id==1,]
  expect_equal(one$distance,truth,tolerance=1e-7)
  expect_equal(one$library_RI,rep(1000,7))
  expect_equal(one$query_total,rep(34,7))
  expect_equal(one$library_total,rep(170,7))
  expect_equal(one$library_mass_coverage,rep(1,7))
  z2<-library_search(a$q,a$lib,ri=c(1002,1018),method="all",n_matches=2,ncores=8,batch_size=1,reference_count="library")
  expect_equal(z$matches,z2$matches)
  z3<-library_search(a$q/rowSums(a$q),a$lib,ri=c(1002,1018),method="all",n_matches=2,pseudocount=c(34,27),reference_count="library")
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
  z<-library_search(x,a$lib,ri_tolerance=Inf,method="all",n_matches=4,reference_count="library")
  zz<-library_search(x*100,a$lib,ri_tolerance=Inf,method="all",n_matches=4,reference_count="library")
  k<-z$matches$method!="multinomial"
  expect_equal(z$matches$distance[k],zz$matches$distance[k],tolerance=1e-7)
  expect_false(isTRUE(all.equal(z$matches$distance[!k],zz$matches$distance[!k])))
  f<-tempfile(fileext=".msp")
  writeLines(c("Name: query","Num Peaks: 3","43 2;57 15;91 10;"),f)
  reverse<-library_search(matrix(c(100,50,20,0),1,dimnames=list(NULL,c(43,57,71,91))),
    read_msp(f),ri_tolerance=Inf,method="multinomial",reference_count="library")
  expect_equal(reverse$matches$distance,z$matches$distance[z$matches$method=="multinomial"&z$matches$library_id==1],tolerance=1e-7)
})

test_that("fit annotation preserves factors, orders compounds and exports all matches", {
  a<-search_fixture();time<-seq(0,5,length.out=501)
  H<-cbind(exp(-.5*((time-3)/.05)^2),exp(-.5*((time-1)/.05)^2))
  H<-sweep(H,2,colSums(H),'/')
  S<-rbind(c(100,50,20,0),c(0,70,0,50))
  model<-list(H=H,S=S,spec=list(name="poisson",dispersion=1,theta=1))
  x<-structure(list(model=model,time=time,mass=c(43,57,71,91),seconds=c(total=1),
    compounds=GCMSprime:::GCMS_compounds(model,time,"test"),cdf_file="test.cdf",sample_name="test"),class="GCMSdeco")
  expect_equal(x$compounds$compound,c(2L,1L))
  unconstrained<-library_search(x,a$lib,ri_tolerance=Inf)
  expect_true(nrow(unconstrained$matches)>0)
  out<-tempfile(fileext=".msp")
  excel<-if(requireNamespace("openxlsx",quietly=TRUE)) tempfile(fileext=".xlsx") else NULL
  z<-identify_compounds(x,a$lib,ri=c(1020,1000),n_matches=2,method="all",msp_file=out,xlsx_file=excel)
  expect_identical(z$model,x$model)
  expect_equal(z$compounds$match_library_id,c(2L,1L))
  expect_equal(z$compounds$match_score,c(1000,1000),tolerance=1e-6)
  expect_equal(z$compounds$retention_index,c(1020,1000))
  expect_match(paste(readLines(out,warn=FALSE),collapse="\n"),"Retention_index: SemiStdNP=1020")
  expect_match(paste(readLines(out,warn=FALSE),collapse="\n"),"Pseudocount: 120")
  if(!is.null(excel)) {
    expect_true(all(c("Compounds","Run","Matches")%in%openxlsx::getSheetNames(excel)))
    expect_equal(nrow(openxlsx::read.xlsx(excel,sheet="Matches")),28L)
  }
  cal<-calibrate_ri(anchors=data.frame(carbon=7:12,retention_time=0:5,narrow_width=.3))
  z<-library_search(x,a$lib,calibration=cal,ri_tolerance=Inf)
  expected<-vapply(c(2,1),function(j) GCMSprime:::ri_crossings(time,H[,j],which.max(H[,j]),.01)[['start']],0.)+.15
  expect_equal(z$queries$retention_coordinate,unname(expected))
  expect_equal(z$queries$query_RI,as.numeric(predict(cal,expected)))
})

test_that("combined deconvolution and identification works without changing the fit", {
  skip_if_not_installed("ncdf4")
  n<-80;kernel<-convkernel("gaussian",11L,list(sigma=1.2));X<-convmatrix(kernel,n)
  Y<-as.matrix(X[,c(20,60)]%*%rbind(c(700,300,50,0),c(0,500,100,400)))+3
  Y<-round(Y);colnames(Y)<-40:43;rownames(Y)<-seq_len(n)/60
  file<-tempfile(fileext=".cdf");write_GCMS_cdf(Y,file,time=seq_len(n)/60,mass=40:43)
  a<-search_fixture();out<-tempfile(fileext=".msp")
  cal<-calibrate_ri(anchors=data.frame(carbon=7:10,retention_time=c(0,.5,1,1.5),narrow_width=.1))
  common<-list(cdf_file=file,msp_file=out,family=poisson("identity"),nwindows=1L,baseline_df=4L,
    precluster=2L,ncores=1L,psf_control=list(kernel_est=kernel,maxiter=1L),
    control=list(broad_repeats=0L,full_sweeps=0L,final_polishes=0L))
  plain<-do.call(GCMSdeco,common)
  z<-do.call(GCMSdeco,c(common,list(library=a$lib,search=list(calibration=cal,ri_tolerance=300))))
  expect_s3_class(z$identification,"GCMSmatches")
  expect_identical(plain$model$H,z$model$H);expect_identical(plain$model$S,z$model$S)
  expect_equal(GIC(plain),GIC(z))
  expect_true("retention_index"%in%names(z$compounds))
  expect_error(GCMSdeco(file,out,search=list(phase="DB5")),"Supply library")
})

test_that("separate column metadata supplies the RI phase", {
  f<-tempfile(fileext=".msp")
  writeLines(c("Name: A","Column: DB-5ms","Retention_index: 1100","RI_DB1: 1090",
    "Num Peaks: 1","43 100;"),f)
  r<-retention_indices(read_msp(f))
  expect_equal(r$phase_class,c("DB5","DB1"))
})

test_that("explicit reference totals survive relative-intensity MSP export", {
  for(backend in c("memory","sqlite")) {
    if(backend=="sqlite"&&!requireNamespace("RSQLite",quietly=TRUE)) next
    f<-tempfile(fileext=".msp")
    writeLines(c("Name: count spectrum","Pseudocount: 5000","Num Peaks: 2","43 999;57 999;"),f)
    lib<-read_msp(f,backend=backend,db_file=if(backend=="sqlite") tempfile(fileext=".sqlite") else NULL)
    expect_equal(library_metadata(lib)$pseudocount,5000)
    expect_equal(sum(spectra(lib)),1998) # source arrays remain exact
    z<-library_search(matrix(c(250,250),1,dimnames=list(NULL,c(43,57))),lib,
      ri_tolerance=Inf,method="multinomial")
    expect_equal(z$matches$query_total,500)
    expect_equal(z$matches$library_total,5000)
    expect_equal(z$matches$distance,0,tolerance=1e-7)
  }
})

test_that("explicit benchmark phase assignments are visible in reports", {
  a<-search_fixture();a$lib$indices$phase_class<-"unknown"
  time<-seq(0,2,length.out=101);H<-matrix(exp(-.5*((time-1)/.1)^2),ncol=1)
  H<-H/sum(H);model<-list(H=H,S=matrix(c(100,50,20,0),1),
    spec=list(name="poisson",dispersion=1,theta=1))
  x<-structure(list(model=model,time=time,mass=c(43,57,71,91),seconds=c(total=1),
    compounds=GCMSprime:::GCMS_compounds(model,time,"phase test"),cdf_file="test.cdf"),class="GCMSdeco")
  msp<-tempfile(fileext=".msp")
  excel<-if(requireNamespace("openxlsx",quietly=TRUE)) tempfile(fileext=".xlsx") else NULL
  z<-identify_compounds(x,a$lib,ri=1000,unknown_phase="DB5",msp_file=msp,xlsx_file=excel)
  expect_true(z$compounds$match_ri_assumed)
  expect_equal(z$compounds$match_ri_phase,"unknown")
  expect_match(paste(readLines(msp,warn=FALSE),collapse="\n"),"reference RI column assumed DB5")
  if(!is.null(excel)) {
    run<-openxlsx::read.xlsx(excel,sheet="Run")
    expect_equal(run$value[run$field=="unlabelled reference RI phase"],"explicitly assumed DB5")
  }
})

test_that("indexed RI batches preserve duplicate values, gaps and unindexed spectra", {
  a<-search_fixture()
  a$lib$indices<-rbind(a$lib$indices,data.frame(a$lib$indices[1,,drop=FALSE]))
  a$lib$indices$value[nrow(a$lib$indices)]<-1004
  for(tolerance in c(Inf,30)) {
    full<-library_search(a$q,a$lib,ri=c(1003,1020),ri_tolerance=tolerance,
      method="all",n_matches=4,ncores=1)
    batched<-library_search(a$q,a$lib,ri=c(1003,1020),ri_tolerance=tolerance,
      method="all",n_matches=4,ncores=2,batch_size=1)
    expect_equal(full$matches,batched$matches)
    hit<-full$matches[full$matches$library_id==1&full$matches$query==1,]
    expect_equal(unique(hit$library_RI),1004)
  }
})
