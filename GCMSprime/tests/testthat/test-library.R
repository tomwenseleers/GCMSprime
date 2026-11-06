msp_fixture <- function() {
  f<-tempfile(fileext=".msp")
  writeLines(c("Name: ethanol", "Synon: alcohol", "Synon: ethyl alcohol",
    "Formula: C2H6O", "CASNO: 64-17-5", "RI_DB5: 560", "Unknown: keep me",
    "Num Peaks: 4", '31 100 "fragment"; 45 60; 46.25 2e1; 31 5;',
    "NAME: compound B", "InChIKey: EXAMPLE", "RI_DBWAX: 1200", "RI_DB1: 700",
    "Num Peaks: 2", "40 1; 50 3;", ""),f)
  f
}
test_that("MSP parsing preserves peaks, duplicate fields and RI provenance", {
  f<-msp_fixture();on.exit(unlink(f))
  x<-readMSP(f,batch_size=1L)
  expect_s3_class(x,"GCMSlibrary");expect_s4_class(x$matrix,"dgCMatrix");expect_equal(x$n,2)
  expect_equal(libraryMetadata(x)$name,c("ethanol","compound B"))
  fields<-libraryMetadata(x,1,raw=TRUE)[[1]]
  expect_equal(fields$value[fields$tag=="Synon"],c("alcohol","ethyl alcohol"))
  expect_true("keep me"%in%fields$value)
  raw<-spectra(x,format="list")
  expect_equal(raw[[1]][,1],c(31,45,46.25,31));expect_equal(raw[[1]][,2],c(100,60,20,5))
  z<-spectra(x);expect_s4_class(z,"dgCMatrix");expect_equal(unname(z["31",1]),105)
  expect_equal(as.numeric(Matrix::colSums(spectra(x,normalize="l1"))),c(1,1))
  expect_equal(as.numeric(Matrix::colSums(spectra(x,normalize="l2")^2)),c(1,1))
  expect_true("46"%in%rownames(spectra(x,nominal=TRUE)))
  ri<-retentionIndices(x);expect_equal(ri$value,c(560,1200,700))
  expect_equal(retentionIndices(x,phase="DB5")$id,1L)
  expect_equal(spectra(x,c(2,1,2),format="list"),raw[c(2,1,2)])
})
test_that("SQLite libraries reproduce memory results and reopen lazily", {
  skip_if_not_installed("RSQLite")
  f<-msp_fixture();db<-tempfile(fileext=".sqlite");on.exit(unlink(c(f,db)))
  a<-readMSP(f,origin="predicted");b<-readMSP(f,backend="sqlite",db_file=db,batch_size=1L,origin="predicted")
  expect_equal(libraryMetadata(a),libraryMetadata(b),ignore_attr=TRUE)
  expect_equal(libraryMetadata(a,raw=TRUE),libraryMetadata(b,raw=TRUE),ignore_attr=TRUE)
  expect_equal(spectra(a),spectra(b));expect_equal(retentionIndices(a),retentionIndices(b),ignore_attr=TRUE)
  expect_equal(retentionIndices(b)$origin,rep("predicted",3))
  expect_equal(spectra(openLibrary(db)),spectra(a))
  expect_equal(spectra(unserialize(serialize(b,NULL)),2),spectra(a,2))
  expect_error(readMSP(f,"sqlite",db_file=db),"new db_file")
  expect_lt(as.numeric(object.size(b)),10000)
})
test_that("invalid MSP records fail without publishing a partial database", {
  f<-tempfile(fileext=".msp");db<-tempfile(fileext=".sqlite");on.exit(unlink(c(f,db)))
  writeLines(c("Name: incomplete","Num Peaks: 2","31 5;"),f)
  expect_error(readMSP(f),"Num Peaks mismatch")
  expect_error(readMSP(f,"sqlite",db_file=db),"Num Peaks mismatch");expect_false(file.exists(db))
  writeLines(c("Name: invalid","Num Peaks: 1","31 -1;"),f)
  expect_error(readMSP(f),"Invalid peak")
})
test_that("BOM, scientific notation and zero-peak records are supported", {
  f<-tempfile(fileext=".msp");on.exit(unlink(f))
  writeLines(c(paste0("\ufeff","Name: test"),"Num Peaks: 1","1.25e2 2.5e1;","", "Name: empty","Num Peaks: 0",""),f,useBytes=TRUE)
  x<-readMSP(f);expect_equal(x$n,2L)
  expect_equal(spectra(x,format="list")[[1]],cbind(mass=125,intensity=25))
  expect_equal(Matrix::colSums(spectra(x)),c(`1`=25,`2`=0))
})
test_that("GCMS-ID CSV predictions retain phase classes and metadata", {
  f<-tempfile(fileext=".csv");on.exit(unlink(f))
  z<-data.frame(compound_name=c("example_SemiStdNP","example_StdNP","example_StdPolar"),
    inchikey="EXAMPLE",compound_formula="C2H6O",retention_index_value=c(560,550,1000),
    exact_mass=46.04,peak_list="31 100; 45 60; 46 20;")
  utils::write.csv(z,f,row.names=FALSE)
  x<-readGCMSID(f,batch_size=1L)
  expect_equal(x$n,3L);expect_equal(retentionIndices(x)$phase_class,c("DB5","DB1","DBWAX"))
  expect_equal(libraryMetadata(x)$origin,rep("predicted",3))
  expect_true("GCMS-ID_peak_list"%in%libraryMetadata(x,raw=TRUE)[[1]]$tag)
  expect_equal(as.numeric(Matrix::colSums(spectra(x))),rep(180,3))
})

test_that("MoNA quoted comments populate searchable chemical metadata", {
  f<-tempfile(fileext=".msp");on.exit(unlink(f))
  text<-'"SMILES=CCO" "cas=64-17-5" "InChI=InChI=1S/C2H6O" "license=CC BY-SA"'
  writeLines(c("Name: ethanol",paste0("Comments: ",text),"Num Peaks: 1","31 100;"),f)
  x<-readMSP(f);z<-libraryMetadata(x)
  expect_equal(z$smiles,"CCO");expect_equal(z$cas,"64-17-5")
  expect_equal(z$inchi,"InChI=1S/C2H6O");expect_equal(z$license,"CC BY-SA")
  expect_equal(libraryMetadata(x,raw=TRUE)[[1]]$value[2],text)
})

test_that("OpenChrom NIST comment indices retain their known column classes", {
  f<-tempfile(fileext=".msp");on.exit(unlink(f))
  comment<-"SemiStdNP=1000/12/3 StdNP=990 StdPolar=1300/20/5 ; NIST MS# 1"
  writeLines(c("Name: test",paste0("COMMENTS: ",comment),"RI: 0","Num Peaks: 1","43 100;"),f)
  a<-readMSP(f);ri<-retentionIndices(a)
  expect_equal(ri$value,c(1000,990,1300))
  expect_equal(ri$phase_class,c("DB5","DB1","DBWAX"))
  expect_equal(ri$raw,c("1000/12/3","990","1300/20/5"))
  expect_equal(libraryMetadata(a,raw=TRUE)[[1]]$value[2],comment)
  skip_if_not_installed("RSQLite")
  db<-tempfile(fileext=".sqlite");on.exit(unlink(db),add=TRUE)
  b<-readMSP(f,backend="sqlite",db_file=db)
  expect_equal(retentionIndices(b),ri,ignore_attr=TRUE)
})

test_that("quoted NIST comment RIs are imported once", {
  f<-tempfile(fileext=".msp");on.exit(unlink(f))
  writeLines(c("Name: test",'Comments: "SemiStdNP=1000" "StdNP=990"',"Num Peaks: 1","43 100;"),f)
  expect_equal(retentionIndices(readMSP(f))$value,c(1000,990))
})
