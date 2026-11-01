search_fixture <- function(backend="memory") {
  f<-tempfile(fileext=".msp")
  writeLines(c("Name: A","CAS: 1","RI_DB5: 1000","RI_DB1: 900","RI_DB5: 1010","Num Peaks: 3",
    "43 100;57 50;71 20;","Name: B","RI_DB5: 1020","Num Peaks: 2","57 70;91 50;",
    "Name: C","RI: 1000","Num Peaks: 1","43 10;","Name: D","RI_DBWAX: 1000","Num Peaks: 1","43 10;"),f)
  lib<-read_msp(f,backend=backend,db_file=if(backend=="sqlite") tempfile(fileext=".sqlite") else NULL)
  list(lib=lib,q=matrix(c(20,10,4,0,2,15,0,10),2,byrow=TRUE,dimnames=list(NULL,c(43,57,71,91))))
}
