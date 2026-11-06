bench_time <- function(expr) {
  t <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value=value,seconds=proc.time()[["elapsed"]]-t)
}

adapter_asset <- function(file) {
  path <- system.file("scripts",file,package="GCMSbenchmark")
  if (!nzchar(path)) stop("Missing adapter script ",file,call.=FALSE)
  path
}

bench_gcmsprime <- function(files,out,ncores,control,timeout,dots) {
  need_package("GCMSprime")
  if (!is.null(control$native)) {
    run <- readRDS(control$native)
    x <- run$result
    if (inherits(x,"error")) stop(conditionMessage(x),call.=FALSE)
    fit <- list(value=x,seconds=run$seconds)
  } else {
    if (length(files)==1L) {
      args <- list(cdf_file=files, msp_file=file.path(out,"peaks.msp"), ncores=ncores)
    } else {
      input <- file.path(out,"input")
      stage_inputs(files,input)
      args <- list(cdf_file=input,output_dir=file.path(out,"analysis"),ncores=ncores)
    }
    fit <- bench_time(do.call(GCMSprime::gcmsDeco,utils::modifyList(args,dots)))
  }
  x <- fit$value
  ans <- if (is.null(x$coupled)) gcmsprime_single(x,sample_id(files[1L])) else
    gcmsprime_batch(x,sample_id(files))
  ans$analysis <- fit$seconds
  ans
}

gcmsprime_single <- function(x,sample) {
  h <- as.matrix(x$model$H)
  s <- as.matrix(stats::coef(x$model))
  n <- nrow(x$compounds)
  stopifnot(nrow(s)==n, ncol(h)==n)
  j <- x$compounds$compound
  h <- h[,j,drop=FALSE]
  s <- s[j,,drop=FALSE]
  peaks <- peak_table(sample,seq_len(n),x$time[x$compounds$apex_scan],
    x$compounds$area,if(n) apply(h,2,max)*rowSums(s) else numeric())
  records(peaks,lapply(seq_len(n),function(j) positive_pairs(x$mass,s[j,])),
    lapply(seq_len(n),function(j) positive_pairs(x$time,h[,j]*sum(s[j,]))),
    area_kind="integrated counts",profile_kind="component TIC",
    native=x$seconds)
}

gcmsprime_batch <- function(x, sample) {
  model <- x$coupled
  ids <- sample_id(x$samples$cdf_file)
  if (!all(ids %in% sample) || anyDuplicated(ids))
    stop("GCMSprime sample names do not match the input files.",call.=FALSE)
  parts <- lapply(seq_along(ids), function(i) {
    f <- logconNMF::sampleNMF(model,i)
    h <- as.matrix(f$H)
    p <- as.matrix(stats::coef(f))
    mass <- as.numeric(colnames(p))
    area <- rowSums(p)
    j <- seq_along(f$compound_id)
    peaks <- peak_table(ids[i],j,f$time[f$apex_scan],area,
      vapply(j,function(k) max(h[,k])*area[k],0),group=f$compound_id)
    records(peaks,lapply(j,function(k) positive_pairs(mass,p[k,])),
      lapply(j,function(k) positive_pairs(f$time,h[,k]*area[k])))
  })
  combined <- combine_records(parts,"integrated counts","component TIC")
  combined$shared <- nrow(model$P)
  combined
}

bench_erah <- function(files,out,ncores,control,timeout,dots) {
  need_package("erah")
  workers <- min(ncores,length(files))
  if (workers>1L) {
    need_package("future")
    old <- future::plan()
    on.exit(future::plan(old))
    future::plan(future::multisession,workers=workers)
  }
  p <- do.call(erah::setDecPar,settings(list(min.peak.width=2,
    avoid.processing.mz=numeric(0)),control$deconvolution))
  alignment_note <- ""
  fit <- bench_time({
    x <- erah::newExp(erah::createInstrumentalTable(files))
    x <- erah::deconvolveComp(x,p)
    if (length(files)>1L) {
      ap <- do.call(erah::setAlPar,settings(list(min.spectra.cor=.85,
        max.time.dist=.15),control$alignment))
      aligned <- tryCatch(erah::alignComp(x,ap),error=identity)
      if (!inherits(aligned,"error")) x <- aligned
      else alignment_note <- conditionMessage(aligned)
    }
    x
  })
  d <- methods::slot(methods::slot(fit$value,"Data"),"FactorList")
  parsed <- lapply(seq_along(files),function(i) {
    a <- d[[i]]
    if (is.null(a) || !nrow(a)) return(records())
    group <- if ("AlignID" %in% names(a)) as.character(a$AlignID) else NA_character_
    peaks <- peak_table(sample_id(files[i]),a$ID,as.numeric(a$RT),a$Area,a[["Peak Height"]],group)
    records(peaks,lapply(a$Spectra,numeric_pairs,pattern=","),
      lapply(a$Profile,numeric_pairs,pattern=","))
  })
  combined <- combine_records(parsed,"eRah reported area","relative component profile")
  if (length(files)>1L) {
    aligned <- methods::slot(methods::slot(fit$value,"Results"),"Alignment")
    groups <- unique(stats::na.omit(combined$peaks$group))
    if (length(groups)) combined$shared <- length(groups)
    else if (nrow(aligned)) combined$shared <- nrow(aligned)
    else combined$note <- if (nzchar(alignment_note)) paste("eRah alignment:",alignment_note)
      else "No eRah alignment table was returned."
  }
  combined$analysis <- fit$seconds
  combined$workers <- workers
  saveRDS(fit$value,file.path(out,"eRah.rds"))
  combined
}

combine_records <- function(x,area_kind,profile_kind) {
  records(do.call(rbind,lapply(x,`[[`,"peaks")),
    unlist(lapply(x,`[[`,"spectra"),recursive=FALSE),
    unlist(lapply(x,`[[`,"profiles"),recursive=FALSE),area_kind,profile_kind)
}

bench_msdial <- function(files,out,ncores,control,timeout,dots) {
  exe <- executable(control$exe,"MSDIALCUI")
  format <- if (is.null(control$format)) "mzML" else control$format
  prep <- bench_time(stage_inputs(files,file.path(out,"input"),format))
  output <- file.path(out,"export")
  temp <- file.path(out,"temp")
  for (d in c(output,temp)) dir.create(d,showWarnings=FALSE)
  template <- file.path(out,"method.txt")
  if (is.null(control$method_file)) {
    run_logged(exe,c("template","gcms","-o",template,"-f"),out,
      file.path(out,"template.log"),timeout)
  } else if (!file.copy(control$method_file,template,overwrite=TRUE))
    stop("Could not copy the MS-DIAL method.",call.=FALSE)
  method <- readLines(template,warn=FALSE)
  overrides <- settings(list("Number of threads"=ncores),control$parameters)
  for (key in names(overrides)) {
    i <- which(tolower(trimws(sub(":.*$","",method)))==tolower(key))
    if (length(i)!=1L) stop("MS-DIAL method setting not found: ",key,call.=FALSE)
    method[i] <- paste0(key,": ",overrides[[key]])
  }
  writeLines(method,template)
  log <- file.path(out,"run.log")
  env <- c(TMP=temp,TEMP=temp)
  if (!is.null(control$netcdf)) {
    env <- c(env,PATH=paste(control$netcdf,Sys.getenv("PATH"),sep=.Platform$path.sep))
  }
  find_dcl <- function() {
    paths <- list.files(out,"\\.dcl$",full.names=TRUE,recursive=TRUE,ignore.case=TRUE)
    vapply(sample_id(files),function(id) {
      candidates <- paths[startsWith(basename(paths),paste0(id,"_"))]
      if (length(candidates)) utils::tail(sort(candidates),1L) else NA_character_
    },"")
  }
  run <- bench_time(run_logged(exe,c("gcms","-i",dirname(prep$value[1L]),"-o",
    output,"-m",template),out,log,timeout,env=env,
    expected=function() !anyNA(find_dcl()) && !any(grepl("Exception|Error:",
      readLines(log,warn=FALSE),ignore.case=TRUE))))
  files_tsv <- file.path(out,"msdial.tsv")
  powershell <- if (is.null(control$powershell) && nzchar(Sys.which("pwsh"))) "pwsh" else "powershell.exe"
  run_logged(executable(control$powershell,powershell),
    c("-NoProfile","-ExecutionPolicy","Bypass","-File",adapter_asset("export-msdial.ps1"),
      dirname(exe),files_tsv,unname(find_dcl())),out,file.path(out,"export.log"),timeout)
  a <- utils::read.delim(files_tsv,check.names=FALSE,na.strings="",quote="",comment.char="")
  peaks <- peak_table(a$sample,a$id,a$rt,a$area,a$height)
  ans <- records(peaks,lapply(a$spectrum,numeric_pairs,pattern=","),
    lapply(a$profile,numeric_pairs,pattern=","),"MS-DIAL model-ion area",
    "MS-DIAL model-ion profile")
  if (length(files)>1L) {
    t <- list.files(output,"\\.mzTab$",full.names=TRUE,ignore.case=TRUE)
    if (length(t)==1L) ans$shared <- sum(startsWith(readLines(t,warn=FALSE),"SML\t"))
  }
  ans$analysis <- run$seconds
  ans$preparation <- prep$seconds
  ans
}
