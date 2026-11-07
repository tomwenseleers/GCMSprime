amdis_field <- function(header, field, pattern="[0-9.eE+-]+") {
  key <- paste0("\\|",field,"(",pattern,")(?=\\||$)")
  pos <- regexpr(key,header,perl=TRUE)
  if (pos[1L]<0) return(NA_character_)
  sub(key,"\\1",regmatches(header,pos),perl=TRUE)
}

read_elu <- function(path, sample, clock) {
  lines <- readLines(path,warn=FALSE,encoding="latin1")
  starts <- which(grepl("^NAME:",lines,ignore.case=TRUE))
  peaks <- peak_table()
  spectra <- profiles <- vector("list",length(starts))
  if (length(starts)) {
    header <- lines[starts]
    id <- vapply(header,amdis_field,"",field="CN")
    rt <- as.numeric(vapply(header,amdis_field,"",field="RT"))
    peaks <- peak_table(sample,id,rt)
    peaks$amdis_base_ion_signal <- as.numeric(vapply(header,amdis_field,"",field="AM"))
    stops <- c(starts[-1L]-1L,length(lines))
    for (j in seq_along(starts)) {
      body <- lines[seq.int(starts[j]+1L,stops[j])]
      np <- grep("^NUM PEAKS:",body,ignore.case=TRUE)[1L]
      if (!is.na(np) && np<length(body)) {
        fragment <- paste(body[seq.int(np+1L,length(body))],collapse=" ")
        spectra[[j]] <- numeric_pairs(fragment,",[[:space:]]*")
      }
      if (is.na(np) || np<=2L || !identical(trimws(body[1L]),"RE")) next
      y <- scan(text=paste(body[2L:(np-1L)],collapse=" "),quiet=TRUE)
      bounds <- amdis_field(header[j],"FR","[0-9]+-[0-9]+")
      if (!length(y) || is.na(bounds)) next
      bounds <- as.integer(strsplit(bounds,"-",fixed=TRUE)[[1L]])
      i <- seq.int(bounds[1L],bounds[2L])
      y <- pmax(y,0)
      peaks$area[j] <- sum(y)
      peaks$height[j] <- max(y)
      # Calibrate 0/1-based scan indexing against the same record's RT.
      scan <- as.integer(amdis_field(header[j],"SC"))
      if (!is.na(scan) && scan>=1L && scan<length(clock) &&
          abs(clock[scan+1L]-rt[j])<abs(clock[scan]-rt[j])) i <- i+1L
      if (length(i)==length(y) && all(i>=1L & i<=length(clock)))
        profiles[[j]] <- cbind(clock[i],y)
    }
  }
  if (!"amdis_base_ion_signal" %in% names(peaks)) peaks$amdis_base_ion_signal <- numeric()
  records(peaks,spectra,profiles,"integrated counts","AMDIS reconstructed TIC")
}

bench_amdis <- function(files,out,ncores,control,timeout,dots) {
  if (!is.null(control$elu)) {
    if (length(control$elu)!=length(files)) stop("Supply one .ELU per input CDF.",call.=FALSE)
    parts <- lapply(seq_along(files),function(i)
      read_elu(control$elu[i],sample_id(files[i]),cdf_clock(files[i])))
    ans <- combine_records(parts,"integrated counts","AMDIS reconstructed TIC")
    ans$analysis <- NA_real_
    return(ans)
  }
  if (.Platform$OS.type!="windows")
    stop("AMDIS execution requires Windows. Existing ELU exports can be imported on any platform.",
      call.=FALSE)
  exe <- executable(control$exe,"AMDIS_32.exe")
  prep <- bench_time(stage_inputs(files,tempfile("input-",tmpdir=out),"mzXML"))
  started <- proc.time()[["elapsed"]]
  elu <- vapply(prep$value,function(path) {
    subdir <- dirname(path)
    # AMDIS 2.71 writes uppercase output names and requires a serial queue.
    existing <- function() {
      all <- list.files(subdir,"\\.elu$",full.names=TRUE,ignore.case=TRUE)
      all[tolower(sample_id(all))==tolower(sample_id(path))]
    }
    if (length(existing())) stop("AMDIS exports already exist in ",subdir,
      ". Use a fresh output directory.",call.=FALSE)
    run_logged(exe,c(normalizePath(path,winslash="\\",mustWork=TRUE),"/A","/GD","/E"),subdir,
      file.path(out,paste0(sample_id(path),".log")),timeout,
      expected=function() length(existing())==1L)
    existing()[1L]
  },"")
  elapsed <- proc.time()[["elapsed"]]-started
  parts <- lapply(seq_along(files),function(i) read_elu(elu[i],sample_id(files[i]),cdf_clock(files[i])))
  ans <- combine_records(parts,"integrated counts","AMDIS reconstructed TIC")
  ans$analysis <- elapsed
  ans$preparation <- prep$seconds
  ans$workers <- 1L
  ans
}
