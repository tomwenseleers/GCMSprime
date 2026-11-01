# Run from the GCMSprime repository root; set the local input path below.
library(GCMSprime)
msp_file <- "path/to/library.msp"
n_spectra <- 10000L
db_file <- tempfile(fileext=".sqlite")

memory_seconds <- system.time(memory <- read_msp(msp_file, backend="memory",
  max_spectra=n_spectra))[["elapsed"]]
sqlite_seconds <- system.time(disk <- read_msp(msp_file, backend="sqlite",
  db_file=db_file, max_spectra=n_spectra))[["elapsed"]]
indices <- unique(c(1L, memory$n %/% 2L, memory$n))
indices <- indices[indices>0]
stopifnot(identical(spectra(memory, indices, format="list"),
                    spectra(disk, indices, format="list")))
data.frame(backend=c("memory","sqlite"), spectra=c(memory$n,disk$n),
  seconds=c(memory_seconds,sqlite_seconds),
  MB=c(as.numeric(object.size(memory)),file.info(db_file)$size)/1e6)
