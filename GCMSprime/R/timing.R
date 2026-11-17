# Nested stage measurements must not force a full collection per fit.
GCMS_time <- function(expr) base::system.time(expr, gcFirst = FALSE)
