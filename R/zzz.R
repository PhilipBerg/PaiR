# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.pair <- list(
#     pair.default_id_col = 'id'
#   )
#   toset <- !(names(op.pair) %in% names(op))
#   if(any(toset)) options(op.pair[toset])
#
#   invisible()
# }
