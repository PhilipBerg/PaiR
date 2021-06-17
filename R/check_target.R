check_target <- function(target) {
  target
  target <- rlang::enquo(target)
  if (rlang::quo_is_null(target)) {
    target <- rlang::expr(where(is.numeric))
  }
  target
}
