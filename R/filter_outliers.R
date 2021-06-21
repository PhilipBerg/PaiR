utils::globalVariables(c("flag", "value", "ref", "all_of"))
#' Filter lowly abundant features
#'
#' @param data data to filter from
#' @param target columns to base the filtering on
#' @param percent A feature gets filtered out if it is lowly abundant in `percent`
#'  columns
#' @param k Parameter for the lower limit of Tukey's fence, any value bellow this
#' will be considered an outlier
#'
#' @return data with outliers removed
#' @export
#'
#' @import utils
#'
#' @examples
filter_outliers <- function(data,
                            target = NULL,
                            percent = 1,
                            k = 1.5,
                            lower_limit = NULL) {
  target_cols <- check_target(target)
  limit <- data %>%
    dplyr::select(!!target_cols) %>%
    ncol()
  limit <- limit * percent
  data %>%
    dplyr::mutate(flag = flag_outliers(
      dplyr::across(!!target_cols),
      k = k,
      lower_limit = lower_limit
    )
    ) %>%
    dplyr::filter(flag < limit) %>%
    dplyr::select(-flag)
}

flag_outliers <- function(..., k = 1.5, lower_limit = NULL) {
  x <- tibble::as_tibble(...)
  if (is.null(lower_limit)) {
    lower_limit <- stats::quantile(
      unlist(..., T, F),
      .25,
      na.rm = T
    ) - k * stats::IQR(
      unlist(..., T, F),
      na.rm = T
    )
  }
  x %>%
    dplyr::mutate(
      dplyr::across(dplyr::everything(), tidyr::replace_na, -Inf),
      dplyr::across(dplyr::everything(), dplyr::between, lower_limit, Inf),
      dplyr::across(dplyr::everything(), magrittr::not)
    ) %>%
    rowSums()
}
