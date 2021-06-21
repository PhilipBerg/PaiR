utils::globalVariables(c("flag", "value", "ref", "where"))
#' Filter lowly abundant features
#'
#' Function for filtering lowly abundant features.
#' By default, it uses all numerical columsn
#'
#' @param data data to filter from
#' @param target columns to base the filtering on, supports \code{\link[tidyselect]{tidyselect-package}}
#' @param percent A feature gets filtered out if it is lowly abundant or missing
#'  in `percent` columns
#' @param k Parameter for the lower limit of Tukey's fence, any value bellow this
#' will be considered an outlier
#' @param lower_limit a user defined lower limit at which a measurement is
#'   considered an outlier
#'
#' @return data with outliers removed
#' @export
#'
#' @import utils
#'
#'
#' @examples
#' # Since Tukey's fences are not ideal for raw proteomics data one could use
#' # the e.g., the tenth percentile as a indicator of lower abundance
#' filter_outliers(yeast, lower_limit = stats::quantile(yeast[-1], .1, na.rm = TRUE))
#' # We recommend normalizing the data before filtering outliers.
#' # This way we ensure that no peptides are considered outliers as an effect
#' # of a set of samples, one average, have lower quantification
#' yeast <- prnn(yeast, 'identifier')
#' filter_outliers(yeast, -1, 1, 1.5)
filter_outliers <- function(data,
                            target = NULL,
                            percent = 1,
                            k = 1.5,
                            lower_limit = NULL) {
  target <- rlang::enquo(target)
  target <- check_target(target)
  limit <- data %>%
    dplyr::select(!!target) %>%
    ncol()
  limit <- limit * percent
  data %>%
    dplyr::mutate(flag = flag_outliers(
      dplyr::across(!!target),
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
      unlist(..., TRUE, FALSE),
      .25,
      na.rm = TRUE
    ) - k * stats::IQR(
      unlist(..., TRUE, FALSE),
      na.rm = TRUE
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
