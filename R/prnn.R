utils::globalVariables(c("where", "value", "ref", "all_of"))
#' Normalize data to a pseudo-reference
#'
#' @param data data.frame
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#' @param load_info logical
#' @param log boolean variable indicating if the data should be log transformed
#'     after normalization
#' @param target target columns to normalize, supports
#'     \code{\link[tidyselect]{tidyselect-package}} syntax. By default, all numerical
#'     columns will be used in the normalization if not specified.
#'
#' @return data.frame
#' @export
#' @importFrom dplyr .data
#' @import utils
#'
#' @examples
#' yeast_prnn <- prnn(yeast, 'identifier')
prnn <- function(data,
                 id_col = "id",
                 log = TRUE,
                 load_info = FALSE,
                 target = NULL) {
  target <- rlang::enquo(target)
  target <- check_target(target)
  data_filtered <- data %>%
    tidyr::drop_na(!!target)
  loading_sizes <- calc_loading_size(data_filtered, target)
  pseudo_reference <- data_filtered %>%
    tidyr::pivot_longer(!!target) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      ref = prod(value^(1 / dplyr::n()))
    )
  scaling_factors <- data_filtered %>%
    tidyr::pivot_longer(!!target, names_to = "sample") %>%
    dplyr::left_join(pseudo_reference, by = id_col) %>%
    dplyr::mutate(
      value = value / ref
    ) %>%
    dplyr::select(-ref) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(rle_factor = stats::median(value)) %>%
    dplyr::left_join(loading_sizes, by = "sample")
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <-
      data[scaling_factors$sample[i]] / scaling_factors$rle_factor[i]
  }
  if (log) {
    data <- data %>%
      dplyr::mutate(
        dplyr::across(!!target, log2)
      )
  }
  if (load_info) {
    return(
      list(
        data = data,
        scaling_factors = scaling_factors
      )
    )
  } else {
    return(data)
  }
}
