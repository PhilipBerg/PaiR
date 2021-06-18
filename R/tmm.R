utils::globalVariables(c("load_size", "lfc", "A", "w"))
#' Normalization by Trimmed m Means
#'
#' @param data data.frame containing the data to normalize
#' @param trim_M percent of fold-change values to trim
#' @param trim_A percent of means to trim
#' @param log Return log2 transformed values?
#' @param load_info Return loading info?
#' @param target Specify a subset of columns to normalize
#' @param reference_sample  Specify a reference sample to normalize to
#'
#' @return data frame with normalized values
#' @export
#' @import utils
#'
#' @examples
tmm <- function(data,
                trim_M = .3,
                trim_A = .05,
                log = TRUE,
                load_info = FALSE,
                target = NULL,
                reference_sample = NULL) {
  data_filtered <- data %>%
    tidyr::drop_na()
  if (is.null(reference_sample) & is.null(target)) {
    target_cols <- rlang::expr(where(is.numeric))
    reference_sample <- calc_cv(data_filtered, target_cols)
    reference_sample <-
      names(reference_sample)[reference_sample == min(reference_sample)]
  } else if (is.null(reference_sample)) {
    target_cols <- rlang::expr(where(is.numeric))
    reference_sample <- calc_cv(data_filtered, target_cols)
    reference_sample <-
      names(reference_sample)[reference_sample == min(reference_sample)]
    target_cols <- c(target, reference_sample)
  } else if (is.null(target)) {
    target_cols <- rlang::expr(where(is.numeric))
  } else {
    target_cols <- c(target, reference_sample)
  }
  loading_sizes <- calc_loading_size(data_filtered, target_cols)
  reference_loading_size <-
    loading_sizes$load_size[loading_sizes$sample == reference_sample]
  reference_sample <- rlang::sym(reference_sample)
  scaling_factors <- data_filtered %>%
    dplyr::select(!!target_cols) %>%
    tidyr::pivot_longer(-!!reference_sample, names_to = "sample") %>%
    dplyr::left_join(loading_sizes, by = "sample") %>%
    dplyr::mutate(
      w = (load_size - value) / (load_size * value) +
        (reference_loading_size - !!reference_sample) /
          (reference_loading_size * !!reference_sample),
      lfc = log2((value / load_size) /
        (!!reference_sample / reference_loading_size)),
      A = .5 * log2((value / load_size) *
        (!!reference_sample / reference_loading_size))
    ) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(
      dplyr::between(
        lfc, stats::quantile(lfc, trim_M), stats::quantile(lfc, 1 - trim_M)
      ) &
        dplyr::between(
          A, stats::quantile(A, trim_A), stats::quantile(A, 1 - trim_A)
        )
    ) %>%
    dplyr::summarise(
      tmm_factor = 2^(sum(lfc * w) / sum(w))
    ) %>%
    dplyr::left_join(loading_sizes, by = "sample") %>%
    dplyr::add_row(
      sample = as.character(reference_sample),
      tmm_factor = 1,
      load_size = reference_loading_size
    )
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <- data[scaling_factors$sample[i]] /
      scaling_factors$tmm_factor[i]
  }
  if (log) {
    data <- data %>%
      dplyr::mutate(
        dplyr::across(!!target_cols, log2)
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

calc_cv <- function(data, targets) {
  data %>%
    dplyr::select(!!targets) %>%
    purrr::map_dbl(
      ~ sd(.x) / mean(.x)
    )
}

calc_loading_size <- function(data, targets) {
  data %>%
    dplyr::select(!!targets) %>%
    colSums() %>%
    tibble::enframe(name = "sample", value = "load_size")
}
