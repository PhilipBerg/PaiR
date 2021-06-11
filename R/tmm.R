#' Trimmed m mean
#'
#' @param data data.frame containing the data to normalize
#' @param trim_M percent of fold-change values to trim
#' @param trim_A percent of means to trim
#' @param id Column of id variables
#'
#' @return dataframe with normalized values
#' @export
#' @import utils
#'
#' @examples
utils::globalVariables(c("load_size", "lfc", "A", "w"))
tmm <- function(data, id = 'id', trim_M = .3, trim_A = .05, load_info = FALSE, target = NULL, reference_sample = NULL){
  if(is.null(target)){
    pivot_cols <- rlang::expr(where(is.numeric))
  }else{
    pivot_cols <- rlang::expr(all_of(target))
  }
  data_filtered <- data %>%
    tidyr::drop_na()
  loading_sizes <- data_filtered %>%
    dplyr::select(!!pivot_cols) %>%
    colSums()
  if(is.null(reference_sample)){
  reference_sample <- data_filtered %>%
    purrr::keep(is.numeric) %>%
    purrr::map_dbl(
      ~sd(.x)/mean(.x)
    )
  }
  reference_sample <- names(reference_sample)[reference_sample == min(reference_sample)]
  reference_loading_size <- loading_sizes[names(loading_sizes) == reference_sample]
  reference_sample <- rlang::sym(reference_sample)
  loading_sizes <- loading_sizes %>%
    tibble::enframe(name = 'sample', value = 'load_size')
  scaling_factors <- data_filtered %>%
    purrr::keep(is.numeric) %>%
    tidyr::pivot_longer(-!!reference_sample, names_to = 'sample') %>%
    dplyr::left_join(loading_sizes, by = "sample") %>%
    dplyr::mutate(
      w = (load_size - value)/(load_size*value) + (reference_loading_size - !!reference_sample)/(reference_loading_size*!!reference_sample),
      lfc = log2((value/load_size)/(!!reference_sample/reference_loading_size)),
      A = .5*log2((value/load_size)*(!!reference_sample/reference_loading_size))
    ) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(
      dplyr::between(lfc, stats::quantile(lfc, trim_M), stats::quantile(lfc, 1 - trim_M)) &
        dplyr::between(A, stats::quantile(A, trim_A), stats::quantile(A, 1 - trim_A))
    ) %>%
    dplyr::summarise(
      tmm_factor = 2^(sum(lfc*w)/sum(w))
    ) %>%
    dplyr::left_join(loading_sizes, by = "sample") %>%
    dplyr::add_row(sample = as.character(reference_sample), tmm_factor = 1, load_size = reference_loading_size)
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <- data[scaling_factors$sample[i]]/scaling_factors$tmm_factor[i]
  }
  data <- data %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), log2)
    )
  if(load_info){
    return(
      list(
        data = data,
        scaling_factors = scaling_factors
      )
    )
  }else{
    return(data)
  }
}
