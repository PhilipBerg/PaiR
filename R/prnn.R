#' Normalize data to a pseudo-reference
#'
#' @param data data.frame
#' @param id character
#' @param load_info logical
#'
#' @return data.frame
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @import utils
#'
#' @examples
utils::globalVariables("where")
prnn <- function(data, id = 'id', load_info = FALSE){
  value <- ref <- NULL
  data_filtered <- data %>%
    tidyr::drop_na()
  loading_sizes <- data_filtered %>%
    purrr::keep(is.numeric) %>%
    colSums() %>%
    tibble::enframe(name = 'sample', value = 'load_size')
  pseudo_reference <- data_filtered %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::summarise(
      ref = prod(value^(1/dplyr::n()))
    )
  scaling_factors <- data_filtered %>%
    tidyr::pivot_longer(where(is.numeric), names_to = 'sample') %>%
    dplyr::left_join(pseudo_reference, by = id) %>%
    dplyr::mutate(
      value = value/ref
    ) %>%
    dplyr::select(-ref) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(rle_factor = stats::median(value)) %>%
    dplyr::left_join(loading_sizes, by = 'sample')
  for (i in seq_len(nrow(scaling_factors))) {
    data[scaling_factors$sample[i]] <- data[scaling_factors$sample[i]]/scaling_factors$rle_factor[i]
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
