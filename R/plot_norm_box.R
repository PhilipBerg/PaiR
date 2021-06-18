utils::globalVariables(c("where", "value", "name", "method", "med", "condition"))
#' Generates a box-plot for `log2` transformed raw, tmm-, and psn- normalized
#'     values.
#'
#' This function can be used to produce a visual aid to for selection
#'     normalization method. Ideally, after the data is normalized all
#'     boxplots should have their median aligned close to the global trend line.
#'
#' @inheritParams tmm
#' @inheritParams prnn
#' @param data
#' @param id_col
#' @param trim_M
#' @param trim_A
#' @param target
#' @param reference_sample Specify a reference sample to normalize to in the
#'     \code{\link[pair]{tmm}} method
#'
#' @return ggplot with samples on the x-axis and observed values on the y-axis,
#'     different colors correspond to the raw or normalized data. Lines
#'     corresponds to the global median across all samples.
#' @export
#'
#' @examples
plot_norm_box <- function(data,
                          id_col = "id",
                          trim_M = .3,
                          trim_A = .05,
                          target = NULL,
                          reference_sample = NULL) {
  prnn <- data %>%
    prnn(id_col = id_col, load_info = F) %>%
    dplyr::rename_with(~ paste0(., "_prnn"), where(is.numeric))
  tmm <- data %>%
    tmm(
      id_col = id_col,
      load_info = F,
      trim_M = trim_M,
      trim_A = trim_A,
      target = target,
      reference_sample = reference_sample
    ) %>%
    dplyr::rename_with(~ paste0(., "_tmm"), where(is.numeric))
  data %>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), log2)
    ) %>%
    dplyr::rename_with(~ paste0(., "_raw"), where(is.numeric)) %>%
    dplyr::left_join(prnn, by = id_col) %>%
    dplyr::left_join(tmm, by = id_col) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    tidyr::extract(name, c("condition", "method"), "^(.*)_(.*)$") %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(
      med = stats::median(value, na.rm = T)
    ) %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(condition, value, fill = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = med, color = method)) +
    ggplot2::theme_bw()
}
