#' Function for plotting the mean-variance gamma regression
#'
#' @param data
#' @param design
#' @param id_col
#'
#' @return
#' @export
#'
#' @importFrom dplyr %>%
#' @import utils
#'
#' @examples
utils::globalVariables(c(".", "sd", "model"))
plot_gamma_regression <- function(data, design, id_col = 'id'){
  cols_to_plot <- design %>%
    get_conditions()
  precision_plot <- data %>%
    tidyr::pivot_longer(tidyr::matches(cols_to_plot)) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE)
    ) %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd)) +
    ggplot2::geom_point(size = 1/10) +
    ggplot2::geom_smooth(
      method = stats::glm,
      formula = y ~ x,
      method.args = list(family = stats::Gamma(log)),
      fullrange = TRUE
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression(hat(mu))) +
    ggplot2::ylab(expression(hat(sigma))) +
    ggplot2::ggtitle('For precision weights')
  imputation_plot <- data %>%
    tidyr::pivot_longer(tidyr::matches(cols_to_plot)) %>%
    dplyr::mutate(name = stringr::str_extract(name, cols_to_plot)) %>%
    dplyr::group_by(name, .data[[id_col]]) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    tidyr::drop_na() %>%
    ggplot2::ggplot(ggplot2::aes(mean, sd)) +
    ggplot2::geom_point(size = 1 / 10) +
    ggplot2::geom_smooth(
      method = stats::glm,
      formula = y ~ x,
      method.args = list(family = stats::Gamma(log)),
      fullrange = TRUE
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression(hat(mu))) +
    ggplot2::ylab(expression(hat(sigma))) +
    ggplot2::facet_wrap(name ~ .) +
    ggplot2::ggtitle('For imputation')
  plots <- cowplot::plot_grid(precision_plot, imputation_plot)
  title <- cowplot::ggdraw() +
    cowplot::draw_label("Mean-Variance trends", fontface = 'bold')
  cowplot::plot_grid(
    title,
    plots,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
}
