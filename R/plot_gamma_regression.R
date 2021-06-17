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
plot_gamma_regression <- function(data, design, id_col = "id") {
  sd_mean <- data %>%
    pivot_data_for_gamma_regression(design)
  precision_plot <- sd_mean %>%
    dplyr::group_by(.data[[id_col]]) %>%
    calc_mean_sd_trend() %>%
    tidyr::drop_na() %>%
    plot_mean_sd_trend() +
    ggplot2::ggtitle("For precision weights")
  imputation_plot <- sd_mean %>%
    dplyr::mutate(name = stringr::str_extract(name, cols_to_plot)) %>%
    dplyr::group_by(name, .data[[id_col]]) %>%
    calc_mean_sd_trend() %>%
    tidyr::drop_na() %>%
    plot_mean_sd_trend() +
    ggplot2::facet_wrap(name ~ .) +
    ggplot2::ggtitle("For imputation")
  plots <- cowplot::plot_grid(precision_plot, imputation_plot)
  title <- cowplot::ggdraw() +
    cowplot::draw_label("Mean-Variance trends", fontface = "bold")
  cowplot::plot_grid(
    title,
    plots,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
}

plot_mean_sd_trend <- function(data){
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
    ggplot2::ylab(expression(hat(sigma)))
}
