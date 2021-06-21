utils::globalVariables(c(".", "sd", "model"))
#' Function for plotting the mean-variance gamma regression
#'
#'
#' @param data The data to use for producing the plots
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param id_col A character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#'
#' @return a `ggplot2` object with the
#' @export
#'
#' @import utils
#'
#' @examples
plot_gamma_regression <- function(data, design, id_col = "id") {
  precision_plot <- data %>%
    prep_data_for_gamma_weight_regression(design, id_col) %>%
    tidyr::drop_na() %>%
    plot_mean_sd_trend() +
    ggplot2::ggtitle("For precision weights")
  imputation_plot <- data %>%
    prep_data_for_gamma_imputation_regression(design, id_col) %>%
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

plot_mean_sd_trend <- function(data) {
  data %>%
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
