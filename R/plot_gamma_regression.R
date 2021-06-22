utils::globalVariables(c(".", "sd", "model"))
#' Function for plotting the mean-variance gamma regression
#'
#'
#' @param data The data to use for producing the plots
#' @param design A design matrix as produced by \code{\link[stats]{model.matrix}}.
#' @param id_col A character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#'
#' @return a `ggplot` object with the mean-variance trend used for the precision
#' weights on the left side, and the trend lines used for imputation on the
#' right side.
#' @export
#'
#' @import utils
#'
#' @examples
#' # Produce a design matrix
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Normalize and log transform the data
#' yeast <- prnn(yeast, "identifier")
#'
#' # Generate the plots
#' plot_gamma_regression(yeast, design, "identifier")
plot_gamma_regression <- function(data, design, id_col = "id") {
  precision_plot <- data %>%
    plot_gamma_precision(design, id_col)
  imputation_plot <- data %>%
    plot_gamma_imputation(design, id_col)
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

plot_gamma_precision <- function(data, design, id_col) {
  data %>%
    prep_data_for_gamma_weight_regression(design, id_col) %>%
    tidyr::drop_na() %>%
    plot_mean_sd_trend() +
    ggplot2::ggtitle("For precision weights")
}

plot_gamma_imputation <- function(data, design, id_col) {
  data %>%
    prep_data_for_gamma_imputation_regression(design, id_col) %>%
    tidyr::drop_na() %>%
    plot_mean_sd_trend() +
    ggplot2::facet_wrap(name ~ .) +
    ggplot2::ggtitle("For imputation")
}
