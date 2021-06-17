#' Title
#'
#' @param hits
#'
#' @return
#' @export
#'
#' @import utils
#'
#' @examples
utils::globalVariables(c("significant", "median_lfc"))
plot_MA <- function(hits) {
  hits %>%
    dplyr::mutate(
      significant = dplyr::case_when(
        binom_p_value < .05 & median_lfc < 0 ~ "Down",
        binom_p_value < .05 & median_lfc > 0 ~ "Up",
        T ~ "Not significant"
      ),
      significant = factor(significant,
        levels = c("Not significant", "Up", "Down")
      )
    ) %>%
    dplyr::arrange(significant) %>%
    ggplot2::ggplot(
      ggplot2::aes(median_mean, median_lfc, color = significant)
    ) +
    ggplot2::geom_point(size = 1 / 2) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(comparison ~ .) +
    ggplot2::xlab(expression(hat(mu))) +
    ggplot2::ylab("|LFC|")
}
