utils::globalVariables(c("significant", "median_lfc"))
#' Generate a MA-plot of the analysis
#'
#' @param hits a `tibble` produced by \code{\link[pair]{extract_results}}
#'
#' @return a `ggplot2` of the distribution of the hits
#' @export
#'
#' @import utils
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_*
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' # Generate the contrast matrix
#' contrast <- limma::makeContrasts(
#'   contrasts = "ng100-ng50",
#'   levels = design
#' )
#'
#' # Normalize and log-transform the data
#' yeast <- prnn(yeast, "identifier")
#'
#' \dontrun{
#' results <- run_pipeline(yeast, design, contrast, 1000, 5, "identifier", TRUE)
#' imputation_summary <- extract_results(yeast, results, .05, 1, "fdr", "identifier")
#' plot_ma(imputation_summary)
#' }
plot_ma <- function(hits) {
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
