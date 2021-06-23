utils::globalVariables(c("significant", "median_lfc"))
#' Generate a MA-plot of the analysis.
#'
#' @param hits a `tibble` produced by \code{\link[pair]{extract_results}}
#' @param data the `data.frame` used in the analysis.
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[pair]{extract_results}}.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[pair]{extract_results}}.
#' @param alpha The alpha cut-off for considering a p-value significant.
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[pair]{extract_results}}.
#' @param abs_lfc If a LFC threshold should also be used in the decision.
#'     This should always be `NULL` if `hits` was produced by
#'     \code{\link[pair]{extract_results}}.
#'
#' @return a `ggplot2` of the distribution of the hits.
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
#' yeast <- psrn(yeast, "identifier")
#' \dontrun{
#'
#' results <- run_pipeline(yeast, design, contrast, 1000, 5, "identifier", TRUE)
#' imputation_summary <- extract_results(yeast, results, .05, 1, "fdr", "identifier")
#' plot_ma(imputation_summary)
#' }
plot_ma <- function(hits, data = NULL, id_col = NULL, alpha = NULL, abs_lfc = NULL) {
  if (is.null(data)) {
    hits <- hits %>%
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
      dplyr::rename(value = median_mean)
  }else {
    tidy_results <- hits %>%
      tidyr::pivot_longer(
        where(is.numeric),
        names_to = c(".value", "comparison"),
        names_pattern = "^(p_val|lfc)_(.*)$"
      )
    comps <- tidy_results$comparison %>%
      unique()
    means <- data %>%
      calc_comp_means(comps, id_col)
    hits <- tidy_results %>%
      dplyr::left_join(
        means, by = c(id_col, 'comparison' = 'name')
      ) %>%
      dplyr::mutate(
        significant = dplyr::case_when(
          p_val < alpha & lfc < abs_lfc ~ "Down",
          p_val < alpha & lfc > -abs_lfc ~ "Up",
          T ~ "Not significant"
        ),
        significant = factor(significant,
                             levels = c("Not significant", "Up", "Down")
        )
      ) %>%
      dplyr::rename(median_lfc = lfc)
  }
  hits %>%
    dplyr::arrange(significant) %>%
    ggplot2::ggplot(
      ggplot2::aes(value, median_lfc, color = significant)
    ) +
    ggplot2::geom_point(size = 1 / 2) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(comparison ~ .) +
    ggplot2::xlab(expression(hat(mu))) +
    ggplot2::ylab("|LFC|")
}
