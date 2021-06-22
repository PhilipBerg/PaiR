utils::globalVariables(
  c(
    "p_val",
    "median_mean.x",
    "median_mean.y",
    "comparison",
    "median_mean"
  )
)
#' Extract the results from running the multiple imputations
#'
#' Does p-value correction and then binomial testing on a given decision.
#' Could be a hybrid decision on both p-value and LFC. Assumes that the decision
#' for each imputation and feature can be considered a bernoulli trial. Therefore,
#' the sum of decisions from all imputations has a binomial distribution. It
#' therefore implements a right-tailed binomial testing with default null
#' hypothesis p = 0.5. It uses the right tail because significant features are
#' expected to have a larger p then the null hypothesis.
#'
#' @param data The data used to generate the results when calling
#'   \code{\link[pair]{run_pipeline}}.
#' @param results The output from \code{\link[pair]{run_pipeline}}.
#' @param alpha The alpha value to decide when a feature is significant.
#' @param abs_lfc If a LFC cut-off values should be used in addition to the
#'     alpha value.
#' @param pcor A p-value correction method, has to be one from
#'     \code{\link[stats]{p.adjust}}.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#' @param null_hyp The value for the null hypothesis of the binomial-test.
#'     Needs to be between 0 and 1.
#'
#' @return A tibble with a summary of the results. From all imputations it
#'     calculates the binomial p-value, the median: LFC, p-value (from all testing),
#'     and the mean in each comparison. It also returns a character column,
#'     comparison, that indicates what comparison the results come from and a
#'     boolean column, imputed, indicating if there was any imputed value or not.
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
#' results <- run_pipeline(yeast, design, contrast, 1000, 5, "identifier", TRUE)
#' extract_results(yeast, results, .05, 1, "fdr", "identifier")
#' }
extract_results <- function(data,
                            results,
                            alpha = 0.05,
                            abs_lfc = 1,
                            pcor = stats::p.adjust.methods,
                            id_col = "id",
                            null_hyp = 0.5) {
  pcor <- match.arg(pcor)
  # Get the number of imputations ran
  imputations <- nrow(results)
  # Count number of success
  pvals_multi_imp <- results %>%
    dplyr::select(imputation, limma_results) %>%
    tidyr::unnest(cols = limma_results)
  if (pcor != "none") {
    pvals_multi_imp <- pvals_multi_imp %>%
      dplyr::group_by(imputation) %>%
      dplyr::mutate(
        dplyr::across(dplyr::contains("p_val"), stats::p.adjust, method = pcor)
      ) %>%
      dplyr::ungroup()
  }
  pvals_multi_imp <- pvals_multi_imp %>%
    dplyr::select(-imputation) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::contains("p_val"),
        ~ . < alpha
      ),
      dplyr::across(
        dplyr::contains("lfc_"),
        ~ abs(.) > abs_lfc
      )
    ) %>%
    tidyr::pivot_longer(
      where(is.logical),
      names_to = c(".value", "comparison"),
      names_pattern = "^(p_val|lfc)_(.*)$"
    ) %>%
    dplyr::group_by(dplyr::across(where(is.character))) %>%
    dplyr::summarise(
      binom_p_value = stats::binom.test(
        x = sum(p_val & lfc),
        n = imputations,
        p = null_hyp,
        "greater"
      )$p.value,
      .groups = "drop"
    )
  imputation_summary <- results %>%
    summarise_imputations(pcor, id_col)
  comps <- imputation_summary$comparison %>%
    unique()
  all_cols <- comps %>%
    stringr::str_replace("_vs_", "|") %>%
    stringr::str_flatten("|")
  non_miss_summary <- data %>%
    tidyr::drop_na(dplyr::matches(all_cols)) %>%
    calc_comp_means(comps, id_col) %>%
    dplyr::rename(comparison = name, median_mean = value)
  complete_summary <- dplyr::left_join(
    imputation_summary,
    non_miss_summary,
    by = c(id_col, "comparison")
  ) %>%
    dplyr::mutate(
      median_mean = dplyr::coalesce(median_mean.x, median_mean.y)
    ) %>%
    dplyr::select(-tidyr::matches("\\.(x|y)$"))
  pvals_multi_imp %>%
    dplyr::left_join(
      complete_summary,
      by = c(id_col, "comparison")
    ) %>%
    dplyr::mutate(
      comparison = stringr::str_replace_all(comparison, "_", " ")
    )
}

summarise_imputations <- function(results,
                                  pcor = stats::p.adjust.methods,
                                  id_col = "id") {
  lfc <- results %>%
    dplyr::select(limma_results) %>%
    tidyr::unnest(limma_results) %>%
    dplyr::select(1, tidyr::contains("lfc")) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), stats::median)) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    dplyr::mutate(
      name = stringr::str_remove(name, "lfc_")
    ) %>%
    dplyr::rename(median_lfc = value, comparison = name)

  comps <- lfc$comparison %>%
    unique()

  mean_values <- results %>%
    dplyr::select(imputed_data) %>%
    tidyr::unnest(imputed_data) %>%
    calc_comp_means(comps, id_col) %>%
    dplyr::group_by(.data[[id_col]], name) %>%
    dplyr::summarise(median_mean = stats::median(value), .groups = "drop") %>%
    dplyr::rename(comparison = name)

  p_val <- results %>%
    dplyr::select(imputation, limma_results) %>%
    tidyr::unnest(limma_results) %>%
    dplyr::select(imputation, {{ id_col }}, dplyr::contains("p_val"))
  if (pcor != "none") {
    p_val <- p_val %>%
      dplyr::group_by(imputation) %>%
      dplyr::mutate(
        dplyr::across(dplyr::contains("p_val"), stats::p.adjust, method = pcor)
      ) %>%
      dplyr::ungroup()
  }
  p_val <- p_val %>%
    dplyr::select(-imputation) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), stats::median)) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    dplyr::mutate(
      name = stringr::str_remove(name, "p_val_")
    ) %>%
    dplyr::rename(median_p_value = value, comparison = name)

  p_val %>%
    dplyr::left_join(
      mean_values,
      by = c(id_col, "comparison")
    ) %>%
    dplyr::left_join(
      lfc,
      by = c(id_col, "comparison")
    ) %>%
    dplyr::mutate(
      imputed = !is.na(median_mean)
    )
}

calc_comp_means <- function(data, comps, id_col = "id") {
  mean_values <- data
  comp_cols <- comps %>%
    stringr::str_replace("_vs_", "|")
  for (i in seq_along(comps)) {
    mean_values <- mean_values %>%
      dplyr::mutate(
        !!comps[i] := rowMeans(dplyr::across(dplyr::matches(comp_cols[i])))
      )
  }
  mean_values %>%
    dplyr::select({{ id_col }}, all_of(comps)) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), mean)) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    dplyr::mutate(
      name = stringr::str_extract(name, comps)
    )
}
