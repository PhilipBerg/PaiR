#' Title
#'
#' @param data
#' @param design
#' @param contrast_matrix
#' @param imputations
#' @param workers
#' @param id_col
#' @param plot_trend
#'
#' @return
#' @export
#'
#' @import utils
#'
#' @examples
utils::globalVariables(
  c("tmp_id", "imputation", "imputed_data", "limma_results")
)
run_pipeline <- function(data,
                         design,
                         contrast_matrix,
                         imputations,
                         workers,
                         id_col = "id",
                         plot_trend = F) {
  # Fit gamma models
  gamma_reg_models <- fit_gamma_regressions(data, design, id_col = id_col)
  gamma_reg_all <- gamma_reg_models$all
  if (plot_trend) {
    plot_gamma_regression(data, design, id_col = id_col)
  }
  # Generate imputation input
  col_order <- names(data)
  missing_data <- data %>%
    dplyr::filter(dplyr::if_any(where(is.numeric), is.na))
  char_cols <- missing_data %>%
    purrr::keep(is.character)
  conditions <- design %>%
    get_conditions()
  impute_nested <- missing_data %>%
    prep_data_for_imputation(conditions, gamma_reg_models$individual)
  # Generate results
  ## Non-missing data
  non_miss_result <- data %>%
    tidyr::drop_na(where(is.numeric)) %>%
    run_limma_and_lfc(design, contrast_matrix, gamma_reg_all, id_col)
  if (workers != 1) {
    cluster <- multidplyr::new_cluster(workers)
    multidplyr::cluster_library(
      cluster,
      c(
        "dplyr",
        "stringr",
        "tidyr",
        "purrr",
        "limma",
        "tibble"
      )
    )
    multidplyr::cluster_copy(
      cluster,
      c(
        "impute",
        "impute_row",
        "char_cols",
        "col_order",
        "run_limma_and_lfc",
        "design",
        "contrast_matrix",
        "gamma_reg_all",
        "calc_weights",
        "calc_lfc",
        "impute_nested",
        "non_miss_result",
        "id_col"
      )
    )
  }
  results <- tibble::tibble(
    imputation = seq_len(imputations)
  )
  if (workers != 1) {
    results <- results %>%
      multidplyr::partition(cluster)
  }
  results <- results %>%
    dplyr::mutate(
      # Run imputation
      imputed_data = purrr::map(
        imputation,
        ~ impute(impute_nested, char_cols, col_order)
      ),
      # Run limma
      limma_results = purrr::map(
        imputed_data,
        run_limma_and_lfc,
        design, contrast_matrix, gamma_reg_all, id_col
      ),
      # Bind non-missing data
      limma_results = purrr::map(
        limma_results,
        dplyr::bind_rows,
        non_miss_result
      )
    )
  if (workers != 1) {
    results <- results %>%
      dplyr::collect() %>%
      dplyr::arrange(imputation)
  }
  return(results)
}
