#' Title
#'
#'
#' @param data a `data.frame` with the samples and feature ids
#' @param design a design or model matrix produced by
#'     `\code{\link[stats]{model.matrix}}`
#' @param contrast_matrix a contrast matrix produced by
#'     `\code{\link[limma]{makeContrasts}}`
#' @param gamma_reg_model the regression model produced by
#'     `\code{\link[pair]{fit_gamma_weights}}` or any `glm` with `formula`
#'     `sd ~ mean`
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#'
#'
#'
#' @return
#' @export
#'
#' @importFrom rlang :=
#'
#' @examples
run_limma_and_lfc <- function(data,
                              design,
                              contrast_matrix,
                              gamma_reg_model,
                              id_col = "id") {
  row_names <- data[[id_col]]
  condi <- design %>%
    get_conditions()
  data <- data %>%
    dplyr::select(tidyr::matches(condi)) %>%
    as.data.frame()
  rownames(data) <- row_names
  # Calc LFC
  lfc <- data %>%
    dplyr::select(tidyr::matches(condi)) %>%
    split.default(stringr::str_remove(names(.), "_\\d+$")) %>%
    purrr::map(rowMeans) %>%
    purrr::imap(~ tibble::enframe(.x, value = .y)) %>%
    purrr::reduce(dplyr::left_join, by = "name")
  calc_lfc_for <- stringr::str_replace(
    dimnames(contrast_matrix)$Contrast,
    "-", "_vs_"
  )
  lfc <- purrr::map(calc_lfc_for, calc_lfc, lfc) %>%
    purrr::map(dplyr::select, !!id_col := name, tidyr::contains("lfc")) %>%
    purrr::reduce(dplyr::left_join, by = id_col)
  # Run LIMMA
  weights <- calc_weights(data, gamma_reg_model) %>%
    as.matrix()
  hits <- limma::lmFit(data, design, weights = weights) %>%
    limma::contrasts.fit(contrast_matrix) %>%
    limma::eBayes(robust = T)
  # Extract p-value from comparisons
  hits$p.value %>%
    tibble::as_tibble(rownames = id_col) %>%
    # Bind the data together
    dplyr::left_join(lfc, by = id_col) %>%
    dplyr::rename_with(
      ~ stringr::str_replace(., "^", "p_val_") %>%
        stringr::str_replace(., "-", "_vs_"),
      tidyr::contains("-")
    )
}

calc_lfc <- function(comparison, means) {
  columns <- stringr::str_split(comparison, "_vs_")[[1]]
  means %>%
    dplyr::mutate(
      !!paste0("lfc_", comparison) := !!dplyr::sym(columns[1]) -
        !!dplyr::sym(columns[2])
    )
}
