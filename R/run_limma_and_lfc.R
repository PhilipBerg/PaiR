#' Run limma and calculate LFC
#'
#' Function to calculate LFC, run limma using \code{\link[limma]{lmFit}}, \code{\link[limma]{contrasts.fit}},
#' and \code{\link[limma]{eBayes}} with the flag robust = TRUE.
#' If neither `weights` or `gamma_reg_model` is provided, then
#' \code{\link[limma]{lmFit}} will be run without precision weights.
#' If both are provided, then it will default to using the `gamma_reg_model`
#' model to produce the weights.
#'
#' @param data a `data.frame` with the samples and feature ids
#' @param design a design or model matrix produced by
#'     \code{\link[stats]{model.matrix}}
#' @param contrast_matrix a contrast matrix produced by
#'     \code{\link[limma]{makeContrasts}}
#' @param gamma_reg_model the regression model produced by
#'     `fit_gamma_weights` (see \code{\link[pair]{Mean-Variance_Gamma_Regressions}}) or any
#'      `glm` with `formula` `sd ~ mean`
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#'
#' @param weights a matrix of precision weights
#'
#' @return a `tibble` with the id_col, then one p_val_* and lfc_* for each
#'     comparison (*) in the contrast matrix
#' @export
#'
#' @importFrom rlang :=
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
#' contrasts = 'ng100-ng50',
#' levels = design
#' )
#'
#' # Normalize and log-transform the data
#' yeast <- prnn(yeast, "identifier")
#'
#' # Fit the gamma regressions
#' gamma_reg_model <- fit_gamma_weights(yeast, design, 'identifier')
#'
#' # Examplify on the non-missing data
#' yeast <- tidyr::drop_na(yeast)
#'
#' results <- run_limma_and_lfc(
#' yeast,
#' design,
#' contrast,
#' gamma_reg_model,
#' 'identifier'
#' )
run_limma_and_lfc <- function(data,
                              design,
                              contrast_matrix,
                              gamma_reg_model = NULL,
                              id_col = "id",
                              weights = NULL) {
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
  if (!is.null(gamma_reg_model)) {
    weights <- calc_weights(data, gamma_reg_model) %>%
      as.matrix()
  } else if (!is.null(weights) & (dim(weights) != dim(data))) {
    msg <- 'Incorect dimensions between data and weights.'
    incorrect_dims <- which(dim(weights) != dim(data))
    for (i in seq_along(incorrect_dims)) {
      msg[i+1] <- dplyr::case_when(
        incorrect_dims[i] == 1 ~ glue::glue('Data has {dim(data)[incorrect_dims[i]]} rows while weights have {dim(weights)[incorrect_dims[i]]}'),
        incorrect_dims[i] == 2 ~ glue::glue('Data has {dim(data)[incorrect_dims[i]]} columns while weights have {dim(weights)[incorrect_dims[i]]}')
      )
    }
    rlang::abort(
      stringr::str_flatten(msg, '\n')
    )
  }
  hits <- limma::lmFit(data, design, weights = weights) %>%
    limma::contrasts.fit(contrast_matrix) %>%
    limma::eBayes(robust = TRUE)
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
