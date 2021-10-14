#' Single imputation
#'
#' Performs a single imputation run and returns the data with NA values replaced
#' by imputed values.
#'
#' @param data a `data.frame` to perform the imputation on, missing values should
#' be `NA`.
#' @param design a design or model matrix as produced by
#'  \code{\link[stats]{model.matrix}} with column names corresponding to the
#'  different conditions.
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.).
#'
#' @return a `data.frame` with `NA` values replaced by imputed values.
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
#'
#' # Set correct colnames, this is important for fit_gamma_weights
#' colnames(design) <- paste0("ng", c(50, 100))
#'
#' yeast %>%
#'   # Normalize and log-transform the data
#'   psrn("identifier") %>%
#'   # Run the imputation
#'   single_imputation(design, "identifier")
single_imputation <- function(data,
                              design,
                              id_col = "id") {
  conditions <- design %>%
    get_conditions()
  gamma_reg <- data %>%
    fit_gamma_imputation(design, id_col)
  char_cols <- data %>%
    purrr::keep(is.character)
  order <- data %>%
    colnames()
  data %>%
    prep_data_for_imputation(conditions, gamma_reg) %>%
    impute(char_cols, order)
}


impute <- function(data, char_cols, order) {
  data %>%
    purrr::map(
      dplyr::mutate,
      data = purrr::pmap(list(mean, sd, data), impute_row)
    ) %>%
    purrr::map(dplyr::select, data) %>%
    purrr::map(tidyr::unnest, data) %>%
    dplyr::bind_cols(char_cols) %>%
    dplyr::select(dplyr::all_of(order))
}

impute_nest <- function(data, condition, gamma_reg_model, LOQ) {
  if (anyNA(data)) {
    data <- data %>%
      tibble::rownames_to_column(var = "tmp_id") %>%
      dplyr::mutate(
        mean = rowMeans(dplyr::across(where(is.numeric)), na.rm = T),
        mean = tidyr::replace_na(mean, LOQ)
      )
    data[["sd"]] <- stats::predict(gamma_reg_model[[condition]], data, type = "response")
    data %>%
      tidyr::nest(data = -c(mean, sd, tmp_id))
  } else {
    return(data)
  }
}

impute_row <- function(mean, sd, data) {
  if (!anyNA(data)) {
    return(data)
  } else {
    data <- as.data.frame(data)
    data[is.na(data)] <- stats::rnorm(n = sum(is.na(data)), mean = mean, sd = sd)
    return(data)
  }
}


prep_data_for_imputation <- function(data, conditions, gamma_reg_imputation) {
  LOQ <- data %>%
    estimate_loq()
  data %>%
    purrr::keep(is.numeric) %>%
    split.default(stringr::str_extract(names(.), conditions)) %>%
    purrr::imap(impute_nest, gamma_reg_imputation, LOQ)
}

estimate_loq <- function(data) {
  data %>%
    purrr::keep(is.numeric) %>%
    unlist(T, F) %>%
    {
      stats::quantile(., .25, na.rm = T) - 1.5 * stats::IQR(., na.rm = T)
    } %>%
    unname()
}
