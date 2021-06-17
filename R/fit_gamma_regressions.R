#' Function for Fitting the Mean-Variance Gamma Regression Models
#'
#'`fit_gamma_regressions` calls both \code{\link[pair]{fit_gamma_regression}}
#'    and \code{\link[pair]{fit_gamma_weights}} and returns a list containing
#'    the models for imputation in `$imputation` and the weights in `$weights`.
#'
#' @param data a `data.frame` to calculate the weights for
#' @param design a design or model matrix
#' @param id_col a character for the name of the column containing the
#'     name of the features in data (e.g., peptides, proteins, etc.)
#'
#' @return a list with the gamma regressions for imputation in
#'         `$imputation$<condition>` and the gamma regression for generating
#'         precision weights in `$weights`.
#' @export
#'
#' @importFrom dplyr %>%
#' @import utils
#'
#' @examples
utils::globalVariables(c(".", "sd", "model"))
fit_gamma_regressions <- function(data, design, id_col = "id") {
  gamma_reg_imp <- data %>%
    fit_gamma_imputation(design, id_col)
  gamma_reg_w <- data %>%
    fit_gamma_weights(design, id_col)
  return(
    list(
      weights = gamma_reg_w,
      imputation = gamma_reg_imp
    )
  )
}

#' Title
#'
#' @param data
#' @param design
#' @param id_col
#'
#' @return
#' @export
#'
#' @examples
fit_gamma_imputation <- function(data, design, id_col = "id") {
  data %>%
    prep_data_for_gamma_imputation_regression(design, id_col) %>%
    dplyr::group_by(name) %>%
    fit_gamma_model() %>%
    split.data.frame(.$name) %>%
    extract_model()
}

#' Title
#'
#' @param data
#' @param design
#' @param id_col
#'
#' @return
#' @export
#'
#' @examples
fit_gamma_weights <- function(data, design, id_col = "id"){
  data %>%
    prep_data_for_gamma_weight_regression(design, id_col) %>%
    fit_gamma_model() %>%
    extract_model()
}

prep_data_for_gamma_weight_regression <- function(data, design, id_col = "id") {
  data %>%
    pivot_data_for_gamma_regression(design) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    calc_mean_sd_trend()
}

calc_mean_sd_trend <- function(data) {
  data %>%
    dplyr::summarise(
      mean = mean(value, na.rm = T),
      sd = stats::sd(value, na.rm = T),
      .groups = "drop"
    ) %>%
    dplyr::filter(sd > 0)
}

fit_gamma_model <- function(data) {
  data %>%
    dplyr::summarise(
      model = list(stats::glm(sd ~ mean, stats::Gamma(log))),
      .groups = "drop"
    )
}

prep_data_for_gamma_imputation_regression <- function(data,
                                                      design,
                                                      id_col = "id") {
  conditions <- design %>%
    get_conditions()
  data %>%
    pivot_data_for_gamma_regression(design) %>%
    dplyr::mutate(
      name = stringr::str_replace(
        name,
        paste0("(", conditions, ")", ".*"), "\\1"
      )
    ) %>%
    dplyr::group_by(name, .data[[id_col]]) %>%
    dplyr::filter(sum(!is.na(value)) >= 2) %>%
    calc_mean_sd_trend()
}

pivot_data_for_gamma_regression <- function(data, design) {
  conditions <- design %>%
    get_conditions()
  data %>%
    tidyr::pivot_longer(dplyr::matches(conditions))
}

extract_model <- function(data) {
  if (is.data.frame(data)) {
    data %>%
      magrittr::use_series(model) %>%
      magrittr::extract2(1)
  }else if(is.list(data)){
    data %>%
      purrr::map(magrittr::use_series, model) %>%
      purrr::map(magrittr::extract2, 1)
  }else {
    error_message <- paste('Data of class',
                           paste0(class(data), collapse = ', '),
                           'not suppported.'
    )
    rlang::abort(error_message)
  }
}
