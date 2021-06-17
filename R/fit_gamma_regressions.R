#' Function for fitting the mean-variance gamma regression
#'
#' @param data
#' @param design
#' @param id_col
#'
#' @return
#' @export
#'
#' @importFrom dplyr %>%
#' @import utils
#'
#' @examples
utils::globalVariables(c(".", "sd", "model"))
fit_gamma_regressions <- function(data, design, id_col = "id") {
  gamma_reg <- data %>%
    fit_gamma_imputation(design, id_col)
  gamma_reg_all <- data %>%

  return(
    list(
      all = gamma_reg_all,
      individual = gamma_reg
    )
  )
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

extract_model <- function(data) {
  if (nrow(data)>1) {
    data %>%
      purrr::map(magrittr::use_series, model) %>%
      purrr::map(magrittr::extract2, 1)
  }else {
    data %>%
      magrittr::use_series(model) %>%
      magrittr::extract2(1)
  }
}
