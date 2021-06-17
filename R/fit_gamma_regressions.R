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
fit_gamma_regressions <- function(data, design, id_col = 'id'){
  conditions <- colnames(design) %>%
    stringr::str_flatten("|")
  sd_mean <- data %>%
    tidyr::pivot_longer(dplyr::matches(conditions)) %>%
    dplyr::mutate(
      name = stringr::str_replace(name, paste0("(", conditions, ")", ".*"), "\\1")
    )
  gamma_reg <- sd_mean %>%
    dplyr::group_by(name, .data[[id_col]]) %>%
    dplyr::filter(sum(!is.na(value)) >= 2) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = T),
      sd = stats::sd(value, na.rm = T),
      .groups = "drop"
    ) %>%
    dplyr::filter(sd>0) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      model = list(stats::glm(sd~ mean, stats::Gamma(log))),
      .groups = "drop"
    ) %>%
    split.data.frame(.$name) %>%
    purrr::map(magrittr::use_series, model) %>%
    purrr::map(magrittr::extract2, 1)
  gamma_reg_all <- sd_mean %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = T),
      sd = stats::sd(value, na.rm = T),
      .groups = "drop"
    ) %>%
    dplyr::summarise(
      model = list(stats::glm(sd ~ mean, stats::Gamma(log))), .groups = "drop"
    ) %>%
    magrittr::use_series(model) %>%
    magrittr::extract2(1)
  return(
    list(
      all = gamma_reg_all,
      individual = gamma_reg
    )
  )
}

fit_gamma_imputation <- funtion(sd_mean){

}
