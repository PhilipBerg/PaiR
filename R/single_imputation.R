single_imputation <- function(data,
                              design,
                              id_col = "id") {
  conditions <- design %>%
    get_conditions()
  sd_mean <- data %>%
    tidyr::pivot_longer(dplyr::matches(conditions)) %>%
    dplyr::mutate(
      name = stringr::str_replace(name,
                                  paste0("(", conditions, ")", ".*"), "\\1"
      )
    )
  gamma_reg <- sd_mean %>%
    fit_gamma_imputation(design, id_col)
  char_cols <- data %>%
    purrr::keep(is.character)
  order <- data %>%
    colnames()
  data %>%
    prep_data_for_imputation(conditions, gamma_reg) %>%
    impute(char_cols, order)
}
