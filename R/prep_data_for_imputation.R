prep_data_for_imputation <- function(data, conditions, gamma_reg_imputation) {
  LOQ <- data %>%
    estimate_loq()
  data %>%
    purrr::keep(is.numeric) %>%
    split.default(stringr::str_extract(names(.), conditions))
  purrr::imap(impute_nest, gamma_reg_imputation, LOQ)
}
