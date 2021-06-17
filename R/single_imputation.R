single_imputation <- function(
  data,
  method = c('gamma', 'quantile', 'cyclicloess'),
  id_col = 'id',
  design = NULL,
  contrast = NULL
){
  method <- match.arg(method)
    if(method == 'gamma'){
      sd_mean <- data %>%
        tidyr::pivot_longer(dplyr::matches(conditions)) %>%
        dplyr::mutate(
          name = stringr::str_replace(name, paste0("(", conditions, ")", ".*"), "\\1")
        )
    }
}



