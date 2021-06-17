estimate_loq <- function(data) {
  data %>%
    purrr::keep(is.numeric) %>%
    unlist(T, F) %>%
    {
      stats::quantile(., .25, na.rm = T) - 1.5 * stats::IQR(., na.rm = T)
    } %>%
    unname()
}
