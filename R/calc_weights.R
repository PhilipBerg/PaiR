#' Title
#'
#' @param data
#' @param gamma_reg_model
#'
#' @return
#' @export
#'
#' @examples
calc_weights <- function(data, gamma_reg_model) {
  data %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ stats::predict.glm(
          gamma_reg_model,
          newdata = data.frame(mean = .),
          type = "response"
        )
      ),
      dplyr::across(where(is.numeric), ~ .^(-2))
    )
}
