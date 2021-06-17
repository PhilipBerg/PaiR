#' Calculate Precision Weights From Gamma Regression
#'
#' \code{calc_weights} returns a data frame with all quantitative values
#' replaced by their corresponding precision weights estimated from the gamma
#' regression.
#'
#' \code{calc_weights} takes as input a data frame and a `glm` object produced
#' by \code{\link[pair]{fit_gamma_regression}} or \code{\link[pair]{fit_gamma_weights}}.
#' For all numeric columns, it predicts the standard deviation using the gamma
#' regression. It then squares and takes the reciprocal of each value to generate the
#' precision weights.
#'
#' @param data A `data.frame` who's quantitative columns should be converted to
#'  precision weights.
#' @param gamma_reg_model a `glm` object as produced by
#'  \code{\link[pair]{fit_gamma_regression}} or \code{\link[pair]{fit_gamma_weights}}.
#'  Any `glm` object with `formula` `sd ~ mean` is valid.
#'
#'
#' @return The same `data.frame` but with all quantitative values replaced by
#'  their precision weights
#' @export
#'
#' @examples
#' # Generate a design matrix for the data
#' design <- model.matrix(~0+factor(rep(1:2, each = 3)))\cr
#'
#' # Set correct colnames, this is important for fit_gamma_weights
#' colnames(design) <- paste0('ng', c(50, 100))
#'
#' # Fit the gamma regression model for the mean-variance trend
#' gamma_model <- fit_gamma_weights(yeast, design, 'identifier')\
#'
#' # Generate the weights for the yeast data
#' calc_weights(yeast, gamma_model)\
#'
#' # Note that, unless data has been log-transformed, it is likely that the
#'  regression model will to not converge
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
