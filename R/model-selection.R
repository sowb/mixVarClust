#' Model selection criteria
#' @description Functions for computing the implemented model selection criteria (BIC, AIC).
#'@describeIn model-selection Compute the BIC criterion
#'
#'@param loglik a numeric value. The log-liklihood of the model.
#'@param nb_params a numeric value. The number of parameters in the model.
#'@param nb_obs a numeric value. The number of observations in the dataset.
#'
#'@export
criterionBIC <- function(loglik, nb_params, nb_obs) {
    -2 * loglik + nb_params * log(nb_obs)
}

#'@describeIn model-selection Compute the AIC criterion
#'
#'@export
criterionAIC <- function(loglik, nb_params) {
    -2 * loglik + (2 * nb_params)
}

