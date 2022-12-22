#' @title  Reverse Fisher Matrix
#'
#' @description This function allows you to calculate Fisher Information Matrix using "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' Calculates Fisher Information Matrix using "lm" and "glm" objects. It uses
#'
#' @return a booster object with below components.
#'  \item{n_train}{Number of cases in the input dataset.}
#'
#' @importFrom Matrix bdiag
#'
#' @export

reverse_fisher <- function(model){
  cov_coef <- vcov(model)
  n <- length(model$residuals)
  lower_left <- 2*var(model$residuals)^2/n
  as.matrix(Matrix::bdiag(cov_coef, lower_left))
}
