#' @title  Generalized Cross-Validation
#'
#' @description Calculates Generalized Cross-Validation (GCV) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' GCV (Koc and Bozdogan, 2015) is calculated as
#'
#' \deqn{RSS/(n(1 - k/n))}
#'
#' RSS is the residual sum of squares.
#'
#' @return GCV measurement of the model
#'
#' @importFrom stats logLik
#' @examples
#' x1 <- rnorm(100, 3, 2)
#' x2 <- rnorm(100, 5, 3)
#' x3 <- rnorm(100, 67, 5)
#' err <- rnorm(100, 0, 4)
#'
#' ## round so we can use it for Poisson regression
#' y <- round(3 + 2*x1 - 5*x2 + 8*x3 + err)
#'
#' m1 <- lm(y~x1 + x2 + x3)
#' m2 <- glm(y~x1 + x2 + x3, family = "gaussian")
#' m3 <- glm(y~x1 + x2 + x3, family = "poisson")
#'
#'GCV(m1)
#'GCV(m2)
#'GCV(m3)
#'
#' @references
#' Koc, E. K., & Bozdogan, H. (2015). Model selection in multivariate adaptive regression splines (MARS) using information complexity as the fitness function. Machine Learning, 101(1), 35-58.
#'
#' @export

GCV <- function(model) {
  LL <- logLik(object = model)
  n <- length(model$residuals)
  if ("glm" %in% class(model)) {
    res <- model$y - model$fitted.values
  } else {
    res <- model$residuals
  }
  MSE <- mean(res^2)
  df <- attr(LL, "df")
  MSE/(1 - df/n)^2
}
