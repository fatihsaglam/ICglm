#' @title Consistent Akaike's Information Criterion and Consistent Akaike's Information Criterion with Fisher Information
#'
#' @description Consistent Akaike's Information Criterion (CAIC) and Consistent Akaike's Information Criterion with Fisher Information (CAICF) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object.
#'
#' @details
#' CAIC (Bozdogan, 1987) is calculated as
#'
#' \deqn{-2LL(theta) + k(log(n) + 1)}
#'
#' CAICF (Bozdogan, 1987) as
#'
#' \deqn{-2LL(theta) + 2k + k(log(n)) + log(|F|)}
#'
#' F is the Fisher information matrix.
#'
#' @return CAIC or CAICF measurement of the model.
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
#'CAIC(m1)
#'CAIC(m2)
#'CAIC(m3)
#'CAICF(m1)
#'CAICF(m2)
#'CAICF(m3)
#'
#' @references
#'Bozdogan, H. (1987). Model selection and Akaike's information criterion (AIC): The general theory and its analytical extensions. Psychometrika, 52(3), 345-370.
#'
#' @rdname CAIC
#' @export

CAIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*(log(n) + 1))
}

#' @rdname CAIC
#' @export
#'
CAICF <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*(log(n) + 2) + log(det(solve(reverse_fisher(model)))))
}
