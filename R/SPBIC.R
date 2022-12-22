#' @title  Scaled Unit Information Prior Bayesian Information Criterion
#'
#' @description Calculates Scaled Unit Information Prior Bayesian Information Criterion (SPBIC) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' SPBIC (Bollen et al., 2012) is calculated as
#'
#' \deqn{-2LL(theta) + k(1 - log(k/(beta^T(Sigma)^{-1}beta)))}
#'
#' beta and Sigma are vector and covariance matrix of regression coefficients.
#'
#' @return SPBIC measurement of the model
#'
#' @importFrom stats logLik var
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
#' SPBIC(m1)
#' SPBIC(m2)
#' SPBIC(m3)
#'
#' @references
#' Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A comparison of Bayes factor approximation methods including two new methods. Sociological Methods & Research, 41(2), 294-324.
#'
#' @export

SPBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  Fisher_M <- solve(reverse_fisher(model))
  coeff <- c(model$coefficients, var(model$residuals))
  E <- coeff%*%Fisher_M%*%coeff
  c(-2*LL + df*(1 - log(df/(E))))
}
