#' @title  Haughton Bayesian information criterion
#'
#' @description Calculates Haughton Bayesian information criterion (HBIC) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' HBIC (Bollen et al., 2014) is calculated as
#'
#' \deqn{-2LL(theta) + klog(n/(2pi))}
#'
#' @return HBIC measurement of the model
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
#' HBIC(m1)
#' HBIC(m2)
#' HBIC(m3)
#'
#' @references
#' Bollen, K. A., Harden, J. J., Ray, S., & Zavisca, J. (2014). BIC and alternative Bayesian information criteria in the selection of structural equation models. Structural equation modeling: a multidisciplinary journal, 21(1), 1-19.
#'
#' @export

HBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n/(2*pi)))
}
