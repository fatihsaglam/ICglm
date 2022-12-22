#' @title  Hannan-Quinn Information Criterion
#'
#' @description Calculates Hannan-Quinn Information Criterion (HQIC) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' HQIC (Hannan and Quinn, 1979) is calculated as
#'
#' \deqn{-2LL(theta) + 2klog(log(n))}
#'
#' @return HQIC measurement of the model
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
#' HQIC(m1)
#' HQIC(m2)
#' HQIC(m3)
#'
#' @references
#' Hannan, E. J., & Quinn, B. G. (1979). The determination of the order of an autoregression. Journal of the Royal Statistical Society: Series B (Methodological), 41(2), 190-195.
#'
#' @export

HQIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + 2*df*log(log(n)))
}
