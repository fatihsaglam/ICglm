#' @title Joint Information Criterion
#'
#' @description Joint Information Criterion (JIC) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' JIC (Rahman and King, 1999) is calculated as
#'
#' \deqn{-2LL(theta) + 1/2*(klog(n) - nlog(1-k/n))}
#'
#' @return JIC measurement of the model
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
#'JIC(m1)
#'JIC(m2)
#'JIC(m3)
#'
#' @references
#'Rahman, M. S., & King, M. L. (1999). Improved model selection criterion. Communications in Statistics-Simulation and Computation, 28(1), 51-71.
#'
#' @export
#'

JIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + 0.5*(df*log(n) - n*log(1 - df/n)))
}
