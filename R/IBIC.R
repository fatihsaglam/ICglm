#' @title Information Matrix-Based Information Criterion
#'
#' @description Calculates Information Matrix-Based Information Criterion (IBIC) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' IBIC (Bollen et al., 2012) is calculated as
#'
#' \deqn{-2LL(theta) + klog(n/(2pi)) + log(|F|)}
#'
#' \eqn{F} is the fisher information matrix.
#'
#' While calculating the Fisher information matrix (\eqn{F}), we used
#' the joint parameters (\eqn{beta,sigma^2}) of the models.
#'
#' @return IBIC measurement of the model
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
#'IBIC(m1)
#'IBIC(m2)
#'IBIC(m3)
#'
#' @references
#' Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A comparison of Bayes factor approximation methods including two new methods. Sociological Methods & Research, 41(2), 294-324.
#'
#' @export

IBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n/(2*pi) + log(det(solve(reverse_fisher(model))))))
}
