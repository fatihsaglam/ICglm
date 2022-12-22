#' @title  Bayesian Information Criterion
#'
#' @description Calculates Bayesian Information Criterion (BIC) and its variants (BICadj, BICQ) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#' @param q adjustment parameter for \code{BICQ}. Default is 0.25.
#'
#' @details
#' BIC (Schwarz, 1978) is calculated as
#'
#' \deqn{-2LL(theta) + klog(n)}
#'
#' Adjusted BIC (Dziak et al., 2020) as
#'
#' \deqn{-2LL(theta) + klog(n/2pi)}
#'
#' and BICQ (Xu, 2010) as
#'
#' \deqn{-2LL(theta) + klog(n) - 2klog(q/(1-q))}.
#'
#' @return BIC, BICadj or BICQ measurement of the model
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
#'BIC(m1)
#'BIC(m2)
#'BIC(m3)
#'BICadj(m1)
#'BICadj(m2)
#'BICadj(m3)
#'
#' @references
#'
#' Dziak, J. J., Coffman, D. L., Lanza, S. T., Li, R., & Jermiin, L. S. (2020). Sensitivity and specificity of information criteria. Briefings in bioinformatics, 21(2), 553-565.
#'
#' Xu, C. (2010). Model Selection with Information Criteria.
#'
#' Schwarz, G. 1978. Estimating the dimension of a model The Annals of Statistics 6 (2), 461–464. <doi:10.1214/aos/1176344136>
#'
#' @rdname BIC
#' @export

BIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n))
}

#' @rdname BIC
#' @export
BICadj <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log((n+2)/24))
}

#' @rdname BIC
#' @export
BICQ <- function(model, q = 0.25) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n) - 2*df*log(q/(1-q)))
}
