#' @title  Akaike Information Criterion
#'
#' @description Calculates Akaike Information Criterion (AIC) and its variants for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' AIC (Akaike, 1973) is calculated as
#'
#' \deqn{-2LL(theta) + 2k}
#'
#' and AIC4 (Bozdogan, 1994) as
#'
#' \deqn{-2LL(theta) + 2klog}
#'
#' @return AIC or AIC4 measurement of the model
#'
#' @importFrom stats logLik
#'
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
#' AIC(m1)
#' AIC(m2)
#' AIC(m3)
#' AIC4(m1)
#' AIC4(m2)
#' AIC4(m3)
#'
#' @references
#' Akaike H., 1973. Maximum likelihood identification of Gaussian
#' autoregressive moving average models. Biometrika, 60(2), 255-265.
#'
#' Bozdogan, H. 1994. Mixture-model cluster analysis using model
#' selection criteria and a new informational measure of complexity.
#' In Proceedings of the first US/Japan conference on the frontiers
#' of statistical modeling: An informational approach, 69–113.
#' Dordrecht: Springer.
#'
#' @rdname AIC
#' @export
AIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 2*df)
}

#' @rdname AIC
#' @export
AIC4 <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 4*df)
}
