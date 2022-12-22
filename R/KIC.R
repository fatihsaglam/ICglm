#' @title  Kullback–Leibler Information Criterion
#'
#' @description Calculates Kullback–Leibler Information Criterion (KIC) and its corrected form (KICC) for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' KIC (Seghouane, 2006) is calculated as
#'
#' \deqn{-2LL(theta) + 3k}
#'
#' and KICC (Seghouane, 2006) is calculated as
#'
#' \deqn{-2LL(theta) + ((k + 1)(3n - k - 2)) + (k/(n-k))}
#'
#' @return KIC measurement of the model
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
#'KIC(m1)
#'KIC(m2)
#'KIC(m3)
#'KICC(m1)
#'KICC(m2)
#'KICC(m3)
#'
#' @references
#' Seghouane, A. K. (2006). A note on overfitting properties of KIC and KICC. Signal Processing, 86(10), 3055-3060.
#'
#' @rdname KIC
#' @export

KIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 3*df)
}

#' @rdname KIC
#' @export
KICC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + ((df + 1)*(3*n - df - 2))/(n - df - 2) + df/(n - df))
}
