#' @title  Information Criteria
#'
#' @description Various information criterias for "lm" and "glm" objects.
#'
#' @param model a "lm" or "glm" object
#'
#' @details
#' \code{AIC} is calculated using
#'
#' -2LL + 2k
#'
#' \code{BIC}
#' asdasdsad
#'
#' @return Information Criteria value
#'
#' @importFrom Matrix bdiag
#'
#'
#' @examples
#'
#' @references
#'
#' @rdname IC
#' @export

AIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 2*df)
}

#' @rdname IC
#' @export

BIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n))
}

#' @rdname IC
#' @export

BIC_Q <- function(model, q = 0.5) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n) - 2*df*log(q/(1-q)))
}

#' @rdname IC
#' @export

BIC_adj <- function(model, q = 0.5) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log((n+2)/24))
}

#' @rdname IC
#' @export

CAIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*(log(n) + 1))
}

#' @rdname IC
#' @export

KIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 3*df)
}
#' @rdname IC
#' @export

AIC_4 <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 4*df)
}
#' @rdname IC
#' @export

HQIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 2*df*log(log(n)))
}
#' @rdname IC
#' @export

FIC <- function(model) {
  LL <- logLik(object = model)
  X <- model.matrix(model)
  c(-2*LL + log(det(crossprod(X))))
}
#' @rdname IC
#' @export

ICOMP_IFIM_CF <- function(model) {
  LL <- logLik(object = model)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  CF <- mean((ev - m)^2)
  c(-2*LL + 2*CF)
}
#' @rdname IC
#' @export

ICOMP_IFIM_C1 <- function(model) {
  LL <- logLik(object = model)
  ev <- eigen(reverse_fisher(model))$values
  p <- length(ev)
  gm <- prod(ev)^(1/p)
  m <- mean(ev)
  C1 <- p/2*log(m/gm)
  c(-2*LL + 2*C1)
}
#' @rdname IC
#' @export

ICOMP_IFIM_C1F <- function(model) {
  LL <- logLik(object = model)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  C1F <- sum((ev - m)^2)/(4*m^2)
  c(-2*LL + 2*C1F)
}
#' @rdname IC
#' @export

ICOMP_IFIM_C1R <- function(model) {
  LL <- logLik(object = model)
  R <- cov2cor(reverse_fisher(model))
  C1R <- -0.5*log(det(R))
  c(-2*LL + 2*C1R)
}


#' @rdname IC
#' @export


ICOMP_PEU_CF <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  CF <- mean((ev - m)^2)
  c(-2*LL + 2*CF + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_C1 <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  ev <- eigen(reverse_fisher(model))$values
  p <- length(ev)
  gm <- prod(ev)^(1/p)
  m <- mean(ev)
  C1 <- p/2*log(m/gm)
  c(-2*LL + 2*C1 + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_C1F <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  C1F <- sum((ev - m)^2)/(4*m^2)
  c(-2*LL + 2*C1F + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_C1R <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  R <- cov2cor(reverse_fisher(model))
  C1R <- -0.5*log(det(R))
  c(-2*LL + 2*C1R + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_LN_CF <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  CF <- mean((ev - m)^2)
  c(-2*LL + 2*log(n)*CF + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_LN_C1 <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  ev <- eigen(reverse_fisher(model))$values
  p <- length(ev)
  gm <- prod(ev)^(1/p)
  m <- mean(ev)
  C1 <- p/2*log(m/gm)
  c(-2*LL + 2*log(n)*C1 + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_LN_C1F <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  C1F <- sum((ev - m)^2)/(4*m^2)
  c(-2*LL + 2*log(n)*C1F + df)
}
#' @rdname IC
#' @export

ICOMP_PEU_LN_C1R <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  R <- cov2cor(reverse_fisher(model))
  C1R <- -0.5*log(det(R))
  c(-2*LL + 2*log(n)*C1R + df)
}
#' @rdname IC
#' @export

CAIC_F <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*(log(n) + 2) + log(det(solve(reverse_fisher(model)))))
}
#' @rdname IC
#' @export

JIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + 0.5*(df*log(n) - n*log(1 - df/n)))
}
#' @rdname IC
#' @export

IBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n/(2*pi) + log(det(solve(reverse_fisher(model))))))
}
#' @rdname IC
#' @export

SPBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  Fisher_M <- solve(reverse_fisher(model))
  coeff <- c(model$coefficients, var(model$residuals))
  E <- coeff%*%Fisher_M%*%coeff
  c(-2*LL + df*(1 - log(df/(E))))
}
#' @rdname IC
#' @export

HBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n/(2*pi)))
}
#' @rdname IC
#' @export

GCV <- function(model) {
  LL <- logLik(object = model)
  res <- model$residuals
  RSS <- mean(res^2)
  df <- attr(LL, "df")
  RSS/(1 - df/n)^2
}
