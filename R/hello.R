n <- 500

x1 <- rnorm(n = n, mean = 34, sd = 2)
x2 <- rnorm(n = n, mean = 36, sd = 3)
err <- rnorm(n = n, mean = 0, sd = 1)

y <- 2 + -x1 + 3*x2 + err

m_lm <- lm(y~x1*x2)

m_glm <- glm(round(y)~x1*x2, family = "poisson")

reverse_fisher <- function(model){
  cov_coef <- vcov(model)
  n <- length(model$residuals)
  lower_left <- 2*var(model$residuals)^2/n
  as.matrix(Matrix::bdiag(cov_coef, lower_left))
}


AIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 2*df)
}

BIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n))
}

BIC_Q <- function(model, q = 0.5) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n) - 2*df*log(q/(1-q)))
}

BIC_adj <- function(model, q = 0.5) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log((n+2)/24))
}

CAIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*(log(n) + 1))
}

KIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 3*df)
}

AIC_4 <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 4*df)
}

HQIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  c(-2*LL + 2*df*log(log(n)))
}

FIC <- function(model) {
  LL <- logLik(object = model)
  X <- model.matrix(model)
  c(-2*LL + log(det(crossprod(X))))
}

ICOMP_IFIM_CF <- function(model) {
  LL <- logLik(object = model)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  CF <- mean((ev - m)^2)
  c(-2*LL + 2*CF)
}

ICOMP_IFIM_C1 <- function(model) {
  LL <- logLik(object = model)
  ev <- eigen(reverse_fisher(model))$values
  p <- length(ev)
  gm <- prod(ev)^(1/p)
  m <- mean(ev)
  C1 <- p/2*log(m/gm)
  c(-2*LL + 2*C1)
}

ICOMP_IFIM_C1F <- function(model) {
  LL <- logLik(object = model)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  C1F <- sum((ev - m)^2)/(4*m^2)
  c(-2*LL + 2*C1F)
}

ICOMP_IFIM_C1R <- function(model) {
  LL <- logLik(object = model)
  R <- cov2cor(reverse_fisher(model))
  C1R <- -0.5*log(det(R))
  c(-2*LL + 2*C1R)
}




ICOMP_PEU_CF <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  CF <- mean((ev - m)^2)
  c(-2*LL + 2*CF + df)
}

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

ICOMP_PEU_C1F <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  C1F <- sum((ev - m)^2)/(4*m^2)
  c(-2*LL + 2*C1F + df)
}

ICOMP_PEU_C1R <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  R <- cov2cor(reverse_fisher(model))
  C1R <- -0.5*log(det(R))
  c(-2*LL + 2*C1R + df)
}

ICOMP_PEU_LN_CF <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  CF <- mean((ev - m)^2)
  c(-2*LL + 2*log(n)*CF + df)
}

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

ICOMP_PEU_LN_C1F <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  ev <- eigen(reverse_fisher(model))$values
  m <- mean(ev)
  C1F <- sum((ev - m)^2)/(4*m^2)
  c(-2*LL + 2*log(n)*C1F + df)
}

ICOMP_PEU_LN_C1R <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  R <- cov2cor(reverse_fisher(model))
  C1R <- -0.5*log(det(R))
  c(-2*LL + 2*log(n)*C1R + df)
}

CAIC_F <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*(log(n) + 2) + log(det(solve(reverse_fisher(model)))))
}

JIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + 0.5*(df*log(n) - n*log(1 - df/n)))
}

IBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n/(2*pi) + log(det(solve(reverse_fisher(model))))))
}

SPBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  Fisher_M <- solve(reverse_fisher(model))
  coeff <- c(model$coefficients, var(model$residuals))
  E <- coeff%*%Fisher_M%*%coeff
  c(-2*LL + df*(1 - log(df/(E))))
}

HBIC <- function(model) {
  LL <- logLik(object = model)
  df <- attr(LL, "df")
  n <- length(model$residuals)
  c(-2*LL + df*log(n/(2*pi)))
}

GCV <- function(model) {
  LL <- logLik(object = model)
  res <- model$residuals
  RSS <- mean(res^2)
  df <- attr(LL, "df")
  RSS/(1 - df/n)^2
}


AIC(m_lm)
BIC(m_lm)
BIC_Q(m_lm)
BIC_adj(m_lm)
CAIC(m_lm)
KIC(m_lm)
AIC_4(m_lm)
HQIC(m_lm)
FIC(m_lm)
ICOMP_IFIM_CF(m_lm)
ICOMP_IFIM_C1(m_lm)
ICOMP_IFIM_C1F(m_lm)
ICOMP_IFIM_C1R(m_lm)
ICOMP_PEU_CF(m_lm)
ICOMP_PEU_C1(m_lm)
ICOMP_PEU_C1F(m_lm)
ICOMP_PEU_C1R(m_lm)
ICOMP_PEU_LN_CF(m_lm)
ICOMP_PEU_LN_C1(m_lm)
ICOMP_PEU_LN_C1F(m_lm)
ICOMP_PEU_LN_C1R(m_lm)
CAIC_F(m_lm)
JIC(m_lm)
IBIC(m_lm)
SPBIC(m_lm)
GCV(m_lm)


IC <- function(model, criterias = NULL, ...) {
  if (is.null(criterias)) {
    criterias <- c("AIC",
                   "BIC",
                   "BIC_Q",
                   "BIC_adj",
                   "CAIC",
                   "KIC",
                   "AIC_4",
                   "HQIC",
                   "FIC",
                   "ICOMP_IFIM",
                   "ICOMP_PEU",
                   "ICOMP_PEU_LN",
                   "CAIC_F")
  }

  ff <- lapply(criterias, function(m) match.fun(m))

}
