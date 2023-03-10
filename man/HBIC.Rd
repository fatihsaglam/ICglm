% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HBIC.R
\name{HBIC}
\alias{HBIC}
\title{Haughton Bayesian information criterion}
\usage{
HBIC(model)
}
\arguments{
\item{model}{a "lm" or "glm" object}
}
\value{
HBIC measurement of the model
}
\description{
Calculates Haughton Bayesian information criterion (HBIC) for "lm" and "glm" objects.
}
\details{
HBIC (Bollen et al., 2014) is calculated as

\deqn{-2LL(theta) + klog(n/(2pi))}
}
\examples{
x1 <- rnorm(100, 3, 2)
x2 <- rnorm(100, 5, 3)
x3 <- rnorm(100, 67, 5)
err <- rnorm(100, 0, 4)

## round so we can use it for Poisson regression
y <- round(3 + 2*x1 - 5*x2 + 8*x3 + err)

m1 <- lm(y~x1 + x2 + x3)
m2 <- glm(y~x1 + x2 + x3, family = "gaussian")
m3 <- glm(y~x1 + x2 + x3, family = "poisson")

HBIC(m1)
HBIC(m2)
HBIC(m3)

}
\references{
Bollen, K. A., Harden, J. J., Ray, S., & Zavisca, J. (2014). BIC and alternative Bayesian information criteria in the selection of structural equation models. Structural equation modeling: a multidisciplinary journal, 21(1), 1-19.
}
