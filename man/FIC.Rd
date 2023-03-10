% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIC.R
\name{FIC}
\alias{FIC}
\title{Fisher Information Criterion}
\usage{
FIC(model)
}
\arguments{
\item{model}{a "lm" or "glm" object}
}
\value{
FIC measurement of the model
}
\description{
Calculates Fisher Information Criterion (FIC) for "lm" and "glm" objects.
}
\details{
FIC (Wei, 1992) is calculated as

\deqn{-2LL(theta) + log(|X^T X|)}
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

FIC(m1)
FIC(m2)
FIC(m3)

}
\references{
Wei, C. Z. (1992). On predictive least squares principles. The Annals of Statistics, 20(1), 1-42.
}
