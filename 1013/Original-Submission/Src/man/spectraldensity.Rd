\name{spectraldensity}
\alias{spectraldensity}
\title{
Spectral density associated with a given covariance model
}
\usage{
spectraldensity(u, model, parameters)
}
\arguments{
  \item{u}{
frequency vector
}
  \item{model}{selected covariance model, chosen from one of \code{"spherical"},
\code{"cubic"}, \code{"penta"}, \code{"exponential"}, \code{"Gaussian"}, \code{"Askey"},
\code{"Wendland"}, \code{"Matern"}, \code{"JBessel"}, \code{"Laguerre"}, \code{"Hypergeometric"}, \code{"Power"}, \code{"Spline"}, \code{"iterated_exponential"}, \code{"Gamma"}, \code{"Cauchy"}, \code{"stable1"}, \code{"stable2"}, \code{"mixed_Power"}
}
  \item{parameters}{
list of parameter values for example:
new.params <- list("C" = [number or matrix], "a" = [number or matrix], "nu1" = [number or matrix], "nu2" = [number or matrix]) 
}
}
\value{
A matrix-valued spectral density value associated with a covariance model, for the chosen frequency vector \eqn{u} in \eqn{R^d}
}
\description{
This function calculates the spectral density associated with a covariance model
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# Bivariate Matern spectral density matrix for a vector u in R^2
u <- matrix(rnorm(1), 1, 2)
model <- "Matern"
# calculate the scale parameters, shape parameters and colocated correlation matrices
a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
spectraldensity(u, model, parameters)
}