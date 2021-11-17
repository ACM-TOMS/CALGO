\name{test}
\alias{test}
\title{Testing the distribution of the simulated random fields with Shapiro-Wilk normality test}
\usage{
test(simu, coord1, coord2, lambda, model, parameters, N, filename)
}
\arguments{
\item{simu}{File with realizations}

\item{coord1}{First set of locations}

\item{coord2}{Second set of locations}

\item{lambda}{Set of weights}

\item{model}{Covariance, variogram or generalized covariance model of the target random field
}
  \item{parameters}{
List of parameter values
}
  \item{N}{
Number of realizations
}
  \item{filename}{
Name of PNG output file
}
}
\value{A probability-probability (P-P) plot: proportions of times when the test is accepted (value 0) and when the test is rejected (value 1) are counted and compared to the significance levels
}
\description{
This function allows testing the distribution of the generalized increments of simulated random fields with Shapiro-Wilk normality test, for selected bivariate models
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# A particular case for a simulated bivariate Matérn random field with L=100 and L=1000
L <- c(100, 1000)
N <- 10000
coord <- matrix(c(1, 1, 11, 1, 1, 11, 11, 11, 6, 1, 1, 6, 11, 6, 6, 11, 6, 6), nrow=9, ncol=2,
byrow=TRUE)
model <- "Matern"
# calculate the scale parameters, shape parameters and colocated correlation matrices
a <- matrix(c(5, (0.5*(5^-2+10^-2))^-0.5, (0.5*(5^-2+10^-2))^-0.5, 10), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
# dimension
d <- 2
# Checking the validity of the parameters of bivariate Matérn model
validationMatern(a, nu, C, d)
# Building a list of parameters to the simulation
parameters.SRF2 <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
# Simulation with L=100 and L=1000, respectively:
simu.SRF2.L100 <- specsim(coord, a0 = 15, nu0 = 0.5, model, parameters.SRF2, N, L[1])
simu.SRF2.L1000 <- specsim(coord, a0 = 15, nu0 = 0.5, model, parameters.SRF2, N, L[2])
# Put the two simulations together in an array
simu.SRF2 <- array(c(simu.SRF2.L100, simu.SRF2.L1000), dim = c(9, 2*N, 2))
# Preparing to the test
coord1 <- coord[c(1,2,3,4,9), ]
coord2 <- coord[5:9, ]
lambda <- c(0.25, 0.25, 0.25, 0.25, -1)
# Test and generate a png file named "test_SRF2"
test(simu.SRF2, coord1, coord2, lambda, model, parameters.SRF2, N, "test_SRF2")
}
