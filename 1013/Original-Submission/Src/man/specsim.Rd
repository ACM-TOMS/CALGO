\name{specsim}
\alias{specsim}
\title{Continuous spectral simulation of random fields}
\usage{
specsim(coord, a0, nu0, model, parameters, N, L, filename)
}
\arguments{
\item{coord}{Coordinates of the locations targeted for simulation}

\item{a0}{Scale parameter of the Matérn spectral density \eqn{g} used to simulate the frequency vectors}

\item{nu0}{Shape parameter of the Matérn spectral density \eqn{g}}

\item{model}{selected covariance model, chosen from \code{"spherical"},
\code{"cubic"}, \code{"penta"}, \code{"exponential"}, \code{"Gaussian"}, \code{"Askey"},
\code{"Wendland"}, \code{"Matern"}, \code{"JBessel"}, \code{"Laguerre"}, \code{"Hypergeometric"}, \code{"Power"}, \code{"Spline"}, \code{"iterated_exponential"}, \code{"Gamma"}, \code{"Cauchy"}, \code{"stable1"}, \code{"stable2"}, \code{"mixed_Power"}. }

\item{parameters}{List of parameter values} 

\item{N}{Number of realizations to generate}

\item{L}{Number of basic random fields to use}

\item{filename}{Name of output txt file}
}

\value{
An \eqn{n * (N * P)} matrix containing the \eqn{N} simulated values for each coordinate (total = \eqn{n}) and for each random field component (\eqn{P})  
}
\description{
This function allows simulating scalar and vector random fields defined in an Euclidian space by a continuous spectral algorithm
}
\examples{
# Simulation of a bivariate Matérn random field on a regular grid with 500x500 nodes
set.seed(9784498)
nx <- 500
ny <- 500
coord <- cbind(c(rep(1, times = ny))\%x\%c(1:nx), (c(1:ny)\%x\%c(rep(1, times = nx))))
model <- "Matern"
# calculating the scale parameters, shape parameters and colocated correlation matrices
a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
# dimension
d <- 2
# Checking the validity of the bivariate Matérn model parameters
validationMatern(a, nu, C, d)
# Building a list of parameters for simulation
parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
# Simulation and generation of a txt file with name "Matern"
simu <- specsim(coord, a0 = 10, nu0 = 0.5, model, parameters, N = 1, L = 1000, "Matern")
# Display of realization of simulated bivariate random field
maps(simu, viewsimu="XY", coord, c(nx,ny), c(1,1), slice=1, "Matern")
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
