\name{variograms}
\alias{variograms}
\title{Comparison of the experimental direct and cross-variograms of \eqn{N} realizations (calculated along one grid axis) with the theoretical direct and cross-variograms}
\usage{
variograms(simu, nnodes, nvar, N, dimsimu, hmax, model, parameters, variogname)
}
\arguments{
\item{simu}{File containing a matrix with realizations}

\item{nnodes}{A vector with the number of grid nodes along X and Y directions (if \eqn{d=2}), or X, Y and Z directions (if \eqn{d=3})}

\item{nvar}{Number of random field components}

\item{N}{Number of realizations}

\item{dimsimu}{Along which axis (between \eqn{1} and \eqn{d}) the variograms will be calculated}

\item{hmax}{The maximum distance (grid meshes) for the calculation of variograms in the direction of the chosen axis}

\item{model}{ Theoretical variogram model}

\item{parameters}{List of parameter values}

\item{variogname}{Name of PNG output file}
}
\value{
A PNG file(s) with direct and cross-variograms (experimental, average and theoretical)
}
\description{This function allows plotting experimental direct and cross-variograms for \eqn{N} realizations (green lines), average of experimental direct and cross-variograms (red lines) and theoretical direct and cross-variograms (black lines).
Valid for simulations over a regular grid in the Euclidean space \eqn{R^d} (with \eqn{d=2} or \eqn{d=3})
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# Direct and cross-variograms of a bivariate Matérn random field on a regular grid 
# with 500x500 nodes 
# Thirty realizations (N = 30) are generated
set.seed(9784498)
nx <- 500
ny <- 500
coord <- cbind(c(rep(1, times = ny))\%x\%c(1:nx), (c(1:ny)\%x\%c(rep(1, times = nx))))
model <- "Matern"
# calculate the scale parameters, shape parameters and colocated correlation matrices
a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
simu <- specsim(coord, a0 = 10, nu0 = 0.5, model, parameters, N = 30, L = 1000, "Matern")
# plot of direct and cross-variograms calculated along the abscissa axis for distances
# ranging from 0 to 120 units
variograms(simu, c(nx,ny), nvar=2, N=30, dimsimu=1, hmax=120, model, parameters,  
"Matern")
}
