\name{maps}
\alias{maps}
\title{Color scale representation of a realization}
\usage{
maps(simu, viewsimu, coord, nnodes, gridmesh, slice, mapname, width.png, height.png)
}
\arguments{
\item{simu}{File with one realization of the simulated random field}

\item{viewsimu}{View at the XY, XZ, or YZ plane: "XY", "XZ" or "YZ"}

\item{coord}{Coordinates of the locations where the random field has been simulated}

\item{nnodes}{A vector with the number of grid nodes along X and Y directions (\eqn{d=2}), or X, Y and Z directions (\eqn{d=3})}

\item{gridmesh}{A vector with the grid mesh along X and Y directions (if \eqn{d=2}: c(dx,dy)), or X, Y and Z directions (if \eqn{d=3}: c(dx,dy,dz))}

\item{slice}{Index of the grid slice to display}

\item{mapname}{Name of PNG output file}

\item{width.png}{Display width depending on whether the grid is square or rectangular. By default, the value is for a square grid}

\item{height.png}{Display height depending on whether the grid is square or rectangular. By default, the value is for a square grid}
}
\value{
A PNG file(s) with a map of one realization of the random field
}
\description{
This function allows displaying a color scale representation of a simulated random field.
Valid for simulation over a regular grid in the Euclidean space \eqn{R^d} with \eqn{d=2} or \eqn{d=3}
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# A map of one realization of a simulated bivariate Matérn random field on a regular grid 
# with 500x500 nodes
set.seed(9784498)
nx <- 500
ny <- 500
coord <- cbind(c(rep(1, times = ny))\%x\%c(1:nx), (c(1:ny)\%x\%c(rep(1, times = nx))))
model <- "Matern"
# Calculating the scale parameters, shape parameters and colocated correlation matrices
a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
# dimension
d <- 2
# Checking the validity of the parameters
validationMatern(a, nu, C, d)
# Building a list of parameters of the simulation
parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
# Generating one realization
simu <- specsim(coord, a0 = 10, nu0 = 0.5, model, parameters, N = 1, L = 1000, "Matern")
# Display of realization of simulated bivariate random field
maps(simu, viewsimu="XY", coord, c(nx,ny), c(1,1), slice=1, "Matern")
}
