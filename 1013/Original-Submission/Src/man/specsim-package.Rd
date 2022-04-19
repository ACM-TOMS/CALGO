\name{specsim-package}
\alias{specsim-package}
\docType{package}
\title{
Continuous spectral simulation of random fields
}
\details{
\tabular{ll}{
Package: \tab specsim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2020-09-01\cr
License: \tab GPL-2\cr
}
}
\description{
This package provides a set of functions that simulate vector Gaussian random fields in \eqn{R^d}. This package is not available from CRAN, but can easily be installed using the following R code: install.packages("specsim_1.0.tar.gz", repos=NULL, type = "source")
}
\section{Overview}{
The package provides eight functions:

1) \code{\link{cova}} and \code{\link{test}} are used for testing the distribution of the simulated random field for particular models and over a small set of points. 

2) \code{\link{spectraldensity}} and \code{\link{Hygeo1F2}} are used for calculating the spectral density for a given frequency vector \eqn{u} in \eqn{R^d} and specific covariance models. 

3) \code{\link{validationMatern}} checks sufficient validity conditions for the multivariate Matérn covariance model.

4) \code{\link{maps}} creates a color scale representation of realizations of a simulated vector random field.

5) \code{\link{variograms}} calculates the experimental variograms of realizations constructed on a regular grid.

6) Finally, \code{\link{specsim}} runs for a given number of realizations and basic random fields, specific parameters and models and coordinates of the target locations in \eqn{R^d}.
}
\author{
Daisy Arroyo and Xavier Emery

Maintainer: Daisy Arroyo <darroyof@udec.cl>

}
\references{
Apanasovich T. V., Genton M. G., Sun Y (2012), A valid Matern class of cross-covariance functions for multivariate random fields with any number of components. \emph{Journal of the American Statistical Association} 107 (497), 180-193.

Arroyo, D. Emery, X. (2017), Spectral simulation of vector random fields with stationary Gaussian increments in \eqn{d}-dimensional Euclidean spaces, \emph{Stochastic Environmental Research and Risk Assessment} 31 (7), 1583-1592.

Arroyo, D. Emery, X. (2018), Simulation of intrinsic random fields of order \eqn{k} with a continuous spectral algorithm, \emph{Stochastic Environmental Researchand Risk Assessment} 32 (11), 3245-3255.

Emery,  X.,  Arroyo,  D.  Porcu,  E.  (2016),  An  improved  spectral  turning-bands algorithm for simulating stationary vector Gaussian random fields, \emph{Stochastic Environmental Research and Risk Assessment} 30  (7), 1863-1873.
}