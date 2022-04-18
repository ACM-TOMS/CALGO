\name{validationMatern}
\alias{validationMatern}
\title{
Check the validity of a multivariate Matérn covariance model
}
\usage{
validationMatern(a, nu, C, d)
}
\arguments{

\item{a}{Matrix of scale factors}

\item{nu}{Matrix of shape factors}

\item{C}{Matrix of colocated correlation}

\item{d}{Space dimension}
}
\value{V: 1 if sufficient conditions are fulfilled, 0 otherwise
}
\description{
Check the validity of a multivariate Matérn covariance model based on the sufficient conditions given by Apanasovich et al. (2012)
}
\references{
Apanasovich T. V., Genton M. G., Sun Y (2012), A valid Matérn class of cross-covariance functions for multivariate random fields with any number of components. \emph{Journal of the American Statistical Association} 107 (497), 180-193.
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# Checking the validity of a bivariate Matérn model with given parameter matrices
a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
# dimension
d <- 2
# Checking the validity of parameters of bivariate Matérn model
validationMatern(a, nu, C, d)
}
