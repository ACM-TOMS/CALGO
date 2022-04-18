\name{cova}
\alias{cova}
\title{Covariance values for two models: Matérn (ordinary covariance) and Power (generalized covariance)
}
\usage{
cova(model, h, b)
}
\description{
This function allows computing covariance values for two models: Matérn and Power
}
\arguments{
  \item{model}{One of the covariance models: \code{"Matern"} or \code{"Power"}
}
  \item{h}{
lag separation distance
}
  \item{b}{
parameter of the model: for Matérn, it is the shape parameter; for power, it is the exponent
}
}
\value{
The covariance value for the model type with the specified lag separation distance and the shape parameter or exponent
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# Compute covariance for a distance h and Matérn covariance
h <- 10
b <- 0.5
model <- "Matern"
cova(model, h, b)
}
