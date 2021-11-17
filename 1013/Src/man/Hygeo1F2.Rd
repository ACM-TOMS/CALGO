\name{Hygeo1F2}
\alias{Hygeo1F2}
\title{
Generalized Hypergeometric function: particular cases for Askey and Wendland models
}
\description{
Generalized Hypergeometric function: particular cases for Askey and Wendland models.
\tabular{ll}{
Askey: \tab with \eqn{d=2} and \eqn{\nu=2}\cr
Wendland: \tab with \eqn{d=2}, \eqn{\nu=3} and \eqn{k=1}\cr
}
}
\usage{
Hygeo1F2(z, model)
}
\arguments{
  \item{z}{
point where the hypergeometric function will be evaluated
}
  \item{model}{
covariance model
}
}
\value{A value of generalized hypergeometric function evaluated for \eqn{z} in \eqn{R}
}
\author{
Daisy Arroyo, \email{darroyof@udec.cl} and Xavier Emery, \email{xemery@ing.uchile.cl}
}
\examples{
# Evaluation of Hypergeometric function at z=1, for Askey covariance
z <- 1
model <- "Askey"
Hygeo1F2(z, model)
}