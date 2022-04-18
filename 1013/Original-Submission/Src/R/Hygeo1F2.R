#' Hygeo1F2 function
#' 
#' Generalized Hypergeometric function: particular cases for Askey and Wendland models
#'    Askey: with \eqn{d=2} and \eqn{\nu=2}
#'    Wendland: with \eqn{d=2}, \eqn{\nu=3} and \eqn{k=1}
#'    
#' @param z point where the hypergeometric function will be evaluated
#' @param model Covariance model: Askey or Wendland
#' @return Value of generalized hypergeometric function evaluated at \eqn{z} in \eqn{R}
#' @examples
#' # Evaluation of Hypergeometric function at z=1, for Askey covariance
#' z <- 1
#' model <- "Askey"
#' Hygeo1F2(z,model)
#' 
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
Hygeo1F2 <- function(z, model) {
  # error checking - model name is correct
  if (model %in% c("Askey", "Wendland") == FALSE) {
    stop("Model name invalid", call. = FALSE)
  }
  y <- 2*sqrt(z)
  struveH0 <- besselJ(y,1)+(7+(20/pi))*((1-cos(y))/y)+(36/pi-12)*((sin(y)-y*cos(y))/(y^2))
  struveH1 <- (2/pi)-besselJ(y,0)+(16/pi-5)*(sin(y)/y)+(12-(36/pi))*(1-cos(y))/(y^2)
  if (model == "Askey") { 
    # Calculate hypergeometric function 1F2(3/2,5/2,3,-z)
    H <- (3*besselJ(y,1)*(-2+pi*sqrt(z)*struveH0))/z^(3/2) - (3*besselJ(y,0)*(-2+pi*struveH1))/z
  }
  if (model == "Wendland") { 
    # Calculate hypergeometric function 1F2(5/2,4,9/2,-z)
    H <- (105*(24*z*besselJ(y,0)-3*pi*(4*z-5)*struveH1-2*y*besselJ(y,1)+3*pi*(4*z-5)*besselJ(y,1)*struveH0))/(16*z^3)
  }
  return(H)
}
