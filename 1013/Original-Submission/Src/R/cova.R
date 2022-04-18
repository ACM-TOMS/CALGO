#' covariance values for two models: Matérn (ordinary covariance) and Power (generalized covariance)
#'
#' This function allows computing covariance values for two models: Matérn and Power 
#'
#' @param model model type:  
#' \code{"Matern"} Matérn model    
#' \code{"Power"} Power model of order k
#' @param h lag separation distance
#' @param b parameter of the model: for Matérn, it is the shape parameter; for power, it is the exponent
#' @return The covariance value for the model type with the specified lag separation distance and the shape parameter or exponent 
#' @examples
#' # Compute covariance for a distance h and Matérn covariance
#' h <- 10
#' b <- 0.5
#' model <- "Matern"
#' cova(model, h, b)
#' 
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#' 
cova <- function(model, h, b){
  # error checking - model name is correct
  if (model %in% c("Matern", "Power") == FALSE) {
    stop("Model name for test is invalid", call. = FALSE)
  }
  epsilon = 1e-12;
  if(model == "Matern"){               
    Cv = ((h+epsilon)/2)^b*2/gamma(b)*besselK((h+epsilon),b) 
  } 
  if(model == "Power"){           
    k = ceiling(b/2)
    Cv = ((-1)^k)*(h^b)
  }
  return(Cv)
}
