#' validationMatern function
#'
#' Check the validity of a multivariate Matérn covariance model based on
#' the sufficient conditions given by Apanasovich et al. (2012)
#'
#' REFERENCE:
#' Tatiyana V. Apanasovich, Marc G. Genton, and Ying Sun. 2012.
#' A valid Matérn class of cross-covariance functions for multivariate random fields
#' with any number of components. J. Amer. Statist. Assoc. 107 (497), 180-193.
#'
#' @param a Matrix of scale parameters
#' @param nu Matrix of shape parameters
#' @param C Matrix of colocated correlation
#' @param d Space dimension
#' @return V: 1 if sufficient conditions are fulfilled, 0 otherwise
#' @examples
#' # Checking the validity of a bivariate Matérn model with given parameter matrices
#' a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
#' nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
#' C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' # dimension
#' d <- 2
#' # Checking the validity of parameters of bivariate Matérn model
#' validationMatern(a, nu, C, d)
#'
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
validationMatern <- function(a, nu, C, d){
  # Check symmetric matrices
  if ((isSymmetric.matrix(a)==FALSE) || (isSymmetric.matrix(nu)==FALSE) || 
      (isSymmetric.matrix(C)==FALSE)){
    stop("At least one of the matrices is not symmetric")
  }
  # Default parameters
  deltamax <- 1e4  # maximum tested value of "Delta" parameters of Apanasovich et al. (2012)
  maxit <- 10000   # number of iterations for bisection method
  tol <- 1e-10     # tolerance for numerical imprecision
  
  # First check on "nu" matrix (condition (i) of Apanasovich et al.)
  #-----------------------------------------------------------------
  a = 1/a
  p <- NROW(a)
  V <- 0
  averagenu <- 0.5*(t(matrix(rep(diag(nu), each = p),ncol=p)) + matrix(rep(diag(nu), each = p),ncol=p))
  deltanu <- nu - averagenu
  Aij = 1-deltanu/deltamax
  eig = eigen(Aij)
  lambda = min(eig$values)
  if(min((deltanu)< -tol) || (lambda < -tol)){
    print('Matrix of shape parameters does not fulfill sufficient condition (i) of Apanasovich et al. (2012)')
    V <- V + 1
  }
  if (max(c(deltanu)) < tol){
    deltaA <- 0
  } 
  else{
    delta1 = 0
    delta2 = deltamax
    for (nit in 1:maxit){
      delta = (delta1+delta2)/2
      Aij = 1-deltanu/delta
      eig = eigen(Aij)
      lambda = min(Re(eig$values))
      if(lambda < -tol){
        delta1 = delta
      }
      else{
        delta2 = delta
      }
    }
    deltaA = delta2
  }
  
  # Second check on "a" matrix (condition (ii) of Apanasovich et al.)
  #------------------------------------------------------------------
  aa <- a*a
  Bij <- - aa + (matrix(1, p, p)*aa[,p]) + (t(matrix(1, p, p)*aa[,p])) - (matrix(1, p, p)*aa[p,p])
  Bij <- Bij[1:(p-1),1:(p-1)]
  eig <- eigen(Bij)
  lambda <- min(eig$values)
  if (lambda < (-tol)){
    print('Matrix of scale parameters does not fulfill sufficient condition (ii) of Apanasovich et al. (2012)')
    V <- V + 1
  }
  
  # Third check on "C" matrix (condition (iii) of Apanasovich et al.)
  #--------------------------------------------------------------------
  R = C*(a^(2*(deltaA+averagenu)))*exp(lgamma(nu+d/2)-lgamma(averagenu+d/2)-lgamma(nu))
  eig = eigen(R)
  lambda = min(eig$values)
  if(lambda < (-tol)){
    print('Matrix of colocated correlation does not fulfill sufficient condition (iii) of Apanasovich et al. (2012)')
    V <- V + 1
  } 
  
  # Conclusion
  # ----------
  if (V == 0){
    print(paste(c("Congratulations. Your multivariate Matern model is valid in R",d), collapse = ""))
    V <- 1
  }
  else{
    V <- 0
  }
  return(V)
}
