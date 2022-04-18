#' Goodness of Fit of Fuzzy Regression Model
#'
#' Calculates mean error rate based on Diamond's distance of two variables representing
#' triangular fuzzy numbers, where one is the response variable and the other is the
#' prediction from a fuzzy regression model. 
#' @param object a \code{fuzzylm} object.
#' @param sc scaling constant used for numerical stability when spreads are equal to zero.
#' @details The Diamond's distance of two triangular fuzzy numbers is sum of squared differences
#'  of the core and both support values of the fuzzy numbers.
#' @return A numeric.
#' @seealso \code{\link{fuzzylm}}
#' @references Diamond P. (1988) Fuzzy least squares. \emph{Information Sciences}
#'  46(3): 141-157.
#'    
#'  Wang N., Zhang W.-X., Mei C.-L. (2007) Fuzzy nonparametric regression
#'  based on local linear smoothing technique. \emph{Information Sciences}
#'  177: 3882-3900.
#' @export
#' @examples
#' data(bats)
#' f <- fuzzylm(temperature ~ MAST, data = bats)
#' GOF(f)


GOF = function(object, sc = 1e-6){
  if(!inherits(object, "fuzzylm")){
    stop(gettextf("'%s' is not a fuzzylm object.", deparse(substitute(object))))
  }
  
  TFN1 = as.matrix(object$y)
  TFN1 = switch(dim(TFN1)[2], cbind(TFN1, rep(sc, nrow(TFN1)), rep(sc, nrow(TFN1))),
                cbind(TFN1, TFN1[,2]),
                TFN1)
  TFN2 = predict.fuzzylm(object)$y
  if(nrow(TFN1) != nrow(TFN2)) stop("Compare the same number of triangular fuzzy numbers")
  if(ncol(TFN1) != 3 | ncol(TFN2) != 3){
    if(nrow(TFN1) != 3 | nrow(TFN2) != 3){
      stop("Both triangular fuzzy numbers must be matrices with three columns")
    } else {
      TFN1 = t(TFN1)
      TFN2 = t(TFN2)
    }
  }
  n = nrow(TFN1)
  lsq = ((TFN1[,1] - TFN1[,2]) - (TFN2[,1] - TFN2[,2])) ^ 2
  csq = (TFN1[,1] - TFN2[,1]) ^ 2
  rsq = ((TFN1[,1] + TFN1[,3]) - (TFN2[,1] + TFN2[,3])) ^ 2
  gof = sum(lsq, csq, rsq) / n
  gof
}
