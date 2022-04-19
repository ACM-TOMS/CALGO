#' Total Error of Fit of Fuzzy Regression Model
#'
#' Calculates total error of fit of a fuzzy regression model based on the
#' concept of difference in membership functions of triangular fuzzy numbers
#' between the estimated and observed fuzzy dependent variables.
#' @param object a \code{fuzzylm} object.
#' @param sc scaling constant used for numerical stability when spreads are equal to zero.
#' @param ... additional arguments passed to the \code{integrate} function.
#' @details Calculates \eqn{\sum{E}}{sum(E)}, where \eqn{E}{E} is the difference in 
#'   membership functions between two triangular fuzzy numbers. Here, between the 
#'   observation and the prediction from a fuzzy regression model \code{fuzzylm}.
#' @note \code{TEF} is not suitable for assessing fuzzy linear regression models that were
#'   fitted from crisp input data. Such data will result in division by zero. The scaling
#'   constant \code{sc} numerically allows the calculation to proceed, but it is not 
#'   advisable. Use \code{\link{GOF}} instead.
#' @return A numeric with sum of pairwise differences between the triangular fuzzy
#'   numbers.
#' @seealso \code{\link{fuzzylm}}, \code{\link{GOF}}
#' @references Kim B. and Bishu R. R. (1998) Evaluation of fuzzy linear regression
#'   models by comparing membership functions. \emph{Fuzzy Sets and Systems}
#'   100: 343-352.
#' @export
#' @examples
#' data(fuzzydat)
#' f <- fuzzylm(y ~ x, fuzzydat$lee)
#' TEF(f)

TEF = function(object, sc = 1e-6, ...){
  if(!inherits(object, "fuzzylm")){
	stop(gettextf("'%s' is not a fuzzylm object.", deparse(substitute(object))))
  }

  TFN1 = as.matrix(object$y)
  TFN1 = switch(dim(TFN1)[2], cbind(TFN1, rep(sc, nrow(TFN1)), rep(sc, nrow(TFN1))),
                              cbind(TFN1, TFN1[,2]),
                              TFN1)
  TFN2 = predict.fuzzylm(object)$y
  if(nrow(TFN1) != nrow(TFN2)) stop("Compare the same number of triangular fuzzy numbers")
  if(ncol(TFN1) != 3 | ncol(TFN2) != 3)
  	if(nrow(TFN1) != 3 | nrow(TFN2) != 3){
  		stop("Both triangular fuzzy numbers must be matrices with three columns")
  	} else {
  		TFN1 = t(TFN1)
  		TFN2 = t(TFN2)
  	}
  e = 0
    for(i in 1:nrow(TFN1)){
    lims = c(TFN1[i,1]-TFN1[i,2], TFN1[i,1]+TFN1[i,3], TFN2[i,1]-TFN2[i,2], TFN2[i,1]+TFN2[i,3])
    D = stats::integrate(function(x) abs(dom(x, TFN1[i,]) - dom(x, TFN2[i,])), min(lims), max(lims),
                           stop.on.error = FALSE)$value
    int = D / stats::integrate(function(x) dom(x, TFN1[i,]), min(lims[1:2]), max(lims[1:2]),
                           stop.on.error = FALSE)$value
    e = e + int
    }
  e
}
