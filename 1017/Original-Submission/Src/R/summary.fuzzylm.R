#' Summarizes Fuzzy Linear Regression 
#'
#' Calculates the summary from the \code{fuzzylm} object. 
#' @param object a \code{fuzzylm} object.
#' @param ... additional parameters passed to and from other methods.
#' @return Returns a list with models for the central tendency and spreads from the fuzzy
#'   linear regression, total error of fit of the model and a goodness-of-fit measure.
#' @seealso \code{\link{TEF}}, \code{\link{GOF}}
#' @keywords fuzzy
#' @export
#' @examples
#' data(fuzzydat)
#' f <- fuzzylm(y ~ x, fuzzydat$lee)
#' sum.f <- summary(f)
#' sum.f

summary.fuzzylm = function(object, ...){
	xvars = colnames(object$x)
	xvars = xvars[2:length(xvars)]
	yvars = colnames(object$y)
	if(is.null(yvars)) yvars = all.vars(object$call)[1]
	n = length(xvars)
	c = round(object$coef[,1], 4)
	l = round(object$coef[,2], 4)
	r = round(object$coef[,3], 4)
	tef = TEF(object, ...)
	gof = GOF(object)
	
	s = list(call = object$call, method = object$method, xvars = xvars, yvars = yvars, 
			 n = n, c = c, l = l, r = r, TEF = tef, GOF = gof)
	class(s) = "summary.fuzzylm"
	s
}

