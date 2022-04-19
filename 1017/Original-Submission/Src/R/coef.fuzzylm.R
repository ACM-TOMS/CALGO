#' Extract Model Coefficients from Fuzzy Linear Model
#'
#' Extracts coefficients of the fuzzy regression model in object \code{fuzzylm}.
#' @param object a \code{fuzzylm} object.
#' @param complete not used for a \code{fuzzylm} object.
#' @param ... other arguments.
#' @export
#' @return matrix with coefficients for the central tendency of the model, upper and lower
#'   boundary.
#' @note The function returns real value numbers that define model predictions at 
#'   \eqn{\mu_{\tilde{y}}(x) = 1}{mu_y(x) = 1} and \eqn{\mu_{\tilde{y}}(x) > 0}{mu_y(x) > 0},
#'   not triangular fuzzy numbers. To extract triangular fuzzy number coefficients of the 
#'   model, use \code{object$coef}.
#' @examples
#'   data(fuzzydat)
#'   f <- fuzzylm(y ~ x, data = fuzzydat$lee)
#'   coef(f)

coef.fuzzylm <- function(object, complete = TRUE, ...){
	xvars = colnames(object$x)
	yvars = colnames(object$y)
	if(is.null(yvars)) yvars = all.vars(object$call)[1]
	n = length(xvars)
	s = object$coef[,1]
	l = object$coef[,2]
	r = object$coef[,3]
	res = matrix(c(s, s - l, s + r), ncol = 3,
	             dimnames = list(xvars, c("central tendency", "lower boundary", "upper boundary")), ...)
	res
}


