#' Prints Fuzzy Linear Regression Result
#'
#' Prints the call and coefficients from the \code{fuzzylm} object. 
#' @param x a \code{fuzzylm} object.
#' @param ... further arguments passed to or from other methods.
#' @keywords fuzzy
#' @export
#' @examples
#' x <- rep(1:3, each = 5)
#' y <- c(rnorm(5, 1), rnorm(5, 2), rnorm(5, 3))
#' dat <- data.frame(x = x, y = y)
#' f <- fuzzylm(y ~ x, dat)
#' f


print.fuzzylm = function(x, ...){
	cat("\nFuzzy linear model using the", x$method, "method\n\nCall:\n")
	print(x$call)
	cat("\nCoefficients in form of", x$fuzzynum, "fuzzy numbers:\n\n")
	print(x$coef)
}
