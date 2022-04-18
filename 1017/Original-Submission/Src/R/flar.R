#' Fuzzy Linear Regression using the Fuzzy Least Absolute Residual Method
#'
#' The function calculates fuzzy regression coeficients using the fuzzy least absolute
#' residual (FLAR) method proposed by Zeng et al. (2017) 
#' for non-symmetric triangular fuzzy numbers.
#' @param x matrix with the second to last columns representing independent variable
#'    observations. The first column is related to the intercept, so it consists of ones.
#'    Missing values not allowed.
#' @param y matrix of dependent variable observations. The first column contains the 
#'    central tendency, the second column the left spread and the third column the right
#'    spread of non-symmetric triangular fuzzy numbers. Missing values not allowed.
#' @details The FLAR method expects real value input for the explanatory variables, and 
#'    non-symmetric triangular fuzzy numbers for the response variable. The prediction 
#'    returns non-symmetric triangular fuzzy numbers.
#' @note Preferred use is through the \code{\link{fuzzylm}} wrapper function with argument
#'    \code{method = "flar"}.
#' @inherit fuzzylm return
#' @inherit plrls seealso
#' @references Zeng, W., Feng, Q. and Li, J. (2017) Fuzzy least absolute linear regression. 
#'    \emph{Applied Soft Computing} 52: 1009-1019.
#' @keywords fuzzy
#' @export
#' @examples
#'    data(fuzzydat)
#'    fuzzylm(y ~ x, fuzzydat$dia, "flar", , , "yl", "yl")

flar <- function(x, y){
	vars <- colnames(x)
	n <- nrow(x)
	p <- ncol(x)
	X <- x

	I <- diag(n)
	Ir <- diag(p)
	Z <- matrix(0, ncol = n, nrow = n)
	ZX <- matrix(0, nrow = n, ncol = p)
	ZXr <- matrix(0, nrow = p, ncol = p)
	Zr <- matrix(0, nrow = p, ncol = n)

	f <- c(rep(1, 6*n), rep(0, 3*p))

	Req <- cbind(I, -I, Z, Z, Z, Z, X, ZX, ZX)
	Req <- rbind(Req, cbind(Z, Z, I, -I, Z, Z, ZX, X, ZX))
	Req <- rbind(Req, cbind(Z, Z, Z, Z, I, -I, ZX, ZX, X))

	leq <- matrix(c(y))

	R <- cbind(-I, Z, Z, Z, Z, Z, ZX, ZX, ZX)
	R <- rbind(R, cbind(Z, -I, Z, Z, Z, Z, ZX, ZX, ZX))
	R <- rbind(R, cbind(Z, Z, -I, Z, Z, Z, ZX, ZX, ZX))
	R <- rbind(R, cbind(Z, Z, Z, -I, Z, Z, ZX, ZX, ZX))
	R <- rbind(R, cbind(Z, Z, Z, Z, -I, Z, ZX, ZX, ZX))
	R <- rbind(R, cbind(Z, Z, Z, Z, Z, -I, ZX, ZX, ZX))
	R <- rbind(R, cbind(Zr, Zr, Zr, Zr, Zr, Zr, ZXr, -Ir, ZXr))
	R <- rbind(R, cbind(Zr, Zr, Zr, Zr, Zr, Zr, ZXr, ZXr, -Ir))
	R <- rbind(R, cbind(Z, Z, Z, Z, Z, Z, ZX, -X, ZX))
	R <- rbind(R, cbind(Z, Z, Z, Z, Z, Z, ZX, ZX, -X))

	l <- matrix(rep(0, 8*n + 2*p))

	sorig <- limSolve::linp(E = Req, F = leq, G = -R, H = -l, Cost = f, ispos = FALSE)
	s <- sorig$X

	coefs <- matrix(c(s[(6*n+1):(6*n+p)],
			s[(6*n+p+1):(6*n+2*p)],
			s[(6*n+2*p+1):(6*n+3*p)]), ncol = 3, 
			dimnames = list(vars, c("center", "left.spread", "right.spread")))
	lims <- t(apply(x, 2, range))
	rownames(lims) <- vars
	colnames(lims) <- c("min", "max")
	fuzzy <- list(call = NULL, x = x, y = y, lims = lims,
		method = "fls", fuzzynum = "non-symmetric triangular", coef = coefs)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}