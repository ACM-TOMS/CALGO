#' Fuzzy Linear Regression using the Fuzzy Least Squares Method
#'
#' The function calculates fuzzy regression coeficients using the fuzzy least squares 
#' (FLS) method proposed by Diamond (1988) for non-symmetric triangular fuzzy numbers.
#' @param x two column matrix with the second column representing independent variable
#'    observations. The first column is related to the intercept, so it consists of ones.
#'    Missing values not allowed.
#' @param y matrix of dependent variable observations. The first column contains the 
#'    central tendency, the second column the left spread and the third column the right
#'    spread of non-symmetric triangular fuzzy numbers. Missing values not allowed.
#' @details The FLS method for the fuzzy linear regression fits a simple model.
#' @note Preferred use is through the \code{\link{fuzzylm}} wrapper function with argument
#'    \code{method = "fls"}.
#' @inherit fuzzylm return
#' @inherit plrls seealso
#' @references Diamond, P. (1988) Fuzzy least squares. \emph{Information Sciences}
#'    46(3): 141-157.
#' @keywords fuzzy
#' @export
#' @examples
#'    data(fuzzydat)
#'    x <- fuzzydat$dia[, 1, drop = FALSE]
#'    x <- cbind(rep(1, nrow(x)), x)
#'    y <- fuzzydat$dia[, c(2,3,3)]
#'    fls(x = x, y = y)


fls <- function(x, y){
	if(ncol(x) != 2) stop("Only one independent variable allowed for the FLS method")
    if(ncol(y) != 3) stop("Only one dependent variable with two spreads allowed for the FLS method")
	vars <- colnames(x)
	n <- nrow(x)
	sx <- sum(x[,2])
	sx2 <- sum(x[,2]^2)
	sy <- sum(y[,1])
	syel <- sum(y[,2])
	syer <- sum(y[,3])
	sxy <- sum(x[,2] * y[,1])
	sxyel <- sum(x[,2] * y[,2])
	sxyer <- sum(x[,2] * y[,3])
	S <- matrix(c(3*n, 3*sx, -n, -sx, n, sx,
				3*sx, 3*sx2, -sx, -sx2, sx, sx2,
				-n, -sx, n, sx, 0, 0,
				-sx, -sx2, sx, sx2, 0, 0,
				n, sx, 0, 0, n, sx,
				sx, sx2, 0, 0, sx, sx2), byrow = TRUE, ncol = 6)
	P <- matrix(c(3*sy-syel+syer,
				3*sxy-sxyel+ sxyer,
				-sy+syel,
				-sxy+sxyel,
				sy+syer,
				sxy+sxyer), ncol = 1)
	C <- solve(S) %*% P
	coefs <- matrix(C, ncol = 3, byrow = FALSE,
		dimnames = list(vars, c("center", "left.spread", "right.spread")))
	lims <- t(apply(x, 2, range))
	rownames(lims) <- vars
	colnames(lims) <- c("min", "max")
	fuzzy <- list(call = NULL, x = x, y = y, lims = lims,
		method = "fls", fuzzynum = "non-symmetric triangular", coef = coefs)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}


