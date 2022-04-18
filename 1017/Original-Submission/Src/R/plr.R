#' Fuzzy Linear Regression Using the Possibilistic Linear Regression Method
#'
#' The function calculates fuzzy regression coeficients using the possibilistic linear
#' regression method (PLR) developed by Tanaka et al. (1989). Specifically, the 
#' \code{min} problem is implemented in this function.
#' @param x matrix of \emph{n} independent variable observations. The first column is
#'    related to the intercept, so it consists of ones. Missing values not allowed.
#' @param y two column matrix of dependent variable values and the respective spread. 
#'    Method assumes symmetric triangular fuzzy input, so the second spread (if present)  
#'    is ignored. Missing values not allowed.
#' @param h a scalar value in interval \code{[0,1]}, specifying the h-level.
#' @details The function input expects the response in form of a symmetric fuzzy
#'    number and the predictors as crisp numbers. The prediction returns 
#'    symmetric triangular fuzzy number coefficients. 
#'    The h-level is a degree of fitting chosen by the decision maker.
#' @note Preferred use is through the \code{\link{fuzzylm}} wrapper function with argument
#'    \code{method = "plr"}.
#' @inherit fuzzylm return
#' @inherit plrls seealso
#' @references Tanaka H., Hayashi I. and Watada J. (1989) Possibilistic linear 
#'    regression analysis for fuzzy data. \emph{European Journal of Operational 
#'     Research} 40: 389-396.
#' @keywords fuzzy
#' @export
#' @import limSolve
#' @examples
#' data(fuzzydat)
#' fuzzylm(y ~ x, fuzzydat$tan, "plr", , , "yl", "yr")


plr = function(x, y, h = 0){
	if(h < 0 | h > 1) stop("h must be a value between 0 and 1.")
	n <- nrow(x)
	k <- ncol(x)
	vars <- colnames(x)
	X <- x
	f <- c(rep(0, k), colSums(abs(X)))
	A1 <- -1 * cbind(X, (1 - h) * abs(X))
	b1 <- -1 * matrix(y[, 1] + (1 - h) * y[, 2], ncol = 1)
	A2 <- cbind(X, -(1 - h) * abs(X))
	b2 <- matrix(y[, 1] - (1 - h) * y[, 2], ncol = 1)
	A3 <- cbind(matrix(0, k, k), -1 * diag(k))
	b3 <- matrix(0, nrow = k, ncol = 1)
	A <- rbind(A1, A2, A3)
	b <- rbind(b1, b2, b3)
	p <- limSolve::linp(E = NULL, F = NULL, G = -A, H = -b, Cost = f, ispos = FALSE)
	res <- matrix(p$X, ncol = 2, byrow = FALSE, dimnames = list(vars, c("center", "left.spread")))
	res <- cbind(res, res[, 2])
	colnames(res)[3] <- "right.spread"
	lims <- t(apply(x, 2, range))
	rownames(lims) <- vars
	colnames(lims) <- c("min", "max")
	fuzzy <- list(call = NULL, method = "PLR", fuzzynum = "symmetric triangular", coef = res, 
		 		  lims = lims, x = x, y = y)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}

