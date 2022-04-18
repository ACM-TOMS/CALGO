#' Fuzzy Linear Regression using the Possibilistic Linear Regression with Least Squares Method
#'
#' The function calculates fuzzy regression coeficients using the possibilistic linear
#' regression with least squares method developed by Lee and Tanaka (1999)
#' that combines the least squares approach (fitting of a central tendency) with the
#' possibilistic approach (fitting of spreads) when approximating an observed linear
#' dependence by a fuzzy linear model.
#' @param x two column matrix with the second column representing independent variable
#'    observations. The first column is related to the intercept, so it consists of ones.
#'    Missing values not allowed.
#' @param y one column matrix of dependent variable values, missing values not allowed.
#' @param h a scalar value in interval \code{[0,1]}, specifying the h-level.
#' @param k1 weight coefficient for the centeral tendency.
#' @param k2 weight coefficient for the spreads.
#' @param epsilon small positive number that supports search for the optimal solution.
#' @details The function input expects crisp numbers of both the explanatory and response
#'    variables, and the prediction returns non-symmetric triangular fuzzy number 
#'    coefficients.
#'    
#'    The h-level is a degree of fitting chosen by the decision maker.
#' @note Preferred use is through the \code{\link{fuzzylm}} wrapper function with argument
#'    \code{method = "plrls"}.
#' @inherit fuzzylm return
#' @references Lee, H. and Tanaka, H. (1999) Fuzzy approximations with non-symmetric fuzzy
#'    parameters in fuzzy regression analysis. \emph{Journal of the Operations Research
#'    Society Japan} 42: 98-112.
#' @keywords fuzzy
#' @seealso \code{\link{fuzzylm}}
#' @export
#' @import quadprog
#' @examples
#' x <- matrix(c(rep(1, 15), rep(1:3, each = 5)), ncol = 2)
#' y <- matrix(c(rnorm(5, 1), rnorm(5, 2), rnorm(5, 3)), ncol = 1)
#' plrls(x = x, y = y)

plrls = function(x, y, h = 0, k1 = 1, k2 = 1, epsilon = 1e-5){
	if(h < 0 | h > 1) stop("h must be a value between 0 and 1.")
	n <- ncol(x)
	m <- nrow(x)
	X <- x
	y <- as.matrix(y)
	vars <- colnames(x)
	j <- matrix(1, ncol = 1, nrow = m)
	X.abs <- abs(X)
	f1 <- -2 * k1 * t(X) %*% y
	f2 <- k2 * (1-h) * t(X.abs) %*% j
	f <- rbind(f1, f2, f2)

	H0 <- matrix(0, nrow = n, ncol = n)
	H1 <- t(X) %*% X
	H2 <- diag(n) * epsilon
	H <- cbind(k1 * H1, H0, H0)
	H <- rbind(H, cbind(H0, H2, H0))
	H <- rbind(H, cbind(H0, H0, H2))
	H <- 2 * H
	X.pos <- X * (X >= 0)
	X.neg <- X * (X < 0)
	L1 <- cbind(X, -(1-h) * X.neg, (1-h) * X.pos)
	L2 <- cbind(X, -(1-h) * X.pos, (1-h) * X.neg)
	L3 <- cbind(H0, H2, H0)
	L4 <- cbind(H0, H0, H2)
	L <- rbind(-L1, L2, -L3, -L4)
	r <- rbind(-y, y, matrix(0, ncol = 1, nrow = 2*n))
	p <- quadprog::solve.QP(Dmat = H, dvec = -f, Amat = -t(L), bvec = -r)

	coefs <- matrix(p$solution, nrow = n, byrow = FALSE,
		dimnames = list(vars, c("center", "left.spread", "right.spread")))
	lims <- t(apply(X, 2, range))
	rownames(lims) <- vars
	colnames(lims) <- c("min", "max")
	fuzzy <- list(call = NULL, x = x, y = y, lims = lims,
		method = "PLRLS", fuzzynum = "non-symmetric triangular", coef = coefs)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}

