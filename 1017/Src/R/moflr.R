#' Fuzzy Linear Regression Using the Multi-Objective Fuzzy Linear Regression Method
#'
#' This function calculates fuzzy regression coeficients using the multi-objective fuzzy
#' linear regression (MOFLR) method developed by Nasrabadi et al. (2005)
#' that combines the least squares approach (fitting of a central tendency) with the
#' possibilistic approach (fitting of spreads) when approximating an observed linear
#' dependence by a fuzzy linear model.
#' @param x matrix of \emph{n} independent variable values, followed by \emph{n} spreads. 
#'    First column is exptected to consist of ones, representing intercept. Missing values 
#'    not allowed.
#' @param y two column matrix of dependent variable values and the respective spread. 
#'    Method assumes symmetric triangular fuzzy input, so the second spread (if present) is 
#'    ignored. Missing values not allowed.
#' @param omega a scalar that specifies weight that determines trade-off of 
#'    between outliers penalization and data fitting in interval \code{[0,1]}, where high 
#'    values of \code{omega} decrease the penalization of outliers.
#' @param sc scaling constant used to input random spreads for the intercept, necessary for
#'    computational stability.
#' @details The function input expects both the response and the predictors in form of
#'    symmetric fuzzy numbers. The 
#'    prediction returns symmetric triangular fuzzy number coefficients.
#'    The Nasrabadi et al.'s method can process datasets with multiple outliers. Values
#'    \code{omega>0.5} decrease weight of outliers on the solution.
#' @inherit fuzzylm return
#' @inherit plrls seealso
#' @note Preferred use is through the \code{\link{fuzzylm}} wrapper function with argument
#'    \code{method = "moflr"}.
#' @references Nasrabadi, M. M., Nasrabadi, E. and Nasrabady, A. R. (2005) Fuzzy linear 
#'    regression analysis: a multi-objective programming approach. \emph{Applied Mathematics
#'    and Computation} 163: 245-251.
#' @keywords fuzzy
#' @export
#' @import quadprog
#' @examples
#' data(fuzzydat)
#' fuzzylm(y~x, fuzzydat$nas, "moflr", "xl", , "yl")

moflr <- function(x, y, omega = 0.5, sc = 1e-5){		
	m <- nrow(x)
	n <- ncol(x)
	n <- (n - 1) / 2
	vars <- colnames(x[, 1:(n + 1)])		
	for(i in (n + 2):ncol(x)){
		if(sum(abs(x[, i])) == 0){
			x[, i] <- abs(stats::runif(m)) * sc
		}
	}
	Xv <- matrix(c(x[, 1:(n + 1)], abs(stats::runif(m)) * sc, x[, (n + 2):ncol(x)]), byrow = FALSE, ncol = 2 * n + 2)
	a <- solve(t(Xv[, 1:(n + 1)]) %*% Xv[, 1:(n + 1)]) %*% t(Xv[, 1:(n + 1)]) %*% y[, 1]
	T <- (a >= 0) - (a < 0)
	Xabs <- Xv[, (n + 2):ncol(Xv)] * matrix(rep(t(T), m), ncol = n + 1, byrow = TRUE)
	Xabs <- cbind(Xabs, abs(Xv[, 1:(n + 1)]))
	Ie <- diag(1, m)
	X2 <- cbind(Xv[, 1:(n + 1)], matrix(0, ncol = n + 1, nrow = m))
	Ze <- matrix(0, m, m)
	Ze2 <- matrix(0, m, 2 * (n + 1))
	XX <- cbind(omega * Xabs, Ze, Ze)
	XX <- rbind(XX, cbind(Ze2, (1-omega) * Ie, Ze))
	XX <- rbind(XX, cbind(Ze2, Ze, (1 - omega) * Ie))
	H <- 0.5 * t(XX) %*% XX
	yy <- rbind(y[,2, drop = FALSE], matrix(0, ncol = 1, nrow = 2*m))
	f <- -t(yy) %*% XX
	L <- cbind(X2, -Ie, Ze)
	L <- rbind(L, cbind(-X2, Ze, -Ie))
	L <- rbind(L, cbind(-Xabs, Ze, Ze))
	L <- rbind(L, cbind(Ze2, -Ie, Ze))
	L <- rbind(L, cbind(Ze2, Ze, -Ie))
	r <- matrix(c(y[, 1], -y[, 1], rep(0, 3 * m)), ncol = 1)
	p <- quadprog::solve.QP(Dmat = H, dvec = -f, Amat = -t(L), bvec = -r)
	res <- p$solution
	res[(n + 2):(2 * (n + 1))] <- round(res[(n + 2):(2 * (n + 1))], digits = nchar(format(sc, scientific = FALSE)) - 3)
	res <- matrix(res[1:(2 * (n + 1))], ncol = 2, byrow = FALSE, 
		          dimnames = list(vars, c("center", "left.spread")))
	res <- cbind(res, res[, 2])
	colnames(res)[3] <- "right.spread"
	lims <- t(apply(x[, 1:(n + 1)], 2, range))
	rownames(lims) <- vars
	colnames(lims) <- c("min", "max")
	fuzzy <- list(call = NULL, method = "MOFLR", fuzzynum = "symmetric triangular", coef = res, 
		 		  lims = lims, x = x, y = y)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}


