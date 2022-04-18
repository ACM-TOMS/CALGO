#' Fuzzy Linear Regression Using the Possibilistic Linear Regression with Omission Method
#'
#' The function calculates fuzzy regression coeficients using the possibilistic linear 
#' regression with an outlier omission approach method (OPLR) developed by Hung and Yang 
#' (2006) that combines the least squares approach (fitting of a central tendency) with 
#' the possibilistic approach (fitting of spreads) when approximating an observed linear
#' dependence by a fuzzy linear model.
#' @param x matrix with the independent variables observations. The first column is
#'    related to the intercept, so it consists of ones. Missing values not allowed.
#' @param y two column matrix of the dependent variable values and the respective spread. 
#'    Method assumes symmetric triangular fuzzy input, so the second spread (if present) is 
#'    ignored. Missing values not allowed.
#' @param h a scalar value in interval \code{[0,1]}, specifying the h-level.
#' @details The function input expects symmetric fuzzy response and crisp predictors. The 
#'    prediction returns symmetric triangular fuzzy number coefficients.
#'    The OPLR method can detect one outlier in the data that is farther than 
#'    \code{1.5 * IQR} from either quartile.
#'    
#'    The h-level is a degree of fitting chosen by the decision maker.
#' @note Preferred use is through the \code{\link{fuzzylm}} wrapper function with argument
#'    \code{method = "oplr"}.
#' @inherit fuzzylm return
#' @inherit plrls seealso
#' @references Hung, W.-L. and Yang, M.-S. (2006) An omission approach for detecting 
#'    outliers in fuzzy regression models. \emph{Fuzzy Sets and Systems} 157: 3109-3122.
#' @keywords fuzzy
#' @export
#' @import limSolve
#' @examples
#' data(fuzzydat)
#' fuzzylm(y ~ x, fuzzydat$hun, "oplr", , , "yl")


oplr <- function(x, y, h = 0){
	n <- nrow(x)
	p <- ncol(x)
	vars <- colnames(x)
	X <- x
	JMi <- matrix(NA, nrow = n, ncol = 1)
	Ai <- matrix(NA, nrow = n, ncol = 2*p)
	
	for(i in 1:n){
		Xi <- X[-i, ]
		yi <- y[-i, ]
		f <- c(rep(0, p), colSums(abs(Xi)))
		A1 <- -1 * cbind(Xi, (1 - h) * abs(Xi))
		b1 <- -1 * matrix(yi[, 1] + (1 - h) * yi[, 2], ncol = 1)
		A2 <- cbind(Xi, -(1 - h) * abs(Xi))
		b2 <- matrix(yi[, 1] - (1 - h) * yi[, 2], ncol = 1)
		A3 <- cbind(matrix(0, p, p), -1 * diag(p))
		b3 <- matrix(0, nrow = p, ncol = 1)
		A <- rbind(A1, A2, A3)
		b <- rbind(b1, b2, b3)
		lpi <- limSolve::linp(E = NULL, F = NULL, G = -A, H = -b, Cost = f, ispos = FALSE)
		JMi[i, ] <- lpi$solutionNorm
		Ai[i, ] <- lpi$X
	}
	
	f <- c(rep(0, p), colSums(abs(X)))
	A1 <- -1 * cbind(X, (1 - h) * abs(X))
	b1 <- -1 * matrix(y[, 1] + (1 - h) * y[, 2], ncol = 1)
	A2 <- cbind(X, -(1 - h) * abs(X))
	b2 <- matrix(y[, 1] - (1 - h) * y[, 2], ncol = 1)
	A3 <- cbind(matrix(0, p, p), -1 * diag(p))
	b3 <- matrix(0, nrow = p, ncol = 1)
	A <- rbind(A1, A2, A3)
	b <- rbind(b1, b2, b3)
	
	lp <- limSolve::linp(E = NULL, F = NULL, G = -A, H = -b, Cost = f, ispos = FALSE)
	lp.obj <- lp$solutionNorm
	lp.coef <- lp$X
	
	ri <- abs(JMi - lp.obj)/lp.obj
	Q <- stats::quantile(ri, probs = c(0.25, 0.75))
	iqr <- stats::IQR(ri)
	ot1 <- ri > (Q[2] + 1.5 * iqr)
	ot2 <- ri < (Q[1] - 1.5 * iqr)
	
	if(sum(ot1) > 1 | sum(ot2) > 1 | sum(ot1, ot2) > 1){ stop("Multiple outliers detected. Revise the data or use method = \"moflr\"") }
	if(max(ot1) == 1){ 
		res <- Ai[ot1, ] 
		warning("Deleting one outlier (ot1)")
	} 
	if(max(ot2) == 1){ 
		res <- Ai[ot2, ] 
		warning("Deleting one outlier (ot2)")
	}
	if(all(!c(ot1, ot2))) { 
		res <- lp.coef 
		warning("No outliers detected")
	}
	
	res <- matrix(res, ncol = 2, byrow = FALSE, dimnames = list(vars, c("center", "left.spread")))
	res <- cbind(res, res[, 2])
	colnames(res)[3] <- "right.spread"
	lims <- t(apply(x, 2, range))
	rownames(lims) <- vars
	colnames(lims) <- c("min", "max")
	fuzzy <- list(call = NULL, method = "OPLR", fuzzynum = "symmetric triangular", coef = res, 
		 		  lims = lims, x = x, y = y)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}


