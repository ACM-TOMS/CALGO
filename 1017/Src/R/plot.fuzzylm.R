#' Plot Fuzzy Linear Regression
#'
#' Plots the data and the central tendency with spreads of a fuzzy linear regression.   
#' For multiple regression, allows choice of which variable to display. Optionally colors 
#' the polygon for the regression.
#' @param x a \code{fuzzylm} object.
#' @param y NULL for plotting a \code{fuzzylm} object.
#' @param which an integer or character string specifying which explanatory variable to 
#'   plot.
#' @param res an integer \code{>= 2} specifying resolution of shading for the regression  
#'    plot. Minimum resolution for shading the plot is 3.
#' @param col.fuzzy color for shading of the regression plot.
#' @param length length of the edges of the arrow head (in inches).
#' @param angle angle from the shaft of the arrow to the edge of the arrow head.
#' @param main a main title for the plot. Default title specifies method used to fit 
#'    the model.
#' @param xlab a label for the x axis, defaults to a description of x.
#' @param ylab a label for the y axis, defaults to a description of y.
#' @param ... additional graphical parameters.
#' @details Silently plots the data. Fuzzy numbers are plotted with points for the central 
#'   value and arrows specifying spreads.
#' @keywords fuzzy
#' @export
#' @examples
#' data(fuzzydat)
#' f = fuzzylm(y ~ x, fuzzydat$lee)
#' plot(f)
#' plot(f, res = 20, col.fuzzy = "red")
plot.fuzzylm = function(x, y = NULL, which = 1, res = 2, col.fuzzy = NA, length = 0.05, angle = 90, main = "method", xlab = NULL, ylab = NULL, ...){
	# assumes intercept in first column
	xc <- ifelse(is.numeric(which), which + 1, which(colnames(x$x) == which))
	coefs <- x$coef
	X <- x$x[, xc]
	y <- as.matrix(x$y)[, 1]
	cf <- function(x) coefs[1,1] + coefs[2,1] * x
	lf <- function(x) coefs[1,1] - coefs[1,2] + (coefs[2,1] - coefs[2,2]) * x
	rf <- function(x) coefs[1,1] + coefs[1,3] + (coefs[2,1] + coefs[2,3]) * x
	x0 <- min(X)
	x1 <- max(X)
	y0c <- cf(x = x0)
	y0l <- lf(x = x0)
	y0r <- rf(x = x0)
	y1c <- cf(x = x1)
	y1l <- lf(x = x1)
	y1r <- rf(x = x1)

	if(res < 2) res <- 2
	y0 <- c(seq(y0l, y0c, length.out = res), seq(y0c, y0r, length.out = res))
	y1 <- c(seq(y1l, y1c, length.out = res), seq(y1c, y1r, length.out = res))
	ylims <- range(y0, y1)
	if(x$method != "PLRLS") ylims <- range(ylims, y - x$y[, 2], y + x$y[, ifelse(ncol(as.matrix(x$y)) == 2, 2, 3)])
	xlims <- range(x0, x1)
	if(x$method == "MOFLR"){
		xc2 <- xc + (ncol(x$x)-1) / 2
		xlims <- range(xlims, X - x$x[, xc2], X + x$x[, xc2])
	}
	if(main == "method") main <- paste("Fuzzy linear regression using the", toupper(x$method), "method")
	
	if(is.null(xlab)) { xlab <- colnames(x$x)[xc] }
	if(is.null(ylab)) { ylab <- all.vars(x$call)[1] }
	
	graphics::plot(1, type = "n", ylim = ylims, xlim = xlims, xlab = xlab, 
		 ylab = ylab, main = main, ...)
	s <- grDevices::dev.size("in")
	
	if(res > 2){
		if(is.na(col.fuzzy)) col.fuzzy = "grey"
		cols = grDevices::colorRampPalette(c("white", col.fuzzy, "white"))(res * 2)
		for(i in 2:(res * 2)){
			graphics::polygon(x = c(x0, x1, x1, x0), y = c(y0[i-1], y1[i-1], y1[i], y0[i]), col = cols[i], border = NA)
		}
	}
	
	graphics::points(X, y, ...)

	if(ncol(as.matrix(x$y)) == 2){
		spread <- x$y[, 2]
		if(any(spread < s[2] / 1000)) spread <- spread + (s[2] / 1000)
		graphics::arrows(X, y, X, y - spread, length = length, angle = angle, ...)
		graphics::arrows(X, y, X, y + spread, length = length, angle = angle, ...)
	}
	if(ncol(as.matrix(x$y)) == 3){
		spread <- x$y[, 2:3]
		if(any(spread < s[2] / 1000)) spread <- spread + (s[2] / 1000)
		graphics::arrows(X, y, X, y - spread[, 1], length = length, angle = angle, ...)
		graphics::arrows(X, y, X, y + spread[, 2], length = length, angle = angle, ...)
	}
	
	if(x$method == "MOFLR"){
		spread <- x$x[, xc2]
		if(any(spread < s[1] / 1000)) { spread <- spread + (s[1] / 1000) }
		graphics::arrows(X, y, X - spread, y, length = length, angle = angle, ...)
		graphics::arrows(X, y, X + spread, y, length = length, angle = angle, ...)
	}
	
	graphics::abline(a = coefs[1,1], b = coefs[2,1], ...)
	graphics::abline(a = coefs[1,1] - coefs[1,2], b = coefs[2,1] - coefs[2,2], lty = 3, ...)
	graphics::abline(a = coefs[1,1] + coefs[1,3], b = coefs[2,1] + coefs[2,3], lty = 3, ...)
}


