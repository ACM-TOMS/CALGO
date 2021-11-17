#' Product of Two Triangular Fuzzy Numbers
#'
#' Calculates product of two triangular fuzzy numbers defined as a central value, left
#' and right spread.
#' @param x a numeric vector of length three, specifying a triangular fuzzy number as its 
#'   central value, left and right spread.
#' @param y a numeric vector of length three, specifying a triangular fuzzy number as its 
#'   central value, left and right spread.
#' @return Returns a numeric vector, representing a triangular fuzzy number.
#' @export
#' @examples x <- c(1, 0.2, 0.2)
#'   y <- c(2, 0.2, 0.2)
#'   prodFuzzy(x = x, y = y)

prodFuzzy <- function(x, y){
	if(length(x) != 3 | length(y) != 3){
		stop("Both numbers need to be vectors of length three, specifying central value, left and right spread of a fuzzy number.")
	}
	cv <- x[1] * y[1]
	mx <- matrix(rep(c(x[1] - x[2], x[1] + x[3]), 2), ncol = 2, byrow = TRUE) 
	my <- matrix(rep(c(y[1] - y[2], y[1] + y[3]), 2), ncol = 2, byrow = FALSE)
	l <- min(mx * my)
	r <- max(mx * my)
	res <- c(cv, cv - l, r - cv)
	names(res) <- c("central value", "left spread", "right spread")
	res
}