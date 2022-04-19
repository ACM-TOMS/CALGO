#' Sum of Two Triangular Fuzzy Numbers
#'
#' Calculates a sum of two triangular fuzzy numbers defined as a 
#' central value, left and right spread.
#' @inheritParams prodFuzzy
#' @return Returns a numeric vector, representing a triangular fuzzy number.
#' @export
#' @examples 
#'   x <- c(1, 0.1, 0.2)
#'   y <- c(2, 0.2, 0.2)
#'   sumFuzzy(x = x, y = y)

sumFuzzy <- function(x, y){
	if(!inherits(x, "matrix")) x <- matrix(x, nrow = 1)
	if(!inherits(y, "matrix")) y <- matrix(y, nrow = 1)
	res <- matrix(c(x[,1] + y[,1],
					x[,2] + y[,2],
					x[,3] + y[,3]), ncol = 3, byrow = FALSE)
	colnames(res) <- c("central value", "left spread", "right spread")
	res
}
