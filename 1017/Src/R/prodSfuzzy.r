#' Product of a Scalar and a Triangular Fuzzy Number
#'
#' Calculates product of a real number scalar and a triangular fuzzy number defined as a 
#' central value, left and right spread.
#' @param x numeric vector of length one.
#' @param y a numeric vector of length three, specifying a triangular fuzzy number as its 
#'   central value, left and right spread.
#' @return Returns a numeric vector, representing a triangular fuzzy number.
#' @details Note that if \code{x < 0} the left and right spread will be reversed.
#' @export
#' @examples 
#'   x <- 2
#'   y <- c(2, 0.2, 0.2)
#'   prodSfuzzy(x = x, y = y)
#'   x <- -2
#'   prodSfuzzy(x = x, y = y)

prodSfuzzy <- function(x, y){
	if(length(x) != 1) stop("x must be a numeric of length one.")
	if(length(y) != 3) stop("y must be a numeric vector of length three, specifying central value, left and right spread of a fuzzy number.")
	if(x >= 0) res <- x * y
	if(x < 0) res <- c(x * y[1], abs(x) * y[3:2])
	names(res) <- c("central value", "left spread", "right spread")
	res
}