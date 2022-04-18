#' Real Value Degree of Membership to a Triangular Fuzzy Number
#'
#' Calculates the degree of membership of a real number to a 
#' triangular fuzzy number. The fuzzy number is defined by its 
#' central value and the left and right spreads.
#' @param x a numeric vector.
#' @param TFN a numeric vector of length 3.
#' @return Returns a numeric in interval \code{[0,1]}.
#' @export
#' @examples
#' x <- seq(from = 0, to = 2, length.out = 10)
#' A <- c(1, 1, 1)
#' dom(x, A)

dom <- function(x, TFN){
  res <- rep(0, length(x))
  TFN[2] <- ifelse(TFN[2] == 0, 1e-6, TFN[2])
  TFN[3] <- ifelse(TFN[3] == 0, 1e-6, TFN[3])
  mn <- TFN[1] - TFN[2]
  ce <- TFN[1]
  mx <- TFN[1] + TFN[3]
  mask1 <- mn < x & x < ce
  mask2 <- ce <= x & x < mx
  res <- res + ((x - mn) / TFN[2]) * mask1
  res <- res + ((mx - x) / TFN[3]) * mask2
  res
}
