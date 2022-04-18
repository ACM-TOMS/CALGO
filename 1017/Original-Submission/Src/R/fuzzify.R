#' Convert Real Value Numbers to Triangular Fuzzy Numbers
#' 
#' Uses naive alternative methods to approximate triangular fuzzy
#' numbers from real value number input data.
#' @param x numeric vector.
#' @param y vector that can be coerced to factor (optional).
#' @param method character vector specifying the conversion method. See Details.
#' @param err numeric vector. Error term for the error method.
#' @param dimnames \code{list} of length 2 giving names of the \code{x} and \code{y} 
#'   variables.
#' @param ... additional parameters passed to other functions.
#' @details Converts crisp numbers in \code{x} to a triangular fuzzy number (TFN). Optionally, 
#'  values in \code{y} can be used as grouping elements and are coerced to a factor.
#'  
#'  Method \code{mean} calculates the central value of a TFN as the mean of \code{x} given
#'  \code{y}, and the left and right spreads as standard deviations. 
#'  
#'  Method \code{median} gives the central values as a median and left and right spreads 
#'  are calculated as distance of the first and third quartile from the median.
#'  
#'  Method \code{zero} inserts zeros to both spreads.
#'  
#'  Method \code{error} uses a user-defined numeric value or vector for the spreads. 
#'  The length of the numeric vector in argument \code{err} must be in 
#'  (1, \code{length(x)}, 2 * \code{length(x)}).
#' @export
#' @return A data.frame with columns representing the central value, left and right spread of
#'  \code{x} and the values in \code{y} coerrced to a factor. Attempt is made to inherit 
#'  names from the input data. Methods 
#'  \code{mean} and \code{zero} will return symmetric TFNs, whereas methods \code{median} and 
#'  \code{error} can return non-symmetric TFNs depending on input data and the data or the
#'  values in the \code{err} argument.
#' @examples 
#' fuzzify(1:5)
#' fuzzify(1:6, c(1,1,1,2,2,2), method = "err", err = runif(6) * 1e-3)


fuzzify = function(x, y = NULL, method = "mean", err = 0, dimnames = list("x", "y"), ...){
  # check input
  metody = c("mean", "median", "zero", "error")
  m = pmatch(method, metody)
  if(is.na(m)) stop(sprintf("%s was not uniquely matched to available methods %s.", method, 
                            paste(metody, collapse = ", ")))
  method = metody[m]
  
  if(is.null(y)){
    y = rep(1, length(x))
  }
  
  if(length(x) != length(y)) stop("x and y must have the same length")
  
  # conversions
  if(method == "mean"){
    res = stats::aggregate(x, by = list(y), FUN = mean, ...)
    res[, ncol(res)+1] = stats::aggregate(x, by = list(y), FUN = stats::sd, ...)$x
    res = res[, c(2,3,3,1)]
  }
  
  if(method == "median"){
    res = stats::aggregate(x, by = list(y), FUN = stats::median, ...)
    res[, ncol(res) + 1] = res[, 2] - stats::aggregate(x, by = list(y), FUN = stats::quantile, probs = 0.25, ...)$x
    res[, ncol(res) + 1] = stats::aggregate(x, by = list(y), FUN = stats::quantile, probs = 0.75, ...)$x - res[, 2]
    res = res[, c(2,3,4,1)]
  }
  
  if(method == "zero"){
    res = data.frame(x = x, xl = rep(0, length(x)), xr = rep(0, length(x)), y = y)
  }
  
  if(method == "error"){
    if(length(err) == 1){
      err = rep(err, 2 * length(x))
    }
    if(length(err) == length(x)){
      err = c(err, err)
    }
    if(all(length(err) != length(x), length(err) != 2 * length(x))){
      stop("spreads included in the err argument length not a multiple of length of x")
    }
    res = data.frame(matrix(c(x, err), ncol = 3, byrow = FALSE))
    res$y = y
  }
  
  # find names
  if(any(sapply(dimnames, FUN = function(x) length(x) != 1))){
    warning("Using first variable name from dimnames")
  }
  name.x = ifelse(is.null(names(x)), dimnames[[1]][1], names(x))
  name.y = ifelse(is.null(names(y)), dimnames[[2]][1], names(y))
  colnames(res) = c(paste0(name.x, "c"),
                    paste0(name.x, "l"),
                    paste0(name.x, "r"),
                    name.y)
  res
}
