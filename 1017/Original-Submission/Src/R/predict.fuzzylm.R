#' Predict Method for Fuzzy Linear Regression
#'
#' Predicts the central tendency and spreads from a fuzzy linear regression model. 
#' @param object a \code{fuzzylm} object.
#' @param newdata an optional data frame in which to look for variables with which to 
#'   predict. If omitted, the fitted values are used.
#' @param ... further arguments passed to or from other methods.
#' @return \code{fuzzylm} object with \code{newdata} replacing the slot \code{x} and 
#'   predictions in triangular fuzzy number format representing the central values
#'   and left and right spreads replacing the slot \code{y}. 
#' @keywords fuzzy
#' @export
#' @examples
#' data(fuzzydat)
#' f <- fuzzylm(y ~ x, data = fuzzydat$lee)
#' predict(f)


predict.fuzzylm <- function(object, newdata, ...){
	if(!inherits(object, "fuzzylm")){
		stop(gettextf("'%s' is not a fuzzylm object.", deparse(substitute(object))))
	}
	if (missing(newdata) || is.null(newdata)) {
		newdata <- object$x[, -1, drop = FALSE]
		colnames(newdata) <- colnames(object$x)[-1]
	}
	lims <- object$lims
	for(i in 2:nrow(lims)){
		v <- rownames(lims)[i]
		if(any(newdata[, v] < lims[i, 1] | newdata[, v] > lims[i, 2])){
			stop(gettextf("'%s' is out of range. Prediction from the fuzzy linear model is defined only for data ranges used to fit the model.", v))
		}
	}
	xvars = colnames(object$x)
	xvars = xvars[2:length(xvars)]
	newdata <- newdata[, match(xvars, colnames(newdata)), drop = FALSE]
	n <- length(xvars)
	ct <- NULL
	fy <- matrix(rep(object$coef[1,], nrow(newdata)), ncol = 3, byrow = TRUE, 
				 dimnames = list(rownames(newdata), colnames(object$coef)))
	if(object$method %in% c("PLRLS", "FLS", "OPLR", "PLR", "FLAR")) {
		for(i in 2:nrow(object$coef)){
			for(j in 1:nrow(newdata)){
				cti <- prodSfuzzy(y = object$coef[i,], x = newdata[j, i-1])
				ct <- c(ct, cti)
			}
		fy <- sumFuzzy(fy, matrix(ct, ncol = 3, byrow = TRUE))
		}
			
	} 
	if(object$method == "MOFLR"){
		n <- n/2
		# sort fuzzy input numbers
		fx <- paste(1:n, (n+1):(2*n), (n+1):(2*n), sep=",")
		for(i in 2:nrow(object$coef)){
			for(j in 1:nrow(newdata)){
				cti <- prodFuzzy(object$coef[i,], newdata[j, as.numeric(unlist(strsplit(fx[i-1], ",")))])
				ct <- c(ct, cti)
			}
		fy <- sumFuzzy(fy, matrix(ct, ncol = 3, byrow = TRUE))
		}
	}
	newdata <- cbind(rep(1, nrow(newdata)), newdata)
	colnames(newdata)[1] <- "(Intercept)"
	res <- list(call = object$call, method = object$method, fuzzynum = object$fuzzynum, 
				coef = object$coef, lims = object$lims, x = newdata, y = fy)
	class(res) <- "fuzzylm"
	res	
}