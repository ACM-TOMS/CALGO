#' Fuzzy Linear Regression
#'
#' A wrapper function that calculates fuzzy regression coeficients using a chosen method.
#' @param formula a model formula.
#' @param data a data.frame, containing the variables in formula.
#' @param method method for fitting of the fuzzy linear model. 
#' @param fuzzy.left.x character string vector specifying column name(s) with the left 
#'   spread of the fuzzy independent variable(s).
#' @param fuzzy.right.x character string vector specifying column name(s) with the right 
#'   spread of the fuzzy independent variable(s).
#' @param fuzzy.left.y character string vector specifying column name(s) with the left 
#'   spread of the fuzzy dependent variable.
#' @param fuzzy.right.y character string vector specifying column name(s) with the right 
#'   spread of the fuzzy dependent variable.
#' @param ... additional parameters used by specific methods.
#' @details The implemented methods include \code{\link{plrls}} for fitting the fuzzy linear
#'   regression from the crisp input data (Lee and Tanaka 1999), and \code{\link{fls}} 
#'   (Diamond 1988), \code{\link{oplr}} (Hung and Yang 2006), \code{\link{moflr}}
#'   (Nasrabadi et al. 2005) and \code{\link{plr}} (Tanaka et al. 1989) methods for
#'   triangular fuzzy numbers. 
#' @return Returns a \code{fuzzylm} object that includes the model coefficients, limits
#'   for data predictions from the model and the input data.
#' @seealso \code{\link[=plot.fuzzylm]{plot}}, \code{\link[=predict.fuzzylm]{predict}}, 
#'   \code{\link[=summary.fuzzylm]{summary}}
#' @references
#'   Diamond, P. (1988) Fuzzy least squares. \emph{Information Sciences}
#'   46(3): 141-157.
#'
#'   Hung, W.-L. and Yang, M.-S. (2006) An omission approach for detecting 
#'   outliers in fuzzy regression models. \emph{Fuzzy Sets and Systems} 157: 3109-3122.
#'
#'   Lee, H. and Tanaka, H. (1999) Fuzzy approximations with non-symmetric fuzzy
#'   parameters in fuzzy regression analysis. \emph{Journal of the Operations Research
#'   Society Japan} 42: 98-112.
#'
#'   Nasrabadi, M. M., Nasrabadi, E. and Nasrabady, A. R. (2005) Fuzzy linear 
#'   regression analysis: a multi-objective programming approach. \emph{Applied Mathematics
#'   and Computation} 163: 245-251.
#'
#'   Tanaka H., Hayashi I. and Watada J. (1989) Possibilistic linear 
#'   regression analysis for fuzzy data. \emph{European Journal of Operational 
#'   Research} 40: 389-396.
#'
#'   Zeng, W., Feng, Q. and Li, J. (2017) Fuzzy least absolute linear regression. 
#'   \emph{Applied Soft Computing} 52: 1009-1019.
#' @export
#' @examples
#' data(fuzzydat)
#' fuzzylm(y ~ x, data = fuzzydat$lee, method = "plrls")
#' \dontrun{
#' # returns an error due to the incorrect number of spreads
#' fuzzylm(y ~ x, data = fuzzydat$dia, method = "fls", fuzzy.left.y = "yl")}
#' # use the same column name for the left and right spread, when the method requests 
#' # non-symmetric fuzzy numbers, but the data specify symmetric fuzzy numbers 
#' fuzzylm(y ~ x, data = fuzzydat$dia, method = "fls", fuzzy.left.y = "yl", fuzzy.right.y = "yl")


fuzzylm = function(formula, data, method = "plrls", fuzzy.left.x = NULL, fuzzy.right.x = NULL, fuzzy.left.y = NULL, fuzzy.right.y = NULL, ...){
	# initiate model.frame
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
    xy = all.vars(mf)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
	x <- attr(mt, "term.labels")
	y <- xy[attr(mt, "response")]
	# check spreads
	if(any(!is.null(c(fuzzy.left.x, fuzzy.right.x, fuzzy.left.y, fuzzy.right.y)))){
		warning("fuzzy spreads detected - assuming same variable order as in formula")
		lhs <- c(y, fuzzy.left.y, fuzzy.right.y)
		rhs <- c(x, fuzzy.left.x, fuzzy.right.x)
		formula <- ifelse(length(lhs) > 1, 
						  sprintf("cbind(%s) ~ %s", toString(lhs), paste(rhs, collapse = " + ")),
						  sprintf("%s ~ %s", toString(lhs), paste(rhs, collapse = " + ")))
		mf <- stats::model.frame(formula, data = data)
	}
	if(attr(stats::terms(mf), "intercept") == 0) stop("fuzzy regression through origin is not supported. Use intercept")
	y <- stats::model.response(mf, "numeric")
	x <- stats::model.matrix(stats::as.formula(formula), data = mf)
	# check method
	methods <- c("plrls", "flar", "fls", "oplr", "moflr", "plr", "diamond", "hung", "lee", "nasrabadi", "tanaka")
	if(!any(grepl(tolower(method), methods)))
		stop(gettextf("method '%s' is not supported.", method))
	method <- methods[methods %in% tolower(sub("-", "", method))]
	coefs <- switch(method, plrls = plrls(x = x, y = y, ...),
	                        flar = flar(x = x, y = y, ...),
							fls = fls(x = x, y = y, ...),
							oplr = oplr(x = x, y = y, ...),
							moflr = moflr(x = x, y = y, ...),
							plr = plr(x = x, y = y, ...),
							diamond = fls(x = x, y = y, ...),
							hung = oplr(x = x, y = x, ...),
							lee = plrls(x = x, y = y, ...),
							nasrabadi = moflr(x = x, y = y, ...),
							tanaka = plr(x = x, y = y, ...))
	fuzzy <- list(call = cl, method = toupper(method), fuzzynum = coefs$fuzzynum, coef = coefs$coef, lims = coefs$lims, x = x, y = y)
	class(fuzzy) <- "fuzzylm"
	fuzzy
}

