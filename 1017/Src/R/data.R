#' Data For Fuzzy Linear Regression
#'
#' Example data reported by the authors of the respective fuzzy linear regression methods
#' for testing model fit performance.
#' 
#' @docType data
#' @usage data(fuzzydat)
#' @format A list of data.frames.
#' @keywords datasets
#' @source Diamond, P. (1988) Fuzzy least squares. \emph{Information Sciences}
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
#' @examples
#' data(fuzzydat)
#' fuzzylm(y ~ x, data = fuzzydat$lee)
#' fuzzylm(y ~ x, data = fuzzydat$dia, method = "fls", fuzzy.left.y = "yl", fuzzy.right.y = "yl")
"fuzzydat"

#' Temperature Data of Hibernating Bats and Climate at Site
#'
#' Body surface temperature of multiple species of hibernating bats and mean annual 
#' surface temperature at the hibernation site.
#'
#' @docType data
#' @usage data(bats)
#' @format A data frame with 528 rows and two variables:
#' \describe{
#'   \item{\code{MAST}}{numeric Mean annual surface temperature at the site in degrees 
#'     Celsius}
#'   \item{\code{temperature}}{numeric Body surface temperature of hibernating bats in 
#'     degrees Celsius}
#' }
#' @keywords datasets
#' @source Martinkova, N., Pikula, J., Zukal, J., Kovacova, V., Bandouchova, H., 
#'   Bartonicka, T., Botvinkin, A.D., Brichta, J., Dundarova, H., Kokurewicz, T., 
#'   Irwin, N.R., Linhart, P., Orlov, O.L., Piacek, V., Skrabanek, P., Tiunov, 
#'   M.P. and Zahradnikova, A., Jr. (2018) Hibernation temperature-dependent 
#'   \emph{Pseudogymnoascus destructans} infection intensity in Palearctic bats. 
#'   \emph{Virulence} 9: 1734-1750.
#' 
#'   Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G. and Jarvis, A. (2005) Very 
#'   high resolution interpolated climate surfaces for global land areas. 
#'   \emph{International Journal of Climatology} 25: 1965-1978.
#' @examples
#' data(bats)
#' # remove outlier
#' dat <- bats[!(bats$MAST < 0 & bats$temperature > 7), ]
#' # fuzzy linear regression model as published
#' fit <- fuzzylm(temperature ~ MAST, data = dat, method = "plrls", h = 0.01, k1 = 5)
#' plot(fit, res = 30, col = "orange")
"bats"
