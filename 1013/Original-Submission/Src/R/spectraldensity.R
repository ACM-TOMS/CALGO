#' spectraldensity function
#'
#' This function allows calculating the spectral density associated with a given isotropic covariance model
#'
#' @param u Frequency vector
#' @param model Model of spectral density associated with the covariance, variogram or generalized covariance of the target random field
#' @param parameters List of parameter values
#' @return A matrix-valued spectral density value associated with a covariance model, for the chosen frequency vector \eqn{u} in \eqn{R^d}
#' @examples
#' # Bivariate Matern spectral density matrix for a vector u in R^2
#' u <- matrix(rnorm(1), 1, 2)
#' model <- "Matern"
#' # calculate the scale parameters, shape parameters and colocated correlation matrices
#' a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
#' nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
#' C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
#' spectraldensity(u, model, parameters)
#'
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
spectraldensity <- function(u, model, parameters = NULL) {
  # error checking - model name is correct
  if (model %in% c("spherical", "cubic", "penta", "exponential", "Gaussian", "Askey", "Wendland",
                   "Matern", "JBessel", "Laguerre", "Hypergeometric",  "Power", "Spline", 
                   "iterated_exponential", "Gamma", "Cauchy", "stable1", "stable2", 
                   "mixed_Power") == FALSE) {
    stop("Model name invalid", call. = FALSE)
  }
  d <- length(u)

  model <- model
  if(is.null(parameters)) {
    parameters <- NULL
  }
  if (model == "spherical") {
    if(all(parameters$a > 0)) {
      f <- parameters$C*gamma(1+d/2)/(pi^(d/2)*(norm(u,"2"))^d)*(besselJ(2*pi*parameters$a*norm(u,"2")/2, d/2))^2
    }
    else {
      stop("Parameter does not satisfy the model validity restriction")
    }
  }
  if (model == "cubic") {
    if(all(parameters$a > 0)) {
      f = parameters$C*2*gamma(3+d/2)/((pi^(2+d/2)*parameters$a^2*(norm(u,"2"))^(d+2)))*besselJ(((2*pi*parameters$a*norm(u,"2"))/2),1+d/2)^2
    }
    else {
      stop("Parameter does not satisfy the model validity restriction")
    }
  }
  if (model == "penta") {
    if(all(parameters$a > 0)) {
      f = parameters$C*8*gamma(5+d/2)/((3*pi^(4+d/2)*parameters$a^4*(norm(u,"2"))^(d+4)))*besselJ(((2*pi*parameters$a*norm(u,"2"))/2),2+d/2)^2
    }
    else {
      stop("Parameter does not satisfy the model validity restriction")
    }
  }
  if (model == "exponential") {
    if(all(parameters$a > 0)) {
      f <- parameters$C*pi^((d-1)/2)*(2*parameters$a)^d*gamma((d+1)/2)*(1+(2*pi*parameters$a)^2*(norm(u,"2")^2))^(-(d+1)/2)
    }
    else {
      stop("Parameter does not satisfy the model validity restriction")
    }
  }
  if (model == "Gaussian") {
    if(all(parameters$a > 0)) {
      f <- parameters$C*pi^(d/2)*parameters$a^d*exp(-((2*pi*parameters$a)^2*norm(u,"2")^2)/4)
    }
    else {
      stop("Parameter does not satisfy the model validity restriction")
    }
  }
  if (model == "Askey") {
    if (all(parameters$a > 0) & (d==2) & all(parameters$nu1 == 2)){
      f <- parameters$C*2*pi^(d/2)*parameters$a^d*gamma(d)*gamma(parameters$nu1+1)/gamma(d/2)/gamma(parameters$nu1+d+1)*Hygeo1F2((pi*parameters$a*norm(u,"2"))^2,model)
    }
    else {
      stop("Sorry, the Askey model is available only for d=2 and nu=2")
    }
  }
  if (model == "Wendland") {
    k <- parameters$nu2
    if (all(parameters$a > 0) & (d==2) & (k==1) & all(parameters$nu1 == (floor(d/2+k)+1))){
      f <- parameters$C*2^(k+d)*(pi^(d/2))*parameters$a^d*gamma(parameters$nu1+1)*gamma((d+1)/2+k)/sqrt(pi)/gamma(parameters$nu1+d+2*k+1)*Hygeo1F2((pi*parameters$a*norm(u,"2"))^2,model)
    }
    else {
      stop("Sorry, the Wendland model is available only for d=2, nu=3 and k=1")
    }
  }
  if (model == "Matern") {
    if(all(parameters$a > 0) & all(parameters$nu1 > 0)) {
      lnf <- d*log(2*pi*parameters$a) + lgamma(parameters$nu1+d/2) - lgamma(parameters$nu1) - (d/2)*log(pi) - (parameters$nu1+d/2)*log(1+((2*pi*parameters$a)^2)*norm(u,"2")^2)
      f <- parameters$C*exp(lnf)
    }
    else {
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "JBessel") {
    if(all(parameters$a > 0) & all(parameters$nu1 > (d/2-1))) {
      if ((2*pi*parameters$a*norm(u,"2")) < 1){
        f <- parameters$C*pi^(d/2)*(2*parameters$a)^d*gamma(parameters$nu1+1)/gamma(parameters$nu1+1-d/2)*(1-(2*pi*parameters$a)^2*norm(u,"2")^2)^(parameters$nu1-d/2)
      }
      else{
        f = matrix(0,nrow(parameters$C),nrow(parameters$C))
      }
    }
    else {
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "Laguerre") {
    if(all(parameters$a > 0) & all(parameters$nu1 > ((d-1)/4))) {
      f <- parameters$C*gamma(d/2)*(pi*parameters$a*norm(u,"2"))^(2*parameters$nu1)/gamma(parameters$nu1)/pi^(d/2)/(norm(u,"2")^d)*exp(-(pi*parameters$a*norm(u,"2"))^2)
    }
    else {
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "Hypergeometric") {
    if(all(parameters$a > 0) & all(max(parameters$nu1,parameters$nu2) > ((d-1)/4))) {
      f <- parameters$C*2*gamma(d/2)*(pi*parameters$a*(norm(u,"2"))^(parameters$nu1+parameters$nu2))/(gamma(parameters$nu1)*gamma(parameters$nu2)*pi^(d/2)*(norm(u,"2"))^d)*besselK(2*pi*parameters$a*norm(u,"2"),parameters$nu1-parameters$nu2)
    }
    else {
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "Power") {
    a <- c(parameters$a)
    alpha <- c(parameters$nu1)
    k <- floor(alpha/2)
    if(all(a > 0) & all(alpha > (2*k)) & all(alpha < (2*k+2))) {
      f <- parameters$C*matrix(gamma((alpha/2)+1)*gamma((alpha+d)/2)/gamma((alpha/2)-k)/gamma(1-(alpha/2)+k)/(pi^(alpha+(d/2)))/(a^(alpha))/((norm(u,"2"))^(alpha+d)), nrow(parameters$C))
    }
    else {
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "Spline") {
    a <- c(parameters$a)
    k <- c(parameters$nu1)
    if(all(a > 0) & all(k > 0)) {
      f <- parameters$C*matrix(gamma(k+1)*gamma((2*k+d)/2)/2/(pi^(2*k+d/2))/a^(2*k)/((norm(u,"2"))^(2*k+d)), nrow(parameters$C))
    }
    else {
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "iterated_exponential") {
    N <- rpois(1, lambda = 1)
    C <- exp(1)/(N+1)*parameters$C
    a <- 1/(N+1)*parameters$a
    if(all(a > 0)) {
      f <- C*pi^((d-1)/2)*(2*a)^d*gamma((d+1)/2)*(1+(2*pi*a)^2*(norm(u,"2")^2))^(-(d+1)/2)
    }
    else{
      stop("Parameter does not satisfy the model validity restriction")
    }
  }
  if (model == "Gamma") {
    C <- parameters$C
    a <- 1/(rgamma(1,parameters$nu1))*parameters$a
    if(all(a > 0) & all(parameters$nu1 > 0)) {
      f <- C*pi^((d-1)/2)*(2*a)^d*gamma((d+1)/2)*(1+(2*pi*a)^2*(norm(u,"2")^2))^(-(d+1)/2)
    }
    else{
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "Cauchy") {
    C <- parameters$C
    a <- 1/(sqrt(rgamma(1,parameters$nu1)))*parameters$a
    if(all(a > 0) & all(parameters$nu1 > 0)) {
      f <- C*pi^(d/2)*a^d*exp(-((2*pi*a)^2*norm(u,"2")^2)/4)
    }
    else{
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "stable1") {
    C <- parameters$C
    alpha <- parameters$nu1
    W <- rexp(1, rate = 1)
    phi <- runif(1,min = -pi/2, max = pi/2)
    phi0 <- -pi/2
    s <- ((cos(alpha*pi/4))^(2/alpha))*sin(alpha/2*(phi-phi0))/(cos(phi)^(2/alpha))*(cos(phi-alpha/2*(phi-phi0))/W)^(2*(1-alpha/2)/alpha)
    a <- 1/(sqrt(s))*parameters$a
    if(all(a > 0) & all(parameters$nu1 > 0) & all(parameters$nu1 < 2)) {
      f <- C*pi^(d/2)*a^d*exp(-((2*pi*a)^2*norm(u,"2")^2)/4)
    }
    else{
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "stable2") {
    C <- parameters$C
    alpha <- parameters$nu1
    W <- rexp(1, rate = 1)
    phi <- runif(1,min = -pi/2, max = pi/2)
    phi0 <- -pi/2
    s <- ((cos(alpha*pi/2))^(1/alpha))*sin(alpha*(phi-phi0))/(cos(phi)^(1/alpha))*(cos(phi-alpha*(phi-phi0))/W)^((1-alpha)/alpha)
    a <- 1/s*parameters$a
    if(all(a > 0) & all(parameters$nu1 > 0) & all(parameters$nu1 <= 1)) {
      f <- C*pi^((d-1)/2)*(2*a)^d*gamma((d+1)/2)*(1+(2*pi*a)^2*(norm(u,"2")^2))^(-(d+1)/2)
    }
    else{
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  if (model == "mixed_Power") {
    a <- c(parameters$a)
    alpha <- c(parameters$nu1)
    for (i in 1:length(alpha)){
      alpha[i] <- runif(1, min = 0, max = alpha[i])
    }
    if(all(a > 0) & all(alpha > 0) & all(alpha < 2)) {
      k <- floor(alpha/2)
      f <- parameters$C*matrix(gamma((alpha/2)+1)*gamma((alpha+d)/2)/gamma((alpha/2)-k)/gamma(1-(alpha/2)+k)/(pi^(alpha+(d/2)))/(a^(alpha))/((norm(u,"2"))^(alpha+d)), nrow(parameters$C))
    }
    else{
      stop("At least one parameter does not satisfy the model validity restrictions")
    }
  }
  return(f)
}
