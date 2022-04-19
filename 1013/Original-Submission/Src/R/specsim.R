#' Continuous spectral simulation of random fields
#'
#' This function allows simulating scalar and vector random fields defined in an Euclidian space by a continuous spectral algorithm
#'
#' @param coord Coordinates of the locations targeted for simulation
#' @param a0 Scale parameter of the Matérn spectral density g used to simulate the frequency vectors
#' @param nu0 Shape parameter of the Matérn spectral density g
#' @param model Spectral density associated with the covariance, variogram or generalized covariance of the target random field
#' @param parameters List of parameter values
#' @param N Number of realizations to generate
#' @param L Number of basic random fields to use
#' @param filename Name of output txt file
#' @return An \eqn{n*(N*P)} matrix containing the N simulated values for each coordinate (total = \eqn{n}) and for each random field component (P)
#' @examples
#' # Simulation of a bivariate Matérn random field on a regular grid with 500x500 nodes
#' set.seed(9784498)
#' nx <- 500
#' ny <- 500
#' coord <- cbind(c(rep(1, times = ny))%x%c(1:nx), (c(1:ny)%x%c(rep(1, times = nx))))
#' model <- "Matern"
#' # calculating the scale parameters, shape parameters and colocated correlation matrices
#' a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
#' nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
#' C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' # dimension
#' d <- 2
#' # Checking the validity of the bivariate Matérn model parameters
#' validationMatern(a, nu, C, d)
#' # Building a list of parameters for simulation
#' parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
#' # Simulation and generation of a txt file with name "Matern"
#' simu <- specsim(coord, a0 = 10, nu0 = 0.5, model, parameters, N = 1, L = 1000, "Matern")
#' # Display of realization of simulated bivariate random field
#' maps(simu, viewsimu="XY", coord, c(nx,ny), c(1,1), slice=1, "Matern_map")
#'
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
specsim <- function(coord, a0, nu0, model, parameters = NULL, N, L, filename="specsim"){
  # error checking - model name is correct
  if (model %in% c("spherical", "cubic", "penta", "exponential", "Gaussian", "Askey", "Wendland",
                   "Matern", "JBessel", "Laguerre", "Hypergeometric",  "Power", "Spline", 
                   "iterated_exponential", "Gamma", "Cauchy", "stable1", "stable2", 
                   "mixed_Power") == FALSE) {
    stop("Model name invalid", call. = FALSE)
  }
  model <- model
  if(is.null(parameters)) {
    parameters <- NULL
  }
  parameters <- parameters
  print("")
  print("Preparing simulation: calculating parameters of the cosine and sine waves")
  print("")
  n <- dim(coord)[1]
  d <- dim(coord)[2]
  f <- spectraldensity(cbind(rep(0, times = d)), model, parameters)
  P <- NROW(f)
  if (max(is.infinite(f))==TRUE){
    vtype = 1}
  else{
    vtype = 0}
  
  # Generate random frequencies and random phases
  # =============================================
  
  S <- ceiling(1e7*runif(3))
  set.seed(S[1])
  G <- matrix(rgamma(P*L*N, nu0, scale = 1),P*L*N,d)
  set.seed(S[2])
  u <- matrix(rnorm(P*L*N*d), P*L*N, d)/sqrt(G*2)/a0/(2*pi)
  set.seed(S[3])
  phi <- 2*pi*runif(P*L*N)
  
  # Calculate alpha and beta weights
  # ================================
  A <- matrix(0, P, P*L*N)
  B <- matrix(0, P, P*L*N)
  k <- 0
  neg1 <- NULL
  parametersg <- list("C" = 1, "a" = a0, "nu1" = nu0, "nu2" = c())
  for(i in 1:(P*L*N)){
    g <- spectraldensity(u[i, ], "Matern", parametersg)
    f <- spectraldensity(u[i, ], model, parameters)
    VD <- eigen(2*f/g)
    neg <- which(VD$values<0)
    if(length(neg)>0){
      k <- k + 1
      neg1[k] <- min(VD$values[neg])
      VD$values[neg] = 0
    }
    if(P==1){
      temp <- crossprod(t(VD$vectors),tcrossprod(sqrt(VD$values),VD$vectors))
    }
    else{
      temp <- crossprod(t(VD$vectors),tcrossprod(sqrt(diag(VD$values)),VD$vectors))
    }
    j <- i - P*floor((i-1)/P)
    A[ ,i] <- Re(temp[ ,j])
    B[ ,i] <- -Im(temp[ ,j])
  }
  
  # Simulation at target locations
  # ==============================
  print("Simulation at target locations")
  progress <- 0
  sequen <- c(seq(0,n-0.5, by = ceiling(1e5/P/N)),n)
  for(n1 in 1:(length(sequen)-1)){
    coord_n <- coord[(sequen[n1]+1):sequen[n1+1], ]
    m <- dim(coord_n)[1]
    simu <- array(rep(0,m*P), c(m,P,N))
    for(k in 1:N){
      index <- cbind(((k-1)*L*P+1):(k*L*P))
      x0 <- 2*pi*tcrossprod(coord_n,u[index, ])
      x <- x0 + tcrossprod(cbind(rep(1,m)),phi[index])
      if(vtype == 0){
        if(P==1){
          simu[ , ,k] <- crossprod(t(cos(x)),A[index]) + crossprod(t(sin(x)),B[index])
        }
        else{
          simu[ , ,k] <- tcrossprod(cos(x),A[ ,index]) + tcrossprod(sin(x),B[ ,index])
        }
      }
      else {
        cosphi <- cos(tcrossprod(cbind(rep(1,m)),phi[index]))
        sinphi <- sin(tcrossprod(cbind(rep(1,m)),phi[index]))
        if(P==1){
          simu[ , ,k] <- crossprod(t(cos(x)-cosphi),A[index]) + crossprod(t(sin(x)-sinphi),B[index])
        }
        else{
          simu[ , ,k] <- tcrossprod(cos(x)-cosphi,A[ ,index]) + tcrossprod(sin(x)-sinphi,B[ ,index])
        }
      }
    }
    simu <- matrix(simu, m, P*N)/sqrt(L)
    if(n1 == 1)
    {simu1 = simu}
    else{
      simu1 <- rbind(simu1, simu)
    }
    
    # Report on progress from time to time
    # ====================================
    progress2 <- 10*floor((10*n1)/(length(sequen)-1))
    if(progress2 > progress){
      report <- paste(c("", as.character(progress2), "% completed"), collapse = " ")
      print(report)
      progress <- progress2
      Sys.sleep(0.001)
    }
  }
  if(length(neg1)>0){
    val_neg <- paste(c("Negative eigenvalues have been found for the spectral density matrix and have been set to zero for simulation. This may be due to numerical imprecision or to the fact that the covariance model is mathematically invalid."), collapse = " ")
    print(val_neg)
    min_value <- paste(c("The minimum eigenvalue is", sprintf(min(neg1), fmt = '%#e') ), collapse = " ")
    print(min_value)
  }
  print("Saving the simulation to an output txt file")
  write.table(simu1, file = paste(filename,".txt",sep =''), sep = " ", dec = ".", row.names = FALSE, col.names=FALSE)
  print("...Done...")
  return(simu1)
}