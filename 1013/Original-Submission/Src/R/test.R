#' test function
#'
#' This function allows testing the distribution of the generalized increments of simulated random fields
#' with Shapiro-Wilk normality test, for selected bivariate models
#'
#' @param simu An array 9x(2xN)x2, containing the two matrix simulations with L=100 and L=1000
#' @param coord1 First set of locations
#' @param coord2 Second set of locations
#' @param lambda Set of weights
#' @param model Covariance, variogram or generalized covariance model of the target random field
#' @param parameters List of parameter values
#' @param N Number of realizations
#' @param filename Name of PNG output file
#' @return A probability-probability (P-P) plot:
#' proportions of times when the test is accepted (value 0) and when
#' the test is rejected (value 1) are counted and compared to the significance levels
#' @examples
#' # A particular case for a simulated bivariate Matérn random field with L=100 and L=1000
#' L <- c(100, 1000)
#' N <- 10000
#' coord <- matrix(c(1, 1, 11, 1, 1, 11, 11, 11, 6, 1, 1, 6, 11, 6, 6, 11, 6, 6), nrow=9, ncol=2, byrow=TRUE)
#' model <- "Matern"
#' # calculate the scale parameters, shape parameters and colocated correlation matrices
#' a <- matrix(c(5, (0.5*(5^-2+10^-2))^-0.5, (0.5*(5^-2+10^-2))^-0.5, 10), 2, 2)
#' nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
#' C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' # dimension
#' d <- 2
#' # Checking the validity of the parameters of bivariate Matérn model
#' validationMatern(a, nu, C, d)
#' # Building a list of parameters to the simulation
#' parameters.SRF2 <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
#' # Simulation with L=100 and L=1000, respectively:
#' simu.SRF2.L100 <- specsim(coord, a0 = 15, nu0 = 0.5, model, parameters.SRF2, N, L[1])
#' simu.SRF2.L1000 <- specsim(coord, a0 = 15, nu0 = 0.5, model, parameters.SRF2, N, L[2])
#' # Put the two simulations together in an array
#' simu.SRF2 <- array(c(simu.SRF2.L100, simu.SRF2.L1000), dim = c(9, 2*N, 2))
#' # Preparing to the test
#' coord1 <- coord[c(1,2,3,4,9), ]
#' coord2 <- coord[5:9, ]
#' lambda <- c(0.25, 0.25, 0.25, 0.25, -1)
#' # Test and generate a png file named "test_SRF2"
#' test(simu.SRF2, coord1, coord2, lambda, model, parameters.SRF2, N, "test_SRF2")
#'
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
test <- function(simu, coord1, coord2, lambda, model, parameters, N, filename="test"){
  rho <- parameters$C[1, 2]
  simu1_100 <- simu[ , seq(1, 2*N, by=2), 1] # First simulated random variable with L=100
  simu2_100 <- simu[ , seq(2, 2*N, by=2), 1] # Second simulated random variable with L=100
  simu1_1000 <- simu[ , seq(1, 2*N, by=2), 2] # First simulated random variable with L=1000
  simu2_1000 <- simu[ , seq(2, 2*N, by=2), 2] # Second simulated random variable with L=1000
  m = nrow(coord1)
  t1 = coord1
  t2 = coord2
  t1 = tcrossprod(t1, y = NULL)
  t2 = tcrossprod(t2, y = NULL)
  h1 = sqrt(pmax(matrix(0, m, m), -2*t1+matrix(diag(t1), m, m)+matrix(t(diag(t1)), m, m, byrow = TRUE))) # distance matrix for coord1
  h2 = sqrt(pmax(matrix(0, m, m), -2*t2+matrix(diag(t2), m, m)+matrix(t(diag(t2)), m, m, byrow = TRUE))) # distance matrix for coord2
  h11 = matrix(coord1[ ,1], m, m)-matrix(t(coord2[ ,1]), m, m, byrow = TRUE) # distance matrix for coord1 and coord2
  h22 = matrix(coord1[ ,2], m, m)-matrix(t(coord2[ ,2]), m, m, byrow = TRUE) # distance matrix for coord1 and coord2
  h12 = sqrt(pmax(matrix(0, m, m),(h11^2+h22^2)))  # distance matrix between coord1 and coord2

  # Variances and covariance considering the set of weights lambda
  # ==============================================================
  a <- parameters$a
  b <- parameters$nu1
  K11 = crossprod(lambda,crossprod(t(cova(model,h1/a[1,1],b[1,1])),cbind(c(lambda))))  # Equation (15)
  K22 = crossprod(lambda,crossprod(t(cova(model,h2/a[2,2],b[2,2])),cbind(c(lambda))))  # Equation (15)
  K12 = crossprod(lambda,crossprod(t(cova(model,h12/a[1,2],b[1,2])),cbind(c(lambda)))) # Equation (16)

  # Set of variances and covariances
  # ================================
  covariances <- matrix(c(K11, K12, K12, K22), 2, 2, byrow = TRUE)

  # Test for the combination of first and second components, with L = 100
  # =====================================================================
  A1 = crossprod(lambda,simu1_100[c(1,2,3,4,9), ])
  A2 = crossprod(lambda,simu2_100[c(5,6,7,8,9), ])
  B1 = matrix(A1,100,(N/100))
  B2 = matrix(A2,100,(N/100))
  C = (B1+B2)/sqrt(covariances[1,1]+covariances[2,2]+2*rho*covariances[1,2])

  t_100_11 <- rep(0,19)
  for (k in 1:19){
    t11_v12 = 0
    for(i in 1:(N/100)){
      h = shapiro.test(C[ ,i])
      if(h$p.value < (k/20)){
        t11_v12 = t11_v12 + 1
      }
    }
    t_100_11[k] = t11_v12
  }

  # Test for first component with  L = 100
  # ======================================
  B1 = B1/sqrt(covariances[1,1])
  t_100_v1_1 = rep(0,19)
  for (k in 1:19){
    t1_v1 = 0
    for(i in 1:(N/100)){
      h = shapiro.test(B1[ ,i])
      if(h$p.value < (k/20)){
        t1_v1 = t1_v1 + 1
      }
    }
    t_100_v1_1[k] = t1_v1
  }

  # Test for second component with L = 100
  # ======================================
  B2 = B2/sqrt(covariances[2,2])
  t_100_v2_1 = rep(0,19)
  for (k in 1:19){
    t1_v2 = 0
    for(i in 1:(N/100)){
      h = shapiro.test(B2[ ,i])
      if(h$p.value < (k/20)){
        t1_v2 = t1_v2 + 1
      }
    }
    t_100_v2_1[k] = t1_v2
  }

  # Test for the combination of first and second components, with L = 1000
  # ======================================================================
  A1 = crossprod(lambda,simu1_1000[c(1,2,3,4,9), ])
  A2 = crossprod(lambda,simu2_1000[c(5,6,7,8,9), ])
  B1 = matrix(A1,100,(N/100))
  B2 = matrix(A2,100,(N/100))
  C = (B1+B2)/sqrt(covariances[1,1]+covariances[2,2]+2*rho*covariances[1,2])

  t_1000_11 <- rep(0,19)
  for (k in 1:19){
    t11_v12 = 0
    for(i in 1:(N/100)){
      h = shapiro.test(C[ ,i])
      if(h$p.value < (k/20)){
        t11_v12 = t11_v12 + 1
      }
    }
    t_1000_11[k] = t11_v12
  }

  # Test for first component with L = 1000
  # ======================================
  B1 = B1/sqrt(covariances[1,1])
  t_1000_v1_1 = rep(0,19)
  for (k in 1:19){
    t1_v1 = 0
    for(i in 1:(N/100)){
      h = shapiro.test(B1[ ,i])
      if(h$p.value < (k/20)){
        t1_v1 = t1_v1 + 1
      }
    }
    t_1000_v1_1[k] = t1_v1
  }

  # Test for second component with L = 1000
  # =======================================
  B2 = B2/sqrt(covariances[2,2])
  t_1000_v2_1 = rep(0,19)
  for (k in 1:19){
    t1_v2 = 0
    for(i in 1:(N/100)){
      h = shapiro.test(B2[ ,i])
      if(h$p.value < (k/20)){
        t1_v2 = t1_v2 + 1
      }
    }
    t_1000_v2_1[k] = t1_v2
  }

  # Plots:
  # =====
  png(paste(filename,".png",sep =''), width=200, height=100,units = 'mm', res=300)
  x <- seq(0.05,0.95, 0.05)
  limsup <- c(0.09,0.15,0.21,0.27,0.32,0.38,0.43,0.48,0.53,0.58,0.63,0.68,0.73,0.77,0.82,0.86,0.91,0.95,0.98)
  liminf <- c(0.02,0.05,0.09,0.14,0.18,0.23,0.27,0.32,0.37,0.42,0.47,0.52,0.57,0.62,0.68,0.73,0.79,0.85,0.91)
  par(mfrow=c(1,3))
  par(pty="s")
  plot(x, t_100_v1_1/((N/100)), col="blue", pch="o",xlab="Nominal level of S-W test",ylab="Proportion of realizations",xlim = c(0,1),ylim = c(0,1))
  lines(x, t_100_v1_1/(N/100),col="blue",lty=1)
  points(x, t_1000_v1_1/(N/100), col="red", pch="o")
  lines(x, t_1000_v1_1/(N/100),col="red",lty=1)
  lines(x,liminf,col="green",lty=1,lwd = 1.5)
  lines(x,limsup,col="green",lty=1,lwd = 1.5)
  abline(a=0, b=1, lwd = 2)
  legend("topleft", inset = 0.01, legend=c("L=100","L=1000"), col=c("blue","red"),pch=c("o","o"),lty=c(1,1),cex=0.85, box.lty=0)

  plot(x, t_100_v2_1/(N/100), col="blue", pch="o", xlab="Nominal level of S-W test",ylab="Proportion of realizations",xlim = c(0,1),ylim = c(0,1))
  lines(x, t_100_v2_1/(N/100),col="blue",lty=1)
  points(x, t_1000_v2_1/(N/100), col="red", pch="o")
  lines(x, t_1000_v2_1/(N/100),col="red",lty=1)
  lines(x,liminf,col="green",lty=1,lwd = 1.5)
  lines(x,limsup,col="green",lty=1,lwd = 1.5)
  abline(a=0, b=1, lwd = 2)
  legend("topleft", inset = 0.01, legend=c("L=100","L=1000"), col=c("blue","red"),pch=c("o","o"),lty=c(1,1),cex=0.85, box.lty=0)

  plot(x, t_100_11/(N/100), col="blue", pch="o",xlab="Nominal level of S-W test",ylab="Proportion of realizations",xlim = c(0,1),ylim = c(0,1))
  lines(x, t_100_11/(N/100),col="blue",lty=1)
  points(x, t_1000_11/(N/100), col="red", pch="o")
  lines(x, t_1000_11/(N/100),col="red",lty=1)
  lines(x,liminf,col="green",lty=1,lwd = 1.5)
  lines(x,limsup,col="green",lty=1,lwd = 1.5)
  abline(a=0, b=1, lwd = 2)
  legend("topleft", inset = 0.01, legend=c("L=100","L=1000"), col=c("blue","red"),pch=c("o","o"),lty=c(1,1),cex=0.85, box.lty=0)
  dev.off()
}
