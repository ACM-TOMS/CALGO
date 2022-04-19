#' variograms function
#'
#' Valid for simulations over a regular grid in the Euclidean space \eqn{R^d} (with \eqn{d=2} or \eqn{d=3})
#'
#' Comparison of the experimental direct and cross-variograms of \eqn{N} realizations
#' (calculated along one grid axis) with the theoretical direct and cross-variograms
#'
#' This function allows plotting experimental direct and cross-variograms for \eqn{N} realizations (green lines),
#' average of experimental direct and cross-variograms (red lines) and
#' theoretical direct and cross-variograms (black lines)
#'
#' @param simu File containing a matrix with realizations
#' @param nnodes A vector with the number of grid nodes along X and Y directions (if \eqn{d=2}), or X, Y and Z directions (if \eqn{d=3})
#' @param nvar Number of random field components
#' @param N Number of realizations
#' @param dimsimu Along which axis (between \eqn{1} and \eqn{d}) the variograms will be calculated
#' @param hmax The maximum distance (grid meshes) for the calculation of variograms in the direction of the chosen axis
#' @param model Theoretical variogram model
#' @param parameters List of parameter values
#' @param variogname Name of output PNG file
#' @return A PNG file(s) with direct and cross-variograms (experimental, average and theoretical)
#' @examples
#' # Direct and cross-variograms of a bivariate Matérn random field on a regular grid with 500x500 nodes
#' # Thirty realizations (N = 30) are generated
#' set.seed(9784498)
#' nx <- 500
#' ny <- 500
#' coord <- cbind(c(rep(1, times = ny))%x%c(1:nx), (c(1:ny)%x%c(rep(1, times = nx))))
#' model <- "Matern"
#' # calculate the scale parameters, shape parameters and colocated correlation matrices
#' a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
#' nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
#' C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
#' simu <- specsim(coord, a0 = 10, nu0 = 0.5, model, parameters, N = 30, L = 1000,"Matern")
#' # plots of direct and cross-variograms calculated along the abscissa axis for distances ranging from 0 to 120 units
#' variograms(simu, c(nx,ny), nvar=2, N=30, dimsimu=1, hmax=120, model, parameters, "Matern")
#'
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
variograms <- function(simu, nnodes, nvar, N, dimsimu, hmax, model, parameters, variogname){
  if (length(nnodes)==2){
    nz <- 1
  }
  else{
    nz <- nnodes[3]
  }
  nx <- nnodes[1]
  ny <- nnodes[2]
  delta <- 0.1
  if (nvar==1){
    simu <- array(simu, dim = c(nx, ny, nz, N))
    variog = matrix(0, nrow = hmax, ncol = N)
    if (dimsimu==1){
      for(h in 1:hmax){
        d = array(simu[1:(nx-h), , , ]-simu[(1+h):nx, , , ],dim = c((nx-h)*ny*nz, N))
        variog[h, ] = 0.5*colMeans(x=d*d, na.rm = TRUE)
      }
    }
    else if (dimsimu==2){
      for(h in 1:hmax){
        d = array(simu[, 1:(ny-h), , ]-simu[, (1+h):ny, , ],dim = c(nx*(ny-h)*nz,N))
        variog[h, ] = 0.5*colMeans(x=d*d, na.rm = TRUE)
      }
    }
    else if (dimsimu==3){
      for(h in 1:hmax){
        d = array(simu[, , 1:(nz-h), ]-simu[, , (1+h):nz, ], dim = c(nx*ny*(nz-h), N))
        variog[h, ] = 0.5*colMeans(x=d*d, na.rm = TRUE)
      }
    }
    filename <- as.character(paste(variogname, "variograms", sep = "_"))
    png(paste(filename,".png",sep =''), width=1300, height=1300, res=250)
    par(mfrow=c(1,1))
    h0 = seq(1:hmax)
    h1 = seq(from=0.01*hmax, to=hmax, by=0.01*hmax)
    plot(h0,variog[ ,1], col="green", xlab="Lag separation distance", ylab="Variogram", type="l", 
         lty=2, lwd = c(1, 2, 2), ylim = c(min(variog[ ,1]), max(variog[ ,1])+delta))
    for (i in 2:N){
      lines(h0, variog[ ,i], col="green", type="l", lty=2)
    }
    lines(h0,colMeans(x=t(variog), na.rm = TRUE), col="red", type = "l", lwd = 2)
    cov = parameters$C*cova(model, h1/(parameters$a), parameters$nu1)
    lines(h1,cov[1]-cov, col="black", type="l", lwd = 2)
    legend("bottomright", inset = 0.01, legend=c("Individual realizations", "Average over realizations", 
                                                 "Theoretical model"), col=c("green", "red", "black"), 
           lty=c(1, 1, 1), lwd = c(1, 2, 2), cex=0.85, box.lty=0)
    dev.off()
  }
  else{
    simus <- array(0, dim=c(nx, ny, nz, N, nvar))
    for (i in 1:nvar){
      simus[ , , , , i] <- array(simu[ , seq(i, nvar*N, by=nvar)], dim = c(nx, ny, nz, N))
    }
    variogs <- array(0, dim=c(hmax, N, nvar+choose(n=nvar, k=2)))
    if (dimsimu==1){
      for(h in 1:hmax){
        d = array(0, dim=c((nx-h)*ny*nz, N, nvar))
        for (i in 1:nvar){
          d[ , ,i] = array(simus[1:(nx-h), , , ,i]-simus[(1+h):nx, , , ,i], dim = c((nx-h)*ny*nz, N))
          variogs[h, ,i] = 0.5*colMeans(x=d[, , i]*d[, , i], na.rm = TRUE)
        }
        aux = 1
        for (k in 1:(nvar-1)){
          for (j in (k+1):nvar){
            variogs[h, ,nvar+aux] = 0.5*colMeans(x=d[, , k]*d[, , j], na.rm = TRUE)
            aux = aux + 1
          }
        }
      }
    }
    else if(dimsimu==2){
      for(h in 1:hmax){
        d = array(0, dim=c(nx*(ny-h)*nz, N, nvar))
        for (i in 1:nvar){
          d[ , ,i] = array(simus[ , 1:(ny-h), , ,i]-simus[ , (1+h):ny, , ,i],dim = c(nx*(ny-h)*nz, N))
          variogs[h, ,i] = 0.5*colMeans(x=d[, , i]*d[, , i], na.rm = TRUE)
        }
        aux = 1
        for (k in 1:(nvar-1)){
          for (j in (k+1):nvar){
            variogs[h, , nvar+aux] = 0.5*colMeans(x=d[, , k]*d[, , j], na.rm = TRUE)
            aux = aux + 1
          }
        }
      }
    }
    else if(dimsimu==3){
      for(h in 1:hmax){
        d = array(0, dim=c(nx*ny*(nz-h), N, nvar))
        for (i in 1:nvar){
          d[ , ,i] = array(simus[ , , 1:(nz-h), ,i]-simus[ , , (1+h):nz, ,i], dim = c(nx*ny*(nz-h), N))
          variogs[h, ,i] = 0.5*colMeans(x=d[, , i]*d[, , i], na.rm = TRUE)
        }
        aux = 1
        for (k in 1:(nvar-1)){
          for (j in (k+1):nvar){
            variogs[h, ,nvar+aux] = 0.5*colMeans(x=d[, , k]*d[, , j], na.rm = TRUE)
            aux = aux + 1
          }
        }
      }
    }
    for (k in 1:nvar){
      h0 = seq(1:hmax)
      h1 = seq(from=0.01*hmax, to=hmax, by=0.01*hmax)
      filename <- as.character(paste(variogname, "direct-variograms", k, sep = "_"))
      png(paste(filename,".png", sep =''), width=1300, height=1300, res=250)
      plot(h0, variogs[ , 1, k], col="green", xlab="Lag separation distance", ylab="Direct variogram", 
           type="l", lty=2, lwd = c(1, 2, 2), ylim = c(0, parameters$C[k, k]+delta))
      for (i in 2:N){
        lines(h0, variogs[ , i, k], col="green", type="l", lty=2)
      }
      lines(h0,colMeans(x=t(variogs[, , k]), na.rm = TRUE),col="red", type = "l", lwd = 2)
      cov = parameters$C[k, k]*cova(model, h1/(parameters$a[k, k]), parameters$nu1[k, k])
      lines(h1, cov[1]-cov, col="black", type="l", lwd = 2)
      legend("bottomright", inset = 0.01, legend=c("Individual realizations", "Average over realizations", 
                                                   "Theoretical model"), col=c("green", "red", "black"), 
             lty=c(1 ,1, 1), lwd = c(1, 2, 2), cex=0.85, box.lty=0)
      dev.off()
    }
    aux = 1
    for (k in 1:(nvar-1)){
      for (j in (k+1):nvar){
        h0 = seq(1:hmax)
        h1 = seq(from=0.01*hmax, to=hmax, by=0.01*hmax)
        filename <- as.character(paste(variogname, "cross-variograms",k, j, sep = "_"))
        png(paste(filename, ".png", sep =''), width=1300, height=1300, res=250)
        if (parameters$C[k, j] > 0){
          loc.legend ="bottomright"
          plot(h0, variogs[ , 1, nvar+aux], col="green", xlab="Lag separation distance", 
               ylab="Cross-variogram", type="l", lty=2, lwd = c(1, 2, 2), 
               ylim = c(0, parameters$C[k, j]+delta))
        }
        else{
          loc.legend ="topright"
          plot(h0, variogs[ , 1, nvar+aux], col="green", xlab="Lag separation distance", 
               ylab="Cross-variogram", type="l", lty=2, lwd = c(1, 2, 2), 
               ylim = c(parameters$C[k, j]-delta,0))
        }
        for (i in 2:N){
          lines(h0, variogs[ , i, nvar+aux], col="green", type="l", lty=2)
        }
        lines(h0, colMeans(x=t(variogs[ , , nvar+aux]), na.rm = TRUE), col="red", type = "l", lwd = 2)
        cov = parameters$C[k, j]*cova(model, h1/(parameters$a[k, j]), parameters$nu1[k, j])
        lines(h1, cov[1]-cov, col="black", type="l", lwd = 2)
        legend(loc.legend, inset = 0.01, legend=c("Individual realizations", "Average over realizations", 
                                                  "Theoretical model"), col=c("green", "red", "black"), 
               lty=c(1, 1, 1), lwd = c(1, 2, 2), cex=0.85, box.lty=0)
        aux = aux + 1
        dev.off()
      }
    }
  }
}
