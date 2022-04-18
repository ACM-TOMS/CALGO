#' Color scale representation of a realization
#'
#' This function allows displaying a color scale representation of a simulated random field.
#' Valid for simulation over a regular grid in the Euclidean space \eqn{R^d} with \eqn{d=2} or \eqn{d=3}
#'
#' @param simu File with one realization of the simulated random field
#' @param viewsimu View at the XY, XZ, or YZ plane: "XY", "XZ" or "YZ"
#' @param coord Coordinates of the locations where the random field has been simulated
#' @param nnodes A vector with the number of grid nodes along X and Y directions (\eqn{d=2}), or X, Y and Z directions (\eqn{d=3})
#' @param gridmesh A vector with the grid mesh along X and Y directions (if \eqn{d=2}: c(dx,dy)), or X, Y and Z directions (if \eqn{d=3}: c(dx,dy,dz))
#' @param slice Index of the grid slice to display
#' @param mapname Name of PNG output file
#' @param width.png Display width depending on whether the grid is square or rectangular. By default, the value is for a square grid.
#' @param height.png Display height depending on whether the grid is square or rectangular. By default, the value is for a square grid.
#' @return A PNG file(s) with a map of one realization of the random field
#' @examples
#' # A map of one realization of a simulated bivariate Matérn random field on a regular grid 
#' # with 500x500 nodes
#' set.seed(9784498)
#' nx <- 500
#' ny <- 500
#' coord <- cbind(c(rep(1, times = ny))%x%c(1:nx), (c(1:ny)%x%c(rep(1, times = nx))))
#' model <- "Matern"
#' # Calculating the scale parameters, shape parameters and colocated correlation matrices
#' a <- matrix(c(20, (0.5*(20^-2+100^-2))^-0.5, (0.5*(20^-2+100^-2))^-0.5, 100), 2, 2)
#' nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
#' C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' # dimension
#' d <- 2
#' # Checking the validity of the parameters
#' validationMatern(a, nu, C, d)
#' # Building a list of parameters of the simulation
#' parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
#' # Generating one realization
#' simu <- specsim(coord, a0 = 10, nu0 = 0.5, model, parameters, N = 1, L = 1000, "Matern")
#' # Display of realization of simulated bivariate random field
#' maps(simu, viewsimu="XY", coord, c(nx,ny), c(1,1), slice=1, "Matern")
#'
#' @author Daisy Arroyo, \email{darroyof@@udec.cl} and Xavier Emery, \email{xemery@@ing.uchile.cl}
#' @export
#'
maps <- function(simu, viewsimu, coord, nnodes, gridmesh, slice, mapname, width.png=1250, height.png=1370){
  if(is.null(ncol(simu))) {
    nvar <- 1
  }
  else{
    nvar <- ncol(simu)
  }
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  if (length(nnodes)==2){
    nx <- nnodes[1]
    ny <- nnodes[2]
    nz <- 1
    dx <- gridmesh[1]
    dy <- gridmesh[2]
    dz <- 1
    coord <- cbind(coord,1)
  }
  else{
    nx <- nnodes[1]
    ny <- nnodes[2]
    nz <- nnodes[3]
    dx <- gridmesh[1]
    dy <- gridmesh[2]
    dz <- gridmesh[3]
  }
  mincoord <- apply(coord, MARGIN=2, FUN=min)
  maxcoord <- apply(coord, MARGIN=2, FUN=max)
  x0 <- mincoord[1]
  y0 <- mincoord[2]
  z0 <- mincoord[3]
  x <- seq(from = x0, to = (dx*nx), by = dx)
  y <- seq(from = y0, to = (dy*ny), by = dy)
  z <- seq(from = z0, to = (dz*nz), by = dz)
  #
  # Displaying a realization of the simulated components over a grid
  #=================================================================
  gridsimu  <- array(simu, dim = c(nx, ny, nz, nvar))
  for (i in 1:nvar){
    if (viewsimu == "XY"){
      filename <- as.character(paste(mapname, "map", viewsimu, paste("variable", i, sep = ""), paste("slice", slice, sep = ""), sep = "_"))
      png(paste(filename, ".png", sep =''), width=width.png, height=height.png, pointsize = 10, res = 300)
      A <- array(gridsimu[ , , slice, i], dim = c(nx, ny))
      image(x, y, A, col = jet.colors(nx*ny), xlim = c(x0, dx*nx), ylim = c(y0, dy*ny), xlab="", ylab="")
      box()
      dev.off()
    }
    if (viewsimu == "XZ"){
      filename <- as.character(paste(mapname, "map",viewsimu, paste("variable", i, sep = ""), paste("slice", slice, sep = ""), sep = "_"))
      png(paste(filename, ".png", sep =''), width = width.png, height = height.png, units = 'in', pointsize = 10, res = 300)
      A <- array(gridsimu[ , slice, , i], dim = c(nx, nz))
      image(x, z, A, col = jet.colors(nx*nz), xlim = c(x0, dx*nx), ylim = c(z0, dy*nz), xlab="", ylab="")
      box()
      dev.off()
    }
    if (viewsimu == "YZ"){
      filename <- as.character(paste(mapname, "map",viewsimu, paste("variable", i, sep = ""), paste("slice", slice, sep = ""), sep = "_"))
      png(paste(filename,".png",sep =''), width = width.png, height = height.png, units = 'in', pointsize = 10, res = 300)
      A <- array(gridsimu[slice, , , i], dim = c(ny, nz))
      image(y, z, A, col = jet.colors(ny*nz), xlim = c(y0, dy*ny), ylim = c(z0, dz*nz), xlab="", ylab="")
      box()
      dev.off()
    }
  }
}
