# ============================================================================#
# specsim: An R Package for Simulating Vector Gaussian Random Fields Using a  #
#          Continuous Spectral Algorithm                                      #
#                                                                             #
# The specsim package is not available from CRAN, but can easily be installed #
# using the following R code:                                                 #
# install.packages("specsim_1.0.tar.gz",repos=NULL, type = "source")          #
# library("specsim")                                                          #
#                                                                             #
# This R script can be used to replicate the results of the paper.            #
# Version: 01 Sep 2020                                                        #
# Authors: D. Arroyo and X. Emery                                             #
# License: GPL-2 (http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)   #
# ============================================================================#

rm(list = ls())

# loading the package
install.packages("specsim_1.0.tar.gz",repos=NULL,type = "source")
library("specsim")

##### ======================================================
##### A worked example of the continuous spectral algorithm
##### Subsection Example
##### ======================================================

# Coordinates of the locations in R^2 targeted for simulation
set.seed(9784498)
nx <- 1000
ny <- 1000
coord <- cbind(c(rep(1, times = ny))%x%c(1:nx), (c(1:ny)%x%c(rep(1, times = nx))))

# Simulation of a trivariate Matérn random field
# ==============================================
# Fix scale and shape parameters for a density g
a0 <- 100
nu0 <- 0.5
# Building parameters of a valid trivariate Matérn model
model <- "Matern"
a <- matrix(c(400, (0.5*(400^-2+100^-2))^-0.5, (0.5*(400^-2+50^-2))^-0.5,
              (0.5*(400^-2+100^-2))^-0.5, 100,  (0.5*(100^-2+50^-2))^-0.5,
              (0.5*(400^-2+50^-2))^-0.5, (0.5*(100^-2+50^-2))^-0.5, 50), 3, 3)
nu <- matrix(c(0.5, 0.75, 1, 0.75, 1, 1.25, 1, 1.25, 1.5), 3, 3)
C <- matrix(c(1, 0.5, 0.33, 0.5, 1, 0.25, 0.33, 0.25, 1), 3, 3)
d <- 2
# Checking the validity of the built parameter matrices
validationMatern(a, nu, C, d)
# Building a list of parameters for the simulation
parameters <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
# Fifty realizations (N=50) will be generated
N <- 50
# Number of basic random fields
L <- 1000
# Name of output txt file with simulations
filename <- "Matern"
# Performing the simulation
simu <- specsim(coord, a0, nu0, model, parameters, N, L, filename)


# A display of one realization of the three simulated components with the "maps" function (Figure 1)
# ==================================================================================================
# Consider the first realizations, that is the three first columns of matrix simu
simu_first <- simu[ ,1:3]
maps(simu_first, viewsimu="XY", coord, c(nx,ny), c(1,1), slice=1, "Figure_1")

# A plot of experimental and theoretical direct and cross-variograms along the abscissa axis with the 
# ===================================================================================================
# "variograms" function (Figure 2)
# ================================
variograms(simu, c(nx,ny), nvar=3, N, dimsimu=1, hmax=500, model, parameters, "Figure_2")


##### ============================================
##### Plot time vs number of grid nodes (Figure 3)
##### ============================================

# Simulation varying number of target nodes, one realization of a trivariate Matérn random field
# ==============================================================================================
root_nodes <- c(1000, 2000, 2500, 3000, 3500, 4000, 4500, 5000) # only number q (root of n=q^2 nodes)
timesimu <- c(0,length(root_nodes))
for(i in 1:length(root_nodes)){
  nx <- ny <-root_nodes[i]
  coord <- cbind(c(rep(1, times = ny))%x%c(1:nx), c(1:ny)%x%c(rep(1, times = nx))) # coord with qxq nodes
  temp <- system.time(specsim(coord, a0=100, nu0=0.5, model, parameters, N = 1, L = 1000, "temp_simu"))
  timesimu[i] <- as.numeric(temp[1])
}
# Plot of times vs nodes
#=======================
png(paste("Figure_3",".png",sep =''), width=1700, height=900, res = 180)
root_nodes2 <- c(1000, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
nodes <- root_nodes2^2 # n = q^2 nodes
plot(nodes, timesimu, col = "blue", type="b", xlab="Number of nodes", ylab="Time in seconds", xaxt='n', yaxt='n', lwd=1.5)
labelnodes <- c("1000*1000","2000*2000", "2500*2500", "3000*3000", "3500*3500", "4000*4000", 
                "4500*4500", "5000*5000")
axis(1, at=nodes, labels=labelnodes)
axis(2, at=timesimu, labels=floor(timesimu), las=2)
dev.off()


##### =====================================================================================
##### A particular case for SRF-2, IRF-0 and IRF-1 models, with L=100 and L=1000 (Figure 5)
##### =====================================================================================
# A vector containing the two values of the number of basic random fields
L <- c(100, 1000)
# Number of realizations for the evaluation of Shapiro-wilk test
N <- 10000
# Nine target locations
coord <- matrix(c(1, 1, 11, 1, 1, 11, 11, 11, 6, 1, 1, 6, 11, 6, 6, 11, 6, 6), nrow=9, ncol=2,
                byrow = TRUE)

# Simulation of Matérn bivariate model with L=100 and L=1000
# ==========================================================
# Building the parameters of the bivariate Matérn model
model <- "Matern"
a <- matrix(c(5, (0.5*(5^-2+10^-2))^-0.5, (0.5*(5^-2+10^-2))^-0.5, 10), 2, 2)
nu <- matrix(c(1.5, 1, 1, 0.5), 2, 2)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
# dimension
d <- 2
# Checking the validity of the built parameter matrices
validationMatern(a, nu, C, d)
# Building a list of parameters for the simulation
parameters.SRF2 <- list("C" = C, "a" = a, "nu1" = nu, "nu2" = c())
# Two simulations considering each value of L:
simu.SRF2.L100 <- specsim(coord, a0 = 5, nu0 = 0.5, model, parameters.SRF2, N, L[1], "SRF2_L100")
simu.SRF2.L1000 <- specsim(coord, a0 = 5, nu0 = 0.5, model, parameters.SRF2, N, L[2], "SRF2_L1000")
# Put the two simulations together in an array
simu.SRF2 <- array(c(simu.SRF2.L100, simu.SRF2.L1000), dim = c(9, 2*N, 2))

# Simulation of Power (IRF-0) bivariate model with L=100 and L=1000
# =================================================================
model <- "Power"
# Defining a list of parameters for the simulation
parameters.IRF0 <- list("C" = matrix(c(1, 0.5, 0.5, 1), 2, 2),"a" = matrix(c(10, 50, 50, 20), 2), 
                        "nu1" = matrix(c(1.5, 1, 1, 0.5),2), "nu2" = c())
# Two simulations considering each value of L:
simu.IRF0.L100 <- specsim(coord, a0 = 5, nu0 = 0.25, model, parameters.IRF0, N, L[1], "IRF0_L100")
simu.IRF0.L1000 <- specsim(coord, a0 = 5, nu0 = 0.25, model, parameters.IRF0, N, L[2], "IRF0_L1000")
# Put the two simulations together in one array
simuI.RF0 <- array(c(simu.IRF0.L100, simu.IRF0.L1000), dim = c(9, 2*N, 2))

# Simulation of Power (IRF-1) bivariate model with L=100 and L=1000
# =================================================================
model <- "Power"
# Defining a list of parameters for the simulation
parameters.IRF1 <- list("C" = matrix(c(1, 0.5, 0.5, 1), 2, 2),"a" = matrix(c(15, 15, 15, 5),2), 
                        "nu1" = matrix(c(2.5, 3, 3, 3.5),2), "nu2" = c())
# Two simulations considering each value of L:
simu.IRF1.L100 <- specsim(coord, a0 = 5, nu0 = 1.25, model, parameters.IRF1, N, L[1], "IRF1_L100")
simu.IRF1.L1000 <- specsim(coord, a0 = 5, nu0 = 1.25, model, parameters.IRF1, N, L[2], "IRF1_L1000")
# Put the two simulations together in one array
simu.IRF1 <- array(c(simu.IRF1.L100, simu.IRF1.L1000), dim = c(9, 2*N, 2))

# Preparing the Shapiro-Wilk test and plots
# =========================================
coord1 <- coord[c(1,2,3,4,9), ]         # subset x^(1)
coord2 <- coord[5:9, ]                  # subset x^(2)
lambda <- c(0.25, 0.25, 0.25, 0.25, -1) # set of weight

# P-P plots with Shapiro-Wilk test of the SRF-2, IRF-0 and IRF-1 realizations
# ===========================================================================
test(simu.SRF2, coord1, coord2, lambda, model="Matern", parameters.SRF2, N, "test_SRF2")
test(simu.IRF0, coord1, coord2, lambda, model="Power", parameters.IRF0, N, "test_IRF0")
test(simu.IRF1, coord1, coord2, lambda, model="Power", parameters.IRF1, N, "test_IRF1")
