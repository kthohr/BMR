
#
# Solve a basic RBC model with BMR using Uhlig's method

rm(list=ls())
library(BMR)

#

beta      <- 0.99
alpha     <- .33
delta     <- .015
eta       <- 1.0
rho       <- 0.95
sigma_T   <- 1

RSS  <- 1/beta
YKSS <- (RSS + delta - 1)/alpha
IKSS <- delta
IYSS <- ((alpha*delta)/(RSS + delta - 1))
CYSS <- 1 - IYSS

#

A <- matrix(c( 0,
               0,
               1,
               0,
               0 ))


B <- matrix(c(                0,
                         -alpha,
                     -(1-delta),
                              0,
               (alpha/RSS)*YKSS))

#              consumption,              output,          labor,    interest,   investment
C <- rbind(c(        -CYSS,                   1,              0,           0,      -IYSS), 
           c(            0,                   1,     -(1-alpha),           0,          0),
           c(            0,                   0,              0,           0,      -IKSS),      
           c(          eta,                  -1,              1,           0,          0),
           c(            0,   -(alpha/RSS)*YKSS,              0,           1,          0))

D <- matrix(c(  0,
               -1,
                0,
                0,
                0))

F <-  matrix(c(0))
G <- matrix(c(0))
H <- matrix(c(0))

J <- t(matrix(c( -1,  0,  0,  (1/eta),  0)))
K <- t(matrix(c( 1,   0,  0,  0,  0)))
L <- matrix(c(0)) 
M <- matrix(c(0)) 
N <- matrix(rho)

#

Sigma <- matrix(sigma_T^2,1,1)

#

dsge_obj <- new(uhlig)

dsge_obj$build(A,B,C,D,F,G,H,J,K,L,M,N)
dsge_obj$solve()

dsge_obj$P_sol
dsge_obj$Q_sol
dsge_obj$R_sol
dsge_obj$S_sol

#

dsge_obj$shocks_cov = Sigma

dsge_obj$state_space()

dsge_obj$F_state
dsge_obj$G_state

#

var_names <- c("Capital","Consumption","Output","Labour","Interest","Investment","Technology")

dsge_irf <- IRF(dsge_obj,30,var_names=var_names)
dsge_sim <- dsge_obj$simulate(800,200)

#
#END
