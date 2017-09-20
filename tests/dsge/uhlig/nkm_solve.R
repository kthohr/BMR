#
# Solve the NK model with Uhlig's method
#

rm(list=ls())
library(BMR)

#

alpha    <- 0.33
vartheta <- 6
beta     <- 0.99
theta    <- 0.6667

eta    <- 1               
phi    <- 1                  
phi_pi <- 1.5             
phi_y  <- 0.5/4
rho_a  <- 0.90
rho_v  <- 0.5

BigTheta <- (1-alpha)/(1-alpha+alpha*vartheta)
kappa    <- (((1-theta)*(1-beta*theta))/(theta))*BigTheta*((1/eta)+((phi+alpha)/(1-alpha)))
psi      <- (eta*(1+phi))/(1-alpha+eta*(phi + alpha))

sigma_T <- 1
sigma_M <- 0.25

#

A <- matrix(0,nrow=0,ncol=0)
B <- matrix(0,nrow=0,ncol=0)
C <- matrix(0,nrow=0,ncol=0)
D <- matrix(0,nrow=0,ncol=0)

#Order:            yg,       y,        pi,        rn,        i,           n
F <- rbind(c(      -1,       0, -0.25*eta,         0,        0,           0),
           c(       0,       0,-0.25*beta,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0))

G33 <- -0.25*phi_pi
G <- rbind(c(       1,       0,         0,    -1*eta, 0.25*eta,           0),
           c(  -kappa,       0,      0.25,         0,        0,           0),
           c(  -phi_y,       0,       G33,         0,     0.25,           0),
           c(       0,       0,         0,         1,        0,           0),
           c(       0,       1,         0,         0,        0,  -(1-alpha)),
           c(       1,      -1,         0,         0,        0,           0))

H <- rbind(c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0),
           c(       0,       0,         0,         0,        0,           0))

J <- matrix(0,nrow=0,ncol=0)
K <- matrix(0,nrow=0,ncol=0)
L <- matrix(0,nrow=6,ncol=2)

M41 <- -(1/eta)*psi*(rho_a - 1)
M <- rbind(c(   0,  0),
           c(   0,  0),
           c(   0, -1),
           c( M41,  0),
           c(  -1,  0),
           c( psi,  0))

#

N <- rbind(c( rho_a,     0),  
           c(     0, rho_v))

Sigma <- rbind(c(sigma_T^2,         0),  
               c(        0, sigma_M^2))

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

var_names <- c("Output Gap","Output","Inflation","Natural Int",
               "Nominal Int","Labour Supply","Technology","MonetaryPolicy")

dsge_irf <- IRF(dsge_obj,12,var_names=var_names)
dsge_sim <- dsge_obj$simulate(800,200)

#
#END
