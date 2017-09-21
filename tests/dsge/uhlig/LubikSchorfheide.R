
# Solve the Lubik-Schorfheide (2007) Model using Uhlig's method

rm(list=ls())
library(BMR)

#
# Setting parameter values as per Lubik and Schorfheide (2007, Table 3, posterior means)

psi1     <- 1.30
psi2     <- 0.23
psi3     <- 0.14

rhoR     <- 0.69
alpha    <- 0.11
rSS      <- 2.52
kappa    <- 0.32
tau      <- 0.31

rhoq     <- 0.31
rhoz     <- 0.42
rhoys    <- 0.97
rhopis   <- 0.46

sigma_r   <- 0.36
sigma_q   <- 1.25
sigma_z   <- 0.84
sigma_ys  <- 1.29
sigma_pis <- 2.00

beta      <- exp(-rSS/400)
tau.alpha <- tau + alpha*(2 - alpha)*(1 - tau)

#

A <- matrix(0,nrow=0,ncol=0)
B <- matrix(0,nrow=0,ncol=0)
C <- matrix(0,nrow=0,ncol=0)
D <- matrix(0,nrow=0,ncol=0)

#

F12 <- tau.alpha
F22 <- beta

#                   y,      pi,         R,         e,        x
F <- rbind(c(       1,     F12,         0,         0,        0),
           c(       0,     F22,         0,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0))

G13 <- -tau.alpha
G21 <- kappa/tau.alpha; G25 <- -kappa/tau.alpha
G31 <- (1-rhoR)*psi2; G32 <- (1-rhoR)*psi1; G34 <- (1-rhoR)*psi3;

G <- rbind(c(      -1,       0,       G13,         0,        0),
           c(     G21,      -1,         0,         0,      G25),
           c(     G31,     G32,        -1,       G34,        0),
           c(       0,       1,         0,        -1,        0),
           c(       0,       0,         0,         0,       -1))

H33 <- rhoR

H <- rbind(c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,       H33,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0))

#

J <- matrix(0,nrow=0,ncol=0)
K <- matrix(0,nrow=0,ncol=0)

#                   q,      ys,       pis,         z,        r
L <- rbind(c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        0))

M11 <- -alpha*tau.alpha*rhoq; M12 <- alpha*(1-alpha)*((1-tau)/tau)*rhoys; M14 <- -rhoz
M21 <- alpha*(beta*rhoq - 1)
M41 <- -(1-alpha)
M52 <- -alpha*(2-alpha)*((1-tau)/tau)

M <- rbind(c(     M11,     M12,         0,       M14,        0),
           c(     M21,       0,         0,         0,        0),
           c(       0,       0,         0,         0,        1),
           c(     M41,       0,        -1,         0,        0),
           c(       0,     M52,         0,         0,        0))

#

N <- rbind(c(    rhoq,       0,         0,         0,        0),
           c(       0,   rhoys,         0,         0,        0),
           c(       0,       0,    rhopis,         0,        0),
           c(       0,       0,         0,      rhoz,        0),
           c(       0,       0,         0,         0,        0))

Sigma <- rbind(c(  sigma_q^2,          0,            0,            0,           0),
               c(          0, sigma_ys^2,            0,            0,           0),
               c(          0,          0,  sigma_pis^2,            0,           0),
               c(          0,          0,            0,    sigma_z^2,           0),
               c(          0,          0,            0,            0,   sigma_r^2))

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

var_names <- c("Output","Inflation","IntRate","ExchangeRate","PotentialY",
               "ToT","OutputF","InflationF","Technology","IntRateShock")

dsge_irf <- IRF(dsge_obj,30,var_names=var_names)
dsge_sim <- dsge_obj$simulate(800,200)$sim_vals

#

y  <- dsge_sim[,1] + dsge_sim[,9]
pi <- dsge_sim[,2]
R  <- dsge_sim[,3]
e  <- dsge_sim[,4]
q  <- dsge_sim[,6]

LSData <- cbind(y,pi,R,e,q)
#save(LSData,file="LSData.RData")
