
# Solve the Lubik-Schorfheide (2007) Model using gensys

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

G0_13 <- tau.alpha
G0_16 <- alpha*tau.alpha*rhoq
G0_17 <- - alpha*(1-alpha)*((1-tau)/tau)*rhoys
G0_111 <- - tau.alpha

G0_21 <- - kappa/tau.alpha
G0_25 <- kappa/tau.alpha
G0_26 <- - alpha*(beta*rhoq - 1)

G0_31 <- - (1-rhoR)*psi2
G0_32 <- - (1-rhoR)*psi1
G0_34 <- - (1-rhoR)*psi3

G0_57 <- alpha*(2-alpha)*((1-tau)/tau)

#                        y,      pi,       R,       e,       x,       q,     y^s,    pi^s,       z,       r,   y_t+1,  pi_t+1,
Gamma0 <- rbind(c(       1,       0,   G0_13,       0,       0,   G0_16,   G0_17,       0,    rhoz,       0,      -1,  G0_111),
                c(   G0_21,       1,       0,       0,   G0_25,   G0_26,       0,       0,       0,       0,       0,   -beta),
                c(   G0_31,   G0_32,       1,   G0_34,       0,       0,       0,       0,       0,      -1,       0,       0),
                c(       0,      -1,       0,       1,       0, 1-alpha,       0,       1,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       1,       0,   G0_57,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0,       0),
                c(       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       1,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0))

Gamma1 <- rbind(c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,    rhoR,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,    rhoq,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,   rhoys,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,  rhopis,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,    rhoz,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       0,       1))

#

C <- matrix(0,ncol(Gamma0),1)

Psi <- matrix(0,ncol(Gamma0),5)
Psi[6:10,] <- diag(5)
     
Pi <- matrix(0,ncol(Gamma0),2)
Pi[11,1] <- 1
Pi[12,2] <- 1

Sigma <- rbind(c(  sigma_q^2,          0,            0,            0,           0),
               c(          0, sigma_ys^2,            0,            0,           0),
               c(          0,          0,  sigma_pis^2,            0,           0),
               c(          0,          0,            0,    sigma_z^2,           0),
               c(          0,          0,            0,            0,   sigma_r^2))

#

dsge_obj <- new(gensys)

dsge_obj$build(Gamma0,Gamma1,C,Psi,Pi)
dsge_obj$solve()

dsge_obj$G_sol
dsge_obj$impact_sol

#

dsge_obj$shocks_cov = Sigma

dsge_obj$state_space()

dsge_obj$F_state
dsge_obj$G_state

#

var_names <- c("Output","Inflation","IntRate","ExchangeRate","PotentialY",
               "ToT","OutputF","InflationF","Technology","IntRateShock")

dsge_irf <- IRF(dsge_obj,30,var_names=var_names)
dsge_sim <- dsge_obj$simulate(1000,200)$sim_vals

#

y  <- dsge_sim[,1] + dsge_sim[,9]
pi <- dsge_sim[,2]
R  <- dsge_sim[,3]
e  <- dsge_sim[,4]
q  <- dsge_sim[,6]

LSData <- cbind(y,pi,R,e,q)
#save(LSData,file="LSData.RData")
