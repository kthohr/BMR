
# Solve the An-Schorfheide (2007) Model using gensys

rm(list=ls())
library(BMR)

#
# Setting parameter values as per An and Schorfheide (2007, Table 2)

tau      <- 2.00
kappa    <- 0.15

psi1     <- 1.50
psi2     <- 1.00

rhoR     <- 0.60
rhog     <- 0.95
rhoz     <- 0.65

rA       <- 0.40
piA      <- 4.00
gammaQ   <- 0.50

sigma_r  <- sqrt(0.002)
sigma_g  <- sqrt(0.008)
sigma_z  <- sqrt(0.0045)

gamma    <- 1 + (gammaQ/100)
beta     <- 1/( 1 + (rA/400) )
pi       <- 1 + (piA/400)

#

G0_16 <- rhoz/tau
G0_42 <- (1 - rhoR)*psi1
G0_43 <- (1 - rhoR)*psi2

#                        c,      pi,       y,       R,       r,       g,       z    c_t+1,  pi_t+1,
Gamma0 <- rbind(c(      -1,       0,       0,  -1/tau,       0,       0,   G0_16,       1,   1/tau),
                c(       0,      -1,   kappa,       0,       0,  -kappa,       0,       0,    beta),
                c(       1,       0,      -1,       0,       0,       1,       0,       0,       0),
                c(       0,  -G0_42,  -G0_43,       1,      -1,   G0_43,       0,       0,       0),
                c(       0,       0,       0,       0,       1,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       1,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       1,       0,       0),
                c(       1,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       1,       0,       0,       0,       0,       0,       0,       0))

Gamma1 <- rbind(c(       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,    rhoR,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,    rhog,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,    rhoz,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       1,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       1))

#

C <- matrix(0,9,1)

Psi <- matrix(0,9,3)
Psi[5,1] <- 1
Psi[6,2] <- 1
Psi[7,3] <- 1
     
Pi <- matrix(0,9,2)
Pi[8,1] <- 1
Pi[9,2] <- 1

Sigma <- rbind(c( sigma_r^2,         0,         0),
               c(         0, sigma_g^2,         0),
               c(         0,         0, sigma_z^2))

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

var_names <- c("Consumption","Inflation","Output","IntRate","IntRateS","GovSpending","Technology")

dsge_irf <- IRF(dsge_obj,30,var_names=var_names)
dsge_sim <- dsge_obj$simulate(1000,200)$sim_vals

#

y  <- gammaQ + 100*(dsge_sim[,3] + dsge_sim[,7])
pi <- piA + 400*dsge_sim[,2]
R  <- piA + rA + 4*gammaQ + 400*dsge_sim[,4]

ASData <- cbind(y,pi,R)
#save(ASData,file="ASData.RData")
