
#

rm(list=ls())
library(BMR)

#

data(BMRVARData)
dsgedata <- USMacroData[24:211,-c(1,3)]
dsgedata <- as.matrix(dsgedata)
for(i in 1:2){
    dsgedata[,i] <- dsgedata[,i] - mean(dsgedata[,i])
}

#

nkm_model_fn <- function(parameters){
    alpha    <- 0.33
    vartheta <- 6
    beta     <- 0.99
    theta    <- 0.6667
    
    eta    <- parameters[1]               
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
    G0_47 <- (1/eta)*psi*(rho_a - 1)
    #Order:                 yg,       y,      pi,      rn,       i,       n,       a,       v,  yg_t+1,  pi_t+1
    Gamma0 <- rbind(c(      -1,       0,       0,     eta,  -eta/4,       0,       0,       0,       1,   eta/4),
                    c(   kappa,       0,    -1/4,       0,       0,       0,       0,       0,       0,  beta/4),
                    c(   phi_y,       0,phi_pi/4,       0,    -1/4,       0,       0,       1,       0,       0),
                    c(       0,       0,       0,      -1,       0,       0,   G0_47,       0,       0,       0),
                    c(       0,      -1,       0,       0,       0, 1-alpha,       1,       0,       0,       0),
                    c(      -1,       1,       0,       0,       0,       0,    -psi,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       1,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       1,       0,       0),
                    c(       1,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       1,       0,       0,       0,       0,       0,       0,       0))
    #
    Gamma1 <- rbind(c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,   rho_a,       0,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,   rho_v,       0,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       1,       0),
                    c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       1))
    #
    C <- matrix(0,10,1)
    #
    Psi <- matrix(0,10,2)
    Psi[7,1] <- 1
    Psi[8,2] <- 1
    #     
    Pi <- matrix(0,10,2)
    Pi[9,1] <- 1
    Pi[10,2] <- 1
    #
    
    shocks <- matrix(c(sigma_T^2,      0,  
                       0, sigma_M^2),nrow=2)
    
    #
    
    ObserveMat <- cbind(c(0,0,1,0,0,0,0,0),
                        c(0,0,0,0,1,0,0,0))
    
    ObsCons  <- matrix(0,nrow=2,ncol=1)
    MeasErrs <- matrix(0,nrow=2,ncol=2)
    #
    return=list(Gamma0=Gamma0,Gamma1=Gamma1,GammaC=C,Psi=Psi,Pi=Pi,shocks_cov_out=shocks,C_out=ObsCons,H_out=ObserveMat,R_out=MeasErrs)
}

obj <- new(dsgevar_gensys)
obj$set_model_fn(nkm_model_fn)

x <- c(1)

obj$eval_model(x)

lrem_obj = obj$get_lrem()
lrem_obj$solve()

lrem_obj$shocks_cov <- matrix(c(1,0,0,0.125),2,2,byrow=TRUE)

sim_data <- lrem_obj$simulate(200,800)$sim_vals

sim_data <- cbind(sim_data[,3],sim_data[,5])

#

prior_pars <- cbind(c(1.0),
                    c(0.05))

prior_form <- c(1)

obj$set_prior(prior_form,prior_pars)

#

par_bounds <- cbind(c(-Inf),
                    c( Inf))

opt_bounds <- cbind(c(0.5),
                    c(3.0))

obj$set_bounds(opt_bounds[,1],opt_bounds[,2])

obj$opt_initial_lb <- opt_bounds[,1]
obj$opt_initial_ub <- opt_bounds[,2]

#

obj$build(sim_data,FALSE,1,1.0);

obj$estim_mode(x)
