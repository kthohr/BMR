
#
# Simple NK Model

nkm_model_fn <- function(parameters)
{
    alpha    <- parameters[1]
    beta     <- parameters[2]
    vartheta <- parameters[3]
    theta    <- parameters[4]
    
    eta    <- parameters[5]               
    phi    <- parameters[6]
    phi_pi <- parameters[7]
    phi_y  <- parameters[8]
    rho_a  <- parameters[9]
    rho_v  <- parameters[10]
    
    sigma_T <- parameters[11]
    sigma_M <- parameters[12]
    
    BigTheta <- (1-alpha)/(1-alpha+alpha*vartheta)
    kappa    <- (((1-theta)*(1-beta*theta))/(theta))*BigTheta*((1/eta)+((phi+alpha)/(1-alpha)))
    psi      <- (eta*(1+phi))/(1-alpha+eta*(phi + alpha))
    
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

nkm_model_simple <- function(parameters)
{
    eta <- parameters[1]
    
    #
    
    pars_inp <- numeric(12)
    
    pars_inp[1]  <- 0.33   # alpha
    pars_inp[2]  <- 0.99   # beta
    pars_inp[3]  <- 6.0    # vartheta
    pars_inp[4]  <- 0.6667 # theta
    
    pars_inp[5]  <- eta
    pars_inp[6]  <- 1      # phi
    pars_inp[7]  <- 1.5    # phi_pi
    pars_inp[8]  <- 0.125  # phi_y
    pars_inp[9]  <- 0.90   # rho_a
    pars_inp[10] <- 0.5    # rho_v
    
    pars_inp[11] <- 1      # sigma_T
    pars_inp[12] <- 0.25   # sigma_M
    
    return(nkm_model_fn(pars_inp))
}
