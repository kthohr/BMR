
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

    ObserveMat <- cbind(c(0,0,1,0,0,0,0,0),
                        c(0,0,0,0,1,0,0,0))
    
    ObsCons  <- matrix(0,nrow=2,ncol=1)
    MeasErrs <- matrix(0,nrow=2,ncol=2)

    return=list(A=A,B=B,C=C,D=D,F=F,G=G,H=H,J=J,K=K,L=L,M=M,N=N,shocks_cov_out=Sigma,C_out=ObsCons,H_out=ObserveMat,R_out=MeasErrs)
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
