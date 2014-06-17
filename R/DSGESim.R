DSGESim.SDSGE <- function(obj,seedval=1122,shocks,sim.periods,burnin=NULL,hpfiltered=FALSE,lambda=1600){
  results <- .dsgesimulation(obj,seedval,shocks,sim.periods,burnin,hpfiltered,lambda)
  return(results)
}

.dsgesimulation <- function(obj,seedval,shocks,sim.periods,burnin=NULL,hpfiltered=FALSE,lambda=1600){
  StateMats <- .DSGEstatespace(obj$N,obj$P,obj$Q,obj$R,obj$S)
  F <- StateMats$F; G <- StateMats$G
  #
  if(class(burnin) != "numeric"){
    burnin <- ceiling(0.5*sim.periods)
  }
  #
  DSim <- matrix(0,nrow=(sim.periods+burnin),ncol=ncol(F))
  #
  set.seed(seedval)
  #
  shockvec1 <- matrix(rep(0,nrow(obj$N)*(nrow(G) - nrow(obj$N))),ncol=nrow(obj$N))
  shockvec2 <- diag(nrow(obj$N))
  for(j in 1:nrow(obj$N)){
    shockvec2[j,j] <- rnorm(1,mean=0,sd=shocks[j])
  }
  shockvec <- rbind(shockvec1,shockvec2)
  #shockvec <- rbind(shockvec1,c(rep(rnorm(1),nrow(obj$N))))
  if(nrow(obj$N) > 1){
    DSim[1,] <- apply(t(G%*%shockvec),c(2),sum)
    for(i in 2:burnin){
      for(j in 1:nrow(obj$N)){
        shockvec2[j,j] <- rnorm(1,mean=0,sd=shocks[j])
      }
      shockvec <- rbind(shockvec1,shockvec2)
      shockvec3 <- apply(t(G%*%shockvec),c(2),sum)
      DSim[i,] <- DSim[i-1,]%*%t(F) + shockvec3
    }
    #
    for(i in 1:sim.periods){
      for(j in 1:nrow(obj$N)){
        shockvec2[j,j] <- rnorm(1,mean=0,sd=shocks[j])
      }
      shockvec <- rbind(shockvec1,shockvec2)
      shockvec3 <- apply(t(G%*%shockvec),c(2),sum)
      DSim[i+burnin,] <- DSim[i-1+burnin,]%*%t(F) + shockvec3
    }
    #
    DSim <- DSim[(burnin+1):(burnin+sim.periods),]
  }else{
    DSim[1,] <- t(G%*%shockvec)
    #
    for(i in 2:burnin){
      shockvec <- rbind(shockvec1,c(rnorm(1,0,shocks[1])))
      DSim[i,] <- DSim[i-1,]%*%t(F) + t(G%*%shockvec)
    }
    #
    for(i in 1:sim.periods){
      shockvec <- rbind(shockvec1,c(rnorm(1,0,shocks[1])))
      DSim[i+burnin,] <- DSim[i-1+burnin,]%*%t(F) + t(G%*%shockvec)
    }
    #
    DSim <- DSim[(burnin+1):(burnin+sim.periods),]
  }
  #
  if(hpfiltered==TRUE){
    hpfilterq <- function(x,lambda=1600){
      eye <- diag(length(x))
      result <- solve(eye+lambda*crossprod(diff(eye,lag=1,d=2)),x)
      return(result)
    }
    for(j in 1:ncol(DSim)){
      DSim[,j]<-hpfilterq(DSim[,j])
    }
  }
  #
  return(DSim)
}