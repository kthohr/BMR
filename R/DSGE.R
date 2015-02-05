SDSGE.default <-function(A,B,C,D,F,G,H,J,K,L,M,N,whichEig=NULL){
  #
  l <- as.numeric(dim(C)[1])
  n <- as.numeric(dim(C)[2])
  k <- as.numeric(min(dim(N)))
  #
  m <- 0
  #
  if (l == 0) {
    m <- as.numeric(dim(F)[2])
  }else if(l >0){
    m <- as.numeric(dim(A)[2])
  }
  #
  pickeig <- 0
  if(class(whichEig) !="NULL"){
    pickeig <- 1
    whichEig <- c(whichEig)
  }else{
    whichEig <- 1
  }
  #
  UhligSol <- .Call("UhligCpp", A,B,C,D,F,G,H,J,K,L,M,N,pickeig,matrix(whichEig,ncol=1),l,n,k,m, PACKAGE = "BMR", DUP = FALSE)
  #
  sdsgeret <- list(N=N,P=UhligSol$P,Q=UhligSol$Q,R=UhligSol$R,S=UhligSol$S,EigenValues=UhligSol$EigValueSorted,EigenVectors=UhligSol$EigVecSorted)
  class(sdsgeret) <- "SDSGE"
  return(sdsgeret)
}

.DSGEstatespace <- function(N,P,Q,R,S){
  KA<-rbind(P,R)
  KB<-rbind(Q,S)
  #
  # The F matrix is (KA,KB*N;0,N)
  F0 <- matrix(0,nrow=nrow(KA),ncol=nrow(R))
  #
  F1<-cbind(KA,F0,KB%*%N)
  F2<-cbind(matrix(rep(0,nrow(KA)*nrow(N)),nrow=ncol(N)),N)
  F<-rbind(F1,F2)
  #
  # G is the shock matrix for the state matrices
  G1<-rbind(KB,diag(ncol(N)))
  G2<-diag(nrow(G1))*0
  G<-cbind(G2[,(ncol(N)+1):ncol(G2)],G1)
  #
  return=list(F=F,G=G)
}