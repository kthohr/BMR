/*################################################################################
  ##
  ##   R package BMR by Keith O'Hara Copyright (C) 2011, 2012, 2013, 2014
  ##   This file is part of the R package BMR.
  ##
  ##   The R package BMR is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 3 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package BMR is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

#include "dsgec.h"
using namespace Rcpp;

SEXP UhligCpp( SEXP mA, SEXP mB, SEXP mC, SEXP mD, SEXP mF, SEXP mG, SEXP mH,
               SEXP mJ, SEXP mK, SEXP mL, SEXP mM, SEXP mN, 
               SEXP mpickEig, SEXP mwhichEig, SEXP ml, SEXP mn, 
               SEXP mk, SEXP mm )
{
  //Suppress warnings:
  std::ostream nullstream(0);
  arma::set_stream_err2(nullstream);
  //
  arma::mat A = as<arma::mat>(mA);
  arma::mat B = as<arma::mat>(mB);
  arma::mat C = as<arma::mat>(mC);
  arma::mat D = as<arma::mat>(mD);
  arma::mat F = as<arma::mat>(mF);
  arma::mat G = as<arma::mat>(mG);
  arma::mat H = as<arma::mat>(mH);
  arma::mat J = as<arma::mat>(mJ);
  arma::mat K = as<arma::mat>(mK);
  arma::mat L = as<arma::mat>(mL);
  arma::mat M = as<arma::mat>(mM);
  arma::mat N = as<arma::mat>(mN);
  //
  int l = as<int>(ml);
  int n = as<int>(mn);
  int k = as<int>(mk);
  int m = as<int>(mm);
  //
  int i;
  //
  arma::mat Xi(2*m,2*m);
  Xi.zeros();
  arma::mat Delta(2*m,2*m);
  Delta.zeros();
  //
  arma::mat Psi = F;
  arma::mat Gamma = -G;
  arma::mat Theta = -H;
  //
  arma::mat Cpinv = C;
  if(l > 0){
    Cpinv = pinv(C);
    Psi = F - J*Cpinv*A;
    Gamma = J*Cpinv*B - G + K*Cpinv*A;
    Theta = K*Cpinv*B - H;
  }
  //
  Xi(arma::span(0,m-1),arma::span(0,m-1)) = Gamma;
  Xi(arma::span(0,m-1),arma::span(m,2*m-1)) = Theta;
  Xi(arma::span(m,2*m-1),arma::span(0,m-1)) = arma::eye(m,m);
  Xi(arma::span(m,2*m-1),arma::span(m,2*m-1)).zeros();
  //
  Delta(arma::span(0,m-1),arma::span(0,m-1)) = Psi;
  Delta(arma::span(0,m-1),arma::span(m,2*m-1)).zeros();
  Delta(arma::span(m,2*m-1),arma::span(0,m-1)).zeros();
  Delta(arma::span(m,2*m-1),arma::span(m,2*m-1)) = arma::eye(m,m);
  //
  //
  //
  arma::cx_vec EigValue;
  arma::cx_mat EigVec;
  // QZ:
  arma::eig_pair(EigValue, EigVec, Xi, Delta);
  //
  int pickEig = as<int>(mpickEig);
  arma::vec whichEig = as<arma::vec>(mwhichEig);
  //
  arma::vec EigValueAbs = abs(EigValue);
  arma::vec EigValueReal = real(EigValue);
  // Deal with infinite values (otherwise sort will return an error):
  float bignum = 10000000;
  //
  arma::uvec infindices = find_nonfinite(EigValueAbs);
  arma::uvec neginfind = arma::find(EigValueReal < -bignum);
  //
  float infin2 = infindices.n_elem;
  float infin3 = neginfind.n_elem;
  //
  if(infin2 > 0){
    arma::cx_vec BigNum(infin2);
    BigNum.ones();
    EigValueAbs.elem( infindices ) = arma::ones(infin2)*bignum;
    EigValue.elem( infindices ) = BigNum*bignum;
  }
  if(infin3 > 0){
    arma::cx_vec BigNum2(infin3);
    BigNum2.ones();
    EigValue.elem( neginfind ) = -BigNum2*bignum;
  }
  // Now sort
  arma::vec EigValueAbsSorted = sort(EigValueAbs);
  arma::cx_vec EigValueSorted = sort(EigValue);
  arma::uvec indices = arma::sort_index(EigValueAbs);
  //
  //arma::fmat egvecinds = arma::conv_to<arma::fmat>::from(indices);
  //
  float egvecind = 1;
  arma::cx_mat EigVecSorted(2*m,2*m);
  EigVecSorted.zeros();
  for(i = 1;i<=(2*m);i++) {
    egvecind = arma::as_scalar(indices(arma::span(i-1,i-1),arma::span(0,0)));
    EigVecSorted(arma::span(),arma::span(i-1,i-1)) = EigVec(arma::span(),arma::span(egvecind,egvecind));
  }
  //
  arma::cx_mat EigVecSortedRet = EigVecSorted;
  arma::cx_vec EigValueSortedRet = EigValueSorted;
  //
  if(pickEig > 0){
    for(i = 1;i<=m;i++) {
      egvecind = arma::as_scalar(whichEig(arma::span(i-1,i-1),arma::span(0,0)));
      EigVecSorted(arma::span(),arma::span(i-1,i-1)) = EigVecSorted(arma::span(),arma::span(egvecind-1,egvecind-1));
      EigValueSorted(arma::span(i-1,i-1),arma::span()) = EigValueSorted(arma::span(egvecind-1,egvecind-1),arma::span());
    }
  }
  //
  arma::cx_vec LambdaVec = arma::conv_to<arma::cx_vec>::from(EigValueSorted(arma::span(0,m-1),arma::span()));
  arma::cx_mat Lambda = arma::diagmat(LambdaVec);
  //
  arma::cx_mat Omega = EigVecSorted(arma::span(m,2*m-1),arma::span(0,m-1));
  //
  arma::cx_mat Pi = Omega*Lambda*pinv(Omega);
  arma::mat P = arma::real(Pi);
  //
  //
  //
  Cpinv = C;
  //
  arma::mat LNM = L*N + M;
  arma::colvec LNMStacked(LNM.begin(), LNM.size(), false);
  //
  arma::mat V = kron(trans(N),F) + kron(arma::eye(k,k),(F*P+G));
  //V.zeros();
  arma::mat V2((k*A.n_rows) + V.n_rows, (k*A.n_rows) + V.n_rows);
  V2.zeros();
  //
  arma::vec QS(1,1);
  QS.zeros();
  arma::mat R(0,P.n_cols);
  R.zeros();
  if(l == 0){
    QS = -inv(V)*LNMStacked;
  }else{
    Cpinv = pinv(C);
    R = - Cpinv*(A*P+B);
    //
    V2(arma::span(0,(k*A.n_rows)-1),arma::span(0,(k*A.n_cols)-1)) = kron(arma::eye(k,k),A);
    V2(arma::span(0,(k*A.n_rows)-1),arma::span(k*A.n_cols,(k*A.n_cols)+(k*C.n_cols)-1)) = kron(arma::eye(k,k),C);
    V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(0,(k*F.n_cols)-1)) = kron(trans(N),F) + kron(arma::eye(k,k),(F*P + J*R + G));
    V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(k*F.n_cols,(k*F.n_cols)+(k*J.n_cols)-1)) = kron(trans(N),J) + kron(arma::eye(k,k),K);
    //
    arma::colvec DStacked(D.begin(), D.size(), false);
    arma::mat DLNM(DStacked.n_rows + LNMStacked.n_rows,1);
    DLNM.zeros();
    DLNM(arma::span(0,DStacked.n_rows-1),arma::span()) = DStacked;
    DLNM(arma::span(DStacked.n_rows,DStacked.n_rows + LNMStacked.n_rows-1),arma::span()) = LNMStacked;
    QS = -inv(V2)*DLNM;
  }
  arma::vec QVec = QS(arma::span(0,(m*k)-1),arma::span());
  arma::mat Q(QVec.begin(),m,k);
  //
  arma::mat S(0,Q.n_cols);
  S.zeros();
  if(l > 0){
    arma::vec Sc = QS(arma::span(m*k,(m+n)*k-1),arma::span());
    arma::mat S2(Sc.begin(),n,k);
    S = S2;
  }
  //
  return Rcpp::List::create(Rcpp::Named("P") = P,Rcpp::Named("Q") = Q,Rcpp::Named("R") = R,Rcpp::Named("S") = S,Rcpp::Named("EigValueSorted") = EigValueSortedRet,Rcpp::Named("EigVecSorted") = EigVecSortedRet);
}

SEXP DSGEKalman( SEXP mdsgedata, SEXP mObserveMat, SEXP mObsCons, SEXP mF, SEXP mG, SEXP mN, SEXP mshocks, 
                 SEXP mR, SEXP mMaxIter )
{
  arma::mat dsgedata = as<arma::mat>(mdsgedata);
  arma::mat ObserveMat = as<arma::mat>(mObserveMat);
  arma::mat ObsCons = as<arma::mat>(mObsCons);
  arma::mat F = as<arma::mat>(mF);
  arma::mat G = as<arma::mat>(mG);
  arma::mat N = as<arma::mat>(mN);
  arma::mat shocks = as<arma::mat>(mshocks);
  arma::mat R = as<arma::mat>(mR);
  //
  double MaxIter = as<int>(mMaxIter);
  //
  arma::mat LogLikelihood(1,1);
  LogLikelihood.zeros();
  //
  int T = dsgedata.n_rows;
  double ndata = dsgedata.n_cols;
  double nlog2pi = ndata*log(2*arma::datum::pi);
  //
  int i;
  //
  arma::mat tObserveMat = trans(ObserveMat);
  arma::mat tF = trans(F);
  //
  arma::mat Q(G.n_rows,G.n_rows);
  Q.zeros();
  Q(arma::span(G.n_rows-N.n_rows,G.n_rows-1),arma::span(G.n_rows-N.n_rows,G.n_rows-1)) = shocks;
  //
  arma::mat GQG = G*Q*trans(G);
  //
  //
  //
  arma::mat InitialState(F.n_cols,1);
  InitialState.zeros();
  //
  arma::mat InitialCov = arma::eye(F.n_cols,F.n_cols);
  //
  arma::mat SSCovOld = GQG;
  arma::mat SSCovNew = GQG;
  arma::mat IMat = F;
  //
  for(i=1;i<=MaxIter;i++){
    SSCovNew = SSCovOld + IMat*SSCovOld*trans(IMat);
    IMat *= IMat;
    //
    SSCovOld = SSCovNew;
  }
  InitialCov = SSCovNew;
  //
  arma::mat StatePredicted = F*InitialState;
  arma::mat StateCovPredicted = F*InitialCov*tF + GQG;
  arma::mat Sigma = tObserveMat*StateCovPredicted*ObserveMat + R;
  arma::mat invSigma = inv_sympd(Sigma);
  //
  arma::mat KalmanGain = StateCovPredicted*ObserveMat*invSigma;
  arma::mat KalmanResid = trans(dsgedata(arma::span(0,0),arma::span())) - tObserveMat*StatePredicted - ObsCons;
  arma::mat StateFiltered = StatePredicted + KalmanGain*KalmanResid;
  arma::mat StateCovFiltered = StateCovPredicted - KalmanGain*tObserveMat*StateCovPredicted;
  LogLikelihood += nlog2pi + log(arma::det(Sigma)) + trans(KalmanResid)*invSigma*KalmanResid;
  #
  for(i=2; i<=T; i++){
    StatePredicted = F*StateFiltered;
    StateCovPredicted = F*StateCovFiltered*tF + GQG;
    Sigma = tObserveMat*StateCovPredicted*ObserveMat + R;
    invSigma = inv_sympd(Sigma);
    //
    KalmanGain = StateCovPredicted*ObserveMat*invSigma;
    KalmanResid = trans(dsgedata(arma::span(i-1,i-1),arma::span())) - tObserveMat*StatePredicted - ObsCons;
    StateFiltered = StatePredicted + KalmanGain*KalmanResid;
    StateCovFiltered = StateCovPredicted - KalmanGain*tObserveMat*StateCovPredicted;
    LogLikelihood += nlog2pi + log(arma::det(Sigma)) + trans(KalmanResid)*invSigma*KalmanResid;
  }
  // Returns the negative of the loglikelihood
  LogLikelihood = 0.5*arma::as_scalar(LogLikelihood);
  //
  return Rcpp::List::create(Rcpp::Named("dsgelike") = LogLikelihood);
}

SEXP DSGECR( SEXP mdsgedata, SEXP mObserveMat, SEXP mObsCons, SEXP mF, SEXP mG, SEXP mN, SEXP mshocks, 
             SEXP mR, SEXP mMaxIter )
{
  arma::mat dsgedata = as<arma::mat>(mdsgedata);
  arma::mat ObserveMat = as<arma::mat>(mObserveMat);
  arma::mat ObsCons = as<arma::mat>(mObsCons);
  arma::mat F = as<arma::mat>(mF);
  arma::mat G = as<arma::mat>(mG);
  arma::mat N = as<arma::mat>(mN);
  arma::mat shocks = as<arma::mat>(mshocks);
  arma::mat R = as<arma::mat>(mR);
  //
  double MaxIter = as<int>(mMaxIter);
  //
  arma::mat LogLikelihood(1,1);
  LogLikelihood.zeros();
  //
  int T = dsgedata.n_rows;
  double ndata = dsgedata.n_cols;
  double nlog2pi = ndata*log(2*arma::datum::pi);
  //
  int i;
  //
  arma::mat tObserveMat = trans(ObserveMat);
  arma::mat tF = trans(F);
  //
  arma::mat Q(G.n_rows,G.n_rows);
  Q.zeros();
  Q(arma::span(G.n_rows-N.n_rows,G.n_rows-1),arma::span(G.n_rows-N.n_rows,G.n_rows-1)) = shocks;
  //
  arma::mat GQG = G*Q*trans(G);
  //
  //
  //
  arma::mat InitialState(F.n_cols,1);
  InitialState.zeros();
  //
  arma::mat InitialCov = arma::eye(F.n_cols,F.n_cols);
  //
  arma::mat SSCovOld = GQG;
  arma::mat SSCovNew = GQG;
  arma::mat IMat = F;
  //
  for(i=1;i<=MaxIter;i++){
    SSCovNew = SSCovOld + IMat*SSCovOld*trans(IMat);
    IMat *= IMat;
    //
    SSCovOld = SSCovNew;
  }
  InitialCov = SSCovNew;
  //
  arma::mat StateFiltered = InitialState;
  arma::mat StateCovPredicted = F*InitialCov*tF + GQG;
  arma::mat Sigma = tObserveMat*StateCovPredicted*ObserveMat + R;
  arma::mat iSigma = inv_sympd(Sigma);
  //
  arma::mat St = F*StateCovPredicted*ObserveMat;
  arma::mat Mt = -iSigma;    
  arma::mat Kt = St*iSigma;
  //
  arma::mat Resid = trans(dsgedata(arma::span(0,0),arma::span())) - tObserveMat*StateFiltered - ObsCons;
  //
  LogLikelihood += nlog2pi + log(arma::det(Sigma)) + trans(Resid)*iSigma*Resid;
  //
  StateFiltered = F*StateFiltered + Kt*Resid;
  //
  arma::mat tObserveMatSt = tObserveMat*St;
  arma::mat MSpZp = Mt*trans(tObserveMatSt);
  arma::mat FSt = F*St;
  //
  arma::mat Sigma1  = Sigma + tObserveMatSt*MSpZp;     
  arma::mat iSigma1 = arma::inv_sympd(Sigma1);
  //
  Kt = (Kt*Sigma + FSt*MSpZp)*iSigma1;
  St = FSt - Kt*tObserveMatSt;
  Mt += MSpZp*iSigma*trans(MSpZp);
  Sigma = Sigma1;
  iSigma = iSigma1;
  //
  for(i=2; i<=T; i++){
    //
    Resid = trans(dsgedata(arma::span(i-1,i-1),arma::span())) - tObserveMat*StateFiltered - ObsCons;
    //
    LogLikelihood += nlog2pi + log(arma::det(Sigma)) + trans(Resid)*iSigma*Resid;
    //
    StateFiltered = F*StateFiltered + Kt*Resid;
    //
    tObserveMatSt = tObserveMat*St;
    MSpZp = Mt*trans(tObserveMatSt);
    FSt = F*St;
    //
    Sigma1 += tObserveMatSt*MSpZp;       
    iSigma1 = inv_sympd(Sigma1);
    //
    Kt = (Kt*Sigma + FSt*MSpZp)*iSigma1;
    St = FSt - Kt*tObserveMatSt;
    Mt += MSpZp*iSigma*trans(MSpZp);
    Sigma = Sigma1;
    iSigma = iSigma1;
  }
  // Returns the negative of the loglikelihood
  LogLikelihood = 0.5*arma::as_scalar(LogLikelihood);
  //
  return Rcpp::List::create(Rcpp::Named("dsgelike") = LogLikelihood);
}