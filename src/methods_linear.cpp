#include <RcppArmadillo.h>
#include "methods_linear.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 1. PCA : Principal Component Analysis
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_pca(arma::mat& psdX){
  // define auxiliary stuffs
  arma::vec eigval;
  arma::mat eigvec;
  // run eigendecomposition
  eig_sym(eigval, eigvec, psdX);
  // return output
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}

// 2. MDS : Multi-dimensional Scaling (Classical)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_mds(arma::mat& centerX){
  const int n = centerX.n_cols;
  // X(centered) = U*S*V.t()
  mat U;
  vec s;
  mat V;
  svd(U, s, V, centerX);
  // output
  return Rcpp::List::create(Rcpp::Named("eigval")=s,
                            Rcpp::Named("eigvec")=V);
}

// 3. MDS given D
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_mdsD(arma::mat& D){
  const int n = D.n_cols;
  // 3-1. H = I-(1/n)*ones(n,n)
  mat D2 = pow(D,2);
  mat I(n,n); I.eye();
  mat oneN(n,n); oneN.ones(); oneN /= n;
  mat H = I - oneN;

  // 3-2. B
  mat B = (-0.5)*(H*D2*H.t());

  // 3-3. eigenvalue decomposition
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, B);

  // 3-4. Return !
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}


// 4. ICA : (Fast) Independent Component Analysis
// tnum = 1 (logcosh), 2 (exp), 3 (poly)
// sym  = FALSE (just run) / TRUE (decorrelate)

typedef vec (*icaPtr)(const vec& x, const double tpar);
vec ica_logcosh(const vec& x, const double tpar){
  vec y = tanh(tpar*x);
  return(y);
}
vec ica_logcoshp(const vec& x, const double tpar){
  vec y = tpar*(1-pow(tanh(tpar*x),2));
  return(y);
}
vec ica_exp(const vec& x, const double tpar){
  const int n = x.n_elem;
  vec y(n,fill::zeros);
  for (int i=0;i<n;i++){
    y(i) = x(i)*exp(-tpar*pow(x(i),2)/2);
  }
  //vec y =x*exp(-tpar*pow(x,2)/2);
  return(y);
}
vec ica_expp(const vec& x, const double tpar){
  const int n = x.n_elem;
  vec y(n,fill::zeros);
  for (int i=0;i<n;i++){
    double x2 = pow(x(i),2);
    y(i) = (1-tpar*x2)*exp((-tpar*x2)/2);
  }
  /*vec ones;// exp is the problem
  ones.copy_size(x);
  ones.ones();
  vec x2 = pow(x,2);
  vec y = (ones-tpar*x2)*exp((-tpar*x2)/2);*/
  return(y);
}
vec ica_poly(const vec& x, const double tpar){
  vec y = pow(x,3);
  return(y);
}
vec ica_polyp(const vec& x, const double tpar){
  vec y = 3*(pow(x,2));
  return(y);
}

XPtr<icaPtr> decideICAg(const int n){
  if (n==1){
    return(XPtr<icaPtr>(new icaPtr(&ica_logcosh)));
  } else if (n==2){
    return(XPtr<icaPtr>(new icaPtr(&ica_exp)));
  } else if (n==3){
    return(XPtr<icaPtr>(new icaPtr(&ica_poly)));
  } else {
    return XPtr<icaPtr>(R_NilValue);
  }
}
XPtr<icaPtr> decideICAgprime(const int n){
  if (n==1){
    return(XPtr<icaPtr>(new icaPtr(&ica_logcoshp)));
  } else if (n==2){
    return(XPtr<icaPtr>(new icaPtr(&ica_expp)));
  } else if (n==3){
    return(XPtr<icaPtr>(new icaPtr(&ica_polyp)));
  } else {
    return XPtr<icaPtr>(R_NilValue);
  }
}

//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_ica(arma::mat& X, const int C, const int maxiter,
                      const double tol, const int tnum, const double tpar,
                      bool sym){
  // 1. setting
  const int N = X.n_rows;
  const int M = X.n_cols;
  mat W       = 0.01*(randu<mat>(N,C));
  for (int i=0;i<C;i++){
    W.col(i) /= norm(W.col(i),2);
  }

  // 2. use function pointer for case branching
  XPtr<icaPtr> xpfun1 = decideICAg(tnum);
  XPtr<icaPtr> xpfun2 = decideICAgprime(tnum);

  icaPtr trf_g  = *xpfun1; // transformation for g
  icaPtr trf_gp = *xpfun2; // transformation for gprime

  vec onesM(M);
  onesM.ones();

  // 3. main computation
  for (int it=0;it<C;it++){
    vec wold = W.col(it);
    vec wnew, wtX(N), wtX1(N), wtX2(N);
    double tgtgap  = 100.0;
    int    tgtiter = 0;

    while (tgtgap > tol){
      wtX = (X.t()*wold);
      wtX1 = trf_g(wtX,tpar);
      wtX2 = trf_gp(wtX,tpar);

      wnew = (X*(wtX1))/M - as_scalar((wtX2.t()*onesM))*wold/M;
      //wnew = randu<vec>(N);
      if (it>0){
        for (int i=0;i<(it-1);i++){
          wnew -= as_scalar((wnew.t())*W.col(i))*W.col(i);
        }
      }
      wnew /= norm(wnew);

      tgtiter += 1;
      tgtgap = norm(wnew-wold);

      wold = wnew;
      if (tgtiter>maxiter){
        break;
      }
    }
    W.col(it) = wold;
  }

  // 4. symmetric decorrelation
  if (sym==false){
    mat S = (W.t())*X;
    return Rcpp::List::create(Rcpp::Named("S")=S,
                              Rcpp::Named("W")=W);
  } else if (sym==true){
    mat WWT = W*W.t();
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, WWT);
    mat preW = eigvec*diagmat(1/sqrt(eigval))*eigvec.t();
    mat Wnew = preW*W;

    mat S = (Wnew.t())*X;
    return Rcpp::List::create(Rcpp::Named("S")=S,
                              Rcpp::Named("W")=Wnew);
  } else {
    Rcpp::stop(" fastica : sym flag is incorrect.");
  }
}


// 5. rpgauss : Random Projection for Gaussian Case
//    Note that in this case, I simply followed R manner of data arrangement
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_rpgauss(arma::mat& X, const int k){
  // 5-1. setting
  const int n = X.n_rows;
  const int d = X.n_cols;
  mat R(d,k);
  mat G = (randu<mat>(d,k));

  // 5-2. generate and iterate
  R.col(0) = G.col(0)/norm(G.col(0));

  for (int i=1;i<k;i++){
    vec tgt = G.col(i);
    for (int j=0;j<i;j++){
      vec u = R.col(j);
      tgt -= (dot(u,tgt)/dot(tgt,tgt))*u;
    }
    R.col(i) = tgt/norm(tgt);
  }

  // 5-3. Compute
  mat Y = X*R;
  return Rcpp::List::create(Rcpp::Named("R")=R,
                            Rcpp::Named("Y")=Y);
}

/*
 * 6. Factor Analysis
 *    Notations are consistent with that of Gharahmani's paper.
 */
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_fa(arma::mat& X, const int k, const int maxiter, const double tolerance){
  // 6-1. basic settings
  const int p = X.n_rows;
  const int n = X.n_cols;

  // 6-2. other defined materials
  mat Ez(k,n);
  mat beta(k,p);
  mat LPhi(p,p);
  mat eyeK(k,k,fill::eye);
  mat Pinv(p,p,fill::zeros);
  mat EzzSum(k,k);
  double inctol = 0;

  // 6-3. initialization
  mat Lold = randn<mat>(p,k);
  mat Lnew(p,k,fill::zeros);
  vec Pold = ones<vec>(p);
  mat Ptmp(p,p,fill::zeros);
  vec Pnew(p);

  // 6-4. main iteration
  for (int it=0;it<maxiter;it++){
    // 6-4-E1. common again
    Pinv = diagmat(1/Pold);
    // 6-4-E2. inverse of LL^T+\Phi : LPhi
    LPhi = Pinv - Pinv*Lold*solve(eyeK+Lold.t()*Pinv*Lold,Lold.t()*Pinv);
    // 6-4-E3. beta
    beta = Lold.t()*LPhi;
    // 6-4-E4. EZ = [EZ(x1), EZ(x2), ... , EZ(xn)]
    Ez = beta*X;
    // 6-4-E5. EZZsum = sum(EZZ(xi))
    EzzSum = n*(eyeK-beta*Lold) + beta*X*X.t()*beta.t();

    // 6-4-M1. update Lambda : Lnew
    Lnew = (X*Ez.t())*pinv(EzzSum);
    // 6-4-M2. update Phi    : Pnew from Ptmp
    Ptmp = (X - Lnew*Ez)*X.t()/n;
    Pnew = Ptmp.diag();

    // 6-4-Update
    inctol = norm(Lnew-Lold,"fro");
    Lold = Lnew;
    Pold = Pnew;
    if (inctol < tolerance){
      break;
    }
  }

  // 6-5. MLE solution for Z
  mat Z(k,n);
  Pinv = diagmat(1/Pold);
  Z = solve(Lold.t()*Pinv*Lold,Lold.t()*sqrt(Pinv)*X);

  // 6-6. return results
  return Rcpp::List::create(Rcpp::Named("L")=Lold,
                            Rcpp::Named("Z")=Z,
                            Rcpp::Named("Pvec")=Pold);
}

/*
 * 7. Locality Preserving Projections (LPP)
 */
// [[Rcpp::export]]
Rcpp::List method_lpp(arma::mat& X, arma::mat& W){
  // 7-1. basic settings
  const int n = X.n_rows;
  const int m = X.n_cols;
  const int w1 = W.n_rows;
  const int w2 = W.n_cols;

  if (w1!=w2){
    Rcpp::stop("ERROR : W is not a square matrix.");
  }
  if (m!=w2){
    Rcpp::stop("ERROR : two inputs are not matching.");
  }

  // 7-2. main computation
  // 7-2-1. Degree matrix D and Laplacian L
  mat D(m,m);
  for (int i=0;i<m;i++){
    D(i,i) = as_scalar(sum(W.row(i)));
  }
  mat L = D-W;

  // 7-2-2. LHS and RHS
  mat LHS = X*L*X.t();
  mat RHS = X*D*X.t();
  mat SOL = solve(RHS, LHS);

  // 7-2-3. eigendecomposition
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, SOL);

  // 7-3. return results
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}


/*
 * 08. Neighborhood Preserving Embedding (NPE)
 */
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List method_npe(arma::mat& X, arma::mat& W){
  // 08-1. basic settings
  const int n = X.n_rows;
  const int m = X.n_cols;
  const int w1 = W.n_rows;
  const int w2 = W.n_cols;

  if (w1!=w2){
    Rcpp::stop("ERROR : W is not a square matrix.");
  }
  if (m!=w2){
    Rcpp::stop("ERROR : two inputs are not matching.");
  }

  // 08-2. preliminary items
  mat Im(m,m,fill::eye);
  mat M = ((Im-W).t())*(Im-W);

  // 08-3. LHS, RHS and SOL
  mat LHS = X*M*X.t();
  mat RHS = X*X.t();
  mat SOL = solve(RHS, LHS);

  // 08-4. main computation
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, SOL);

  // 08-5. return results
  return Rcpp::List::create(Rcpp::Named("eigval")=eigval,
                            Rcpp::Named("eigvec")=eigvec);
}
