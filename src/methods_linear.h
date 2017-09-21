#ifndef _Rdimtools_METHODS_LINEAR_H
#define _Rdimtools_METHODS_LINEAR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// TO DO : method_mdsSD - sparse distance matrix !
// 1. PCA
Rcpp::List method_pca(arma::mat& psdX);
// 2. MDS
Rcpp::List method_mds(arma::mat& centerX);
// 3. MDS given D
Rcpp::List method_mdsD(arma::mat& D);
// 4. ICA
Rcpp::List method_ica(arma::mat& X, const int C, const int maxiter, const double tol, const int tnum, const double tpar, bool sym);
// 5. RNDPROJ
Rcpp::List method_rpgauss(arma::mat& X, const int k);
// 6. FA
Rcpp::List method_fa(arma::mat& X, const int k, const int maxiter, const double tolerance);
// 7. LPP
Rcpp::List method_lpp(arma::mat& X, arma::mat& W);
// 8. NPE
Rcpp::List method_npe(arma::mat& X, arma::mat& W);

#endif
