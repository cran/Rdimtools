#ifndef _Rdimtools_METHODS_NONLINEAR_H
#define _Rdimtools_METHODS_NONLINEAR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// 1. SNE
arma::mat method_sne(arma::mat& P, const int ndim, const double eta,
                     const int maxiter, double jitter, double decay,
                     const double momentum);
// 2. Symmetric SNE
arma::mat method_snesym(arma::mat& P, const int ndim, const double eta,
                     const int maxiter, double jitter, double decay,
                     const double momentum);
// 3. t-SNE
arma::mat method_tsne(arma::mat& P, const int ndim, const double eta,
                        const int maxiter, double jitter, double decay,
                        const double momentum);
// 4. Eigenmaps
Rcpp::List method_eigenmaps(arma::mat& W);
// 5. Sammon Mapping
arma::mat method_sammon(arma::mat& X, arma::mat& init);
// 6. LLE
arma::vec method_lleW(arma::mat& mat_tgt, arma::vec& vec_tgt, const double regparam);
// 7. LLE with automatic choice
Rcpp::List method_lleWauto(arma::mat& mat_tgt, arma::vec& vec_tgt);
// 8. LLE M
Rcpp::List method_lleM(arma::mat& W);

#endif
