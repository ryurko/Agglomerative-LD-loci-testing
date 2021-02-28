#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_cor_matrix(arma::mat X) {
  int ncols = X.n_cols;
  
  // Initialize the correlation matrix
  arma::mat cor_matrix(ncols, ncols);
  
  // Calculate it:
  cor_matrix = arma::cor(X, 1);
  return cor_matrix;
}

// [[Rcpp::export]]
arma::vec compute_eig_vals(arma::mat cor_matrix) {
  arma::vec eig_val; // vector of eigenvalues
  arma::eig_sym(eig_val, cor_matrix);
  
  return eig_val;
}


// [[Rcpp::export]]
arma::vec sim_ev_quad_stats(arma::vec eig_val,
                             int n_sims) {
  
  int n_snps = eig_val.n_rows;
  
  // Initialize the vector to hold the quadratic test statistics
  // for the given number of simulations and blocks:
  arma::vec quad_stats = arma::zeros(n_sims);
  
  // Now loop through the number of blocks to generate the simulations for the
  // test statistics:
  for (int i = 0; i < n_sims; i++) {
    arma::vec z = arma::randn(n_snps);
    quad_stats.row(i) = eig_val.t() * square(z);
    
  }
  
  return quad_stats;
}
