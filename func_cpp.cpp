#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export()]]
arma::mat new_weight_cts (arma::mat inbag_counts,
                          arma::mat terminal_nodes,
                          arma::mat new_terminal_nodes) {
  int n_old = terminal_nodes.n_rows;
  int n_new = new_terminal_nodes.n_rows;

  arma::mat new_wt_cts(n_new, n_old); new_wt_cts.fill(0);
  arma::uvec sampleids, i_vec, k_vec;
  for (int i=0; i<=(n_new-1); i++) {
    i_vec = i;
    for (int k=0; k<=(terminal_nodes.n_cols-1); k++) {
      k_vec = k;
      sampleids = arma::find(terminal_nodes.col(k) == new_terminal_nodes(i,k));
      new_wt_cts.submat(i_vec, sampleids) = new_wt_cts.submat(i_vec, sampleids) +
        trans(inbag_counts.submat(sampleids, k_vec) / arma::accu(inbag_counts.submat(sampleids, k_vec)));
    }
  }
  return new_wt_cts;
}
