#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export()]] //OOB weight counts used in 'PhyEnRF_fit'
arma::mat OOB_weight_cts (arma::mat inbag_counts,
                          arma::mat terminal_nodes) {
  int n_sub = inbag_counts.n_rows; // number of subjects
  arma::uvec RFindex, sampleids, i_vec, j_vec; int j;
  arma::mat weight_counts(n_sub, n_sub); weight_counts.fill(0); // initialize distances to 0
  for (int i=0; i<=(n_sub-1); i++) {
    i_vec = i;
    RFindex = arma::find(inbag_counts.row(i) == 0); //which trees such that i-th sample in data_x is not used
    if (RFindex.n_elem >= 1) { // 09/29/2021: Add this "if"
      for (int k=0; k<=(RFindex.n_elem-1); k++) {
        j_vec = RFindex(k); j = RFindex(k);
        sampleids = arma::find(terminal_nodes.col(j) == terminal_nodes(i,j)); // contributed old sample ids to predict the i-th newsample
        weight_counts.submat(i_vec, sampleids) = weight_counts.submat(i_vec, sampleids) +
          trans(inbag_counts.submat(sampleids, j_vec) / arma::accu(inbag_counts.submat(sampleids, j_vec)));
      }
    }
  }
  return weight_counts;
}


// [[Rcpp::export()]] //OOB weight counts used in 'PhyEnRF_fit'
arma::mat new_weight_cts (arma::mat inbag_counts,
                          arma::mat terminal_nodes,
                          arma::mat new_terminal_nodes) {
  int n_old = terminal_nodes.n_rows; // number of old subjects
  int n_new = new_terminal_nodes.n_rows; // number of new subjects

  arma::mat new_wt_cts(n_new, n_old); new_wt_cts.fill(0);
  arma::uvec sampleids, i_vec, k_vec;
  // int n_tree = terminal_nodes.n_cols;
  for (int i=0; i<=(n_new-1); i++) {
    i_vec = i;
    for (int k=0; k<=(terminal_nodes.n_cols-1); k++) { //for each tree
      k_vec = k;
      sampleids = arma::find(terminal_nodes.col(k) == new_terminal_nodes(i,k));//contributed old sample ids to predict the i-th newsample
      new_wt_cts.submat(i_vec, sampleids) = new_wt_cts.submat(i_vec, sampleids) +
        trans(inbag_counts.submat(sampleids, k_vec) / arma::accu(inbag_counts.submat(sampleids, k_vec)));
    }
  } //note: sum(weight_counts[1,]) == n_tree
  return new_wt_cts;
}



