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

// [[Rcpp::export]]
NumericMatrix UniFrac_UN (NumericMatrix OTUtab, NumericVector br_len)
{
  int n = OTUtab.nrow();
  NumericMatrix dist_mat(n, n);
  for(int i=0; i<=(n-2); i++){
    for(int j=i+1; j<=(n-1); j++) {
      dist_mat(i,j) = sum(br_len * abs((OTUtab(i,_) - OTUtab(j,_))))/
        sum(br_len * (OTUtab(i,_) + OTUtab(j,_) - OTUtab(i,_) * OTUtab(j,_)));
      dist_mat(j,i) = dist_mat(i,j);
    }
  }
  return dist_mat;
}

// [[Rcpp::export()]]
arma::mat UniFrac_G (arma::mat OTUtab, arma::rowvec br_len, double alpha) {
  const int n_sub = OTUtab.n_rows;
  arma::mat dist_mat(n_sub,n_sub); dist_mat.fill(0);
  int p = OTUtab.n_cols; 
  arma::rowvec sub_j(p), sub_i(p), temp_sum(p), temp_diff(p); arma::uvec pos_ids;
  for (int i=0; i<=(n_sub-2); i++) {
    sub_i = OTUtab.row(i);
    for (int j=(i+1); j<=(n_sub-1); j++) {
      sub_j = OTUtab.row(j);
      temp_sum = sub_i+sub_j; temp_diff = abs(sub_i-sub_j);
      pos_ids = arma::find(temp_sum > 0);
      dist_mat(i,j) = arma::accu((br_len.elem(pos_ids) % arma::pow(temp_sum.elem(pos_ids), alpha-1)) %
        temp_diff.elem(pos_ids)) / arma::accu(br_len.elem(pos_ids) % arma::pow(temp_sum.elem(pos_ids),alpha));
    }
  }
  dist_mat = dist_mat+trans(dist_mat);
  return dist_mat;
}

