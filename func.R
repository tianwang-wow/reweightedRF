Rcpp::sourceCpp('func_cpp.cpp')

reweighted_RF = function(RF, # a fitted random forest from the R package "ranger"
                         data_x, # a matrix of input features: row: samples x column: features
                         data_y # a vector of responses: either continuous or binary
)
{
  n_sam = RF$num.samples
  inbag_counts = matrix(unlist(RF$inbag.counts), nrow = n_sam, byrow = FALSE)
  terminal_nodes = predict(object = RF, data = data.frame(data_x), type = 'terminalNodes', predict.all = TRUE)$predictions
  
  out = list(RF = RF,
             inbag_counts = inbag_counts, # inbag counts
             terminal_nodes = terminal_nodes, # terminal nodes
             data_x = data_x,
             data_y = data_y)
  return(out)
}

predict_reweighted_RF = function(fit, # a fitted reweightedRF
                                 new_data_x, # a matrix of input features: row: new samples x column: features
                                 new_Kmat_ls # a list of kernel values: row testing x column training samples
)
{
  if (is.factor(fit$data_y)) { data_y = as.numeric(as.character(fit$data_y)) }
  new_terminal_nodes = predict(object = fit$RF,
                               data = data.frame(new_data_x), 
                               type = 'terminalNodes', 
                               predict.all = TRUE)$predictions
  new_wt = new_weight_cts(fit$inbag_counts, 
                          fit$terminal_nodes, 
                          new_terminal_nodes)
  Kmat_RF = new_wt / rowSums(new_wt) # kernel from RF
  
  pred_mat = Kmat_RF %*% fit$data_y
  colnames(pred_mat)[1] = 'Original'
  for (l_id in 1:length(new_Kmat_ls)) {
    new_Kmat = new_Kmat_ls[[l_id]]
    temp_weight = new_Kmat * Kmat_RF
    temp_weight = temp_weight / rowSums(temp_weight)
    temp_predictions = temp_weight %*% fit$data_y
    pred_mat = cbind(pred_mat, temp_predictions)
  }
  colnames(pred_mat)[-1] = names(new_Kmat_ls)
  
  out = list(pred_mat = pred_mat)
  
  return(out)
}

btm_ind = function(tree)
{
  n_node = tree$Nnode
  n_tip = length(tree$tip.label)
  subtree_ls = castor::get_subtrees_at_nodes(tree = tree, nodes = 1:n_node)
  
  btm_ind_mat = matrix(0, nrow = n_tip, ncol = n_node)
  colnames(btm_ind_mat) = paste0('Inter', (n_tip+1):(n_tip+n_node))
  rownames(btm_ind_mat) = tree$tip.label
  for (t_id in 1:n_node) { btm_ind_mat[subtree_ls$new2old_tips[[t_id]], t_id] = 1 }
  btm_ind_mat = Matrix::Matrix(btm_ind_mat, sparse = T)  
  
  return(btm_ind_mat)
}

extend_OTUtab = function(OTUtab,
                         tree,
                         btm_ind_mat,
                         scale = T)
{
  if (missing(btm_ind_mat)) { btm_ind_mat = btm_ind(tree) }
  if (sum(!colnames(OTUtab) %in% tree$tip.label) != 0) {
    stop('OTU table column names do not match the phylogenetic tree labels!')
  }
  OTUtab = as.matrix(OTUtab)
  OTUtab = OTUtab[, tree$tip.label, drop = F]
  
  if (scale) {
    scaled_OTUtab = OTUtab / rowSums(OTUtab)
  } else {
    scaled_OTUtab = OTUtab
  }
  extra_part = scaled_OTUtab %*% btm_ind_mat
  exp_OTUtab = cbind(scaled_OTUtab, extra_part)
  
  return(as.matrix(exp_OTUtab))
}

brl = function(tree) {
  n_OTUs = length(tree$tip.label)
  son_ind = tree$edge[,2]
  mom_ind = tree$edge[,1]
  root_id = unique(mom_ind[which(!mom_ind %in% son_ind)])
  taxa_names = c(tree$tip.label, paste0('Inter', (n_OTUs+1):max(son_ind)))
  
  br_len = setNames(c(tree$edge.length, 0),
                    c(taxa_names[son_ind], taxa_names[root_id]))
  br_len = br_len[taxa_names]
  return(br_len)
}

UniFrac_W = function(exp_OTUtab,
                     br_len)
{
  n = nrow(exp_OTUtab)
  exp_OTUtab_br = t(t(exp_OTUtab) * br_len)
  top = as.matrix(stats::dist(exp_OTUtab_br, method = 'manhattan'))
  summation = rowSums(exp_OTUtab_br)
  btm = outer(summation, summation, function(x, y) { x+y })
  btm_inv = 1 / btm
  btm_inv[btm_inv == Inf] = 0
  dist_mat = top * btm_inv
  return(dist_mat)
}

UniFrac_BC = function(OTUtab = OTUtab, w = 2)
{ return(as.matrix(stats::dist(OTUtab, method = 'manhattan')) / w) }

UniFrac_dist = function(OTUtab,
                        type = c('BC', 'Unweighted', 0, 0.5, 'Weighted'),
                        tree)
{
  if (length(rownames(OTUtab)) != 0) { sample_names = rownames(OTUtab) } else { sample_names = 1:nrow(OTUtab) }
  btm_ind_mat = btm_ind(tree = tree)
  abu_exp_OTUtab = extend_OTUtab(OTUtab = OTUtab, tree = tree, btm_ind_mat = btm_ind_mat, scale = T) 
  cts_exp_OTUtab = extend_OTUtab(OTUtab = OTUtab, tree = tree, btm_ind_mat = btm_ind_mat, scale = F) 
  br_len = brl(tree)
  
  n_type = length(type)
  distmat_ls = setNames(vector('list', length = n_type), type)
  OTU_names = tree$tip.label
  OTUtab = OTUtab[, OTU_names, drop = F]
  abun_OTUtab = abu_exp_OTUtab[, OTU_names]
  
  for (k in 1:length(type)) {
    if (type[k] == 'Weighted') {
      distmat_ls[[k]] = UniFrac_W(abu_exp_OTUtab, br_len) 
    } else if (type[k] == 'Unweighted') {
      temp_OTUtab = cts_exp_OTUtab
      temp_OTUtab[temp_OTUtab > 0] = 1
      distmat_ls[[k]] = UniFrac_UN(temp_OTUtab, br_len)
    } else if (! is.na(suppressWarnings(as.numeric(type[k])))) {
      alpha = as.numeric(type[k])
      if (alpha == 1) {
        distmat_ls[[k]] = UniFrac_W(abu_exp_OTUtab, br_len) 
      } else {
        distmat_ls[[k]] = UniFrac_G(abu_exp_OTUtab, br_len, alpha)
      }
    } else if (type[k] == 'BC') {
      distmat_ls[[k]] = UniFrac_BC(abun_OTUtab)
    }
    rownames(distmat_ls[[k]]) = colnames(distmat_ls[[k]]) = sample_names
  }
  
  return(distmat_ls)
}

ensemble_weights = function(stacked_CV_pred_mat,
                            outcome_type,
                            selected_methods,
                            outcome = 'Y'
)
{
  n_methods = length(selected_methods)
  X_mat = stacked_CV_pred_mat[, selected_methods, drop = F]
  Y_vec = stacked_CV_pred_mat[, outcome]
  if (n_methods > 1) {
    if (outcome_type == 'binary') {
      rows_to_remove = which(rowMeans(X_mat) == 0 | rowMeans(X_mat) == 1)
      if (length(rows_to_remove) > 0) {
        X_mat = X_mat[-rows_to_remove, ]; Y_vec = Y_vec[ -rows_to_remove]
      }
      f = function(w) {
        w_prob = X_mat %*% w; return( - sum( Y_vec * log(w_prob) + (1-Y_vec) * log(1-w_prob)) )
      }
      eqfun = function(w) { sum(w) }
      opt_results = Rsolnp::solnp(pars = rep(1/n_methods, n_methods),
                                  fun = f,
                                  eqfun = eqfun,
                                  eqB = 1,
                                  LB = rep(0, n_methods),
                                  control = list(trace = 0))
      w_ens_vec = setNames(opt_results$pars, selected_methods)
    } else if (outcome_type == 'continuous') {
      dvec = t(X_mat) %*% Y_vec
      Dmat = t(X_mat) %*% X_mat
      Amat = t(rbind(rep(1, n_methods), diag(n_methods)))
      bvec = c(1, rep(0, n_methods))
      corrected_Dmat = Matrix::nearPD(Dmat)$mat
      QP_results = quadprog::solve.QP(Dmat = corrected_Dmat,
                                      dvec = dvec,
                                      Amat = Amat,
                                      bvec = bvec,
                                      meq = 1) 
      w_ens_vec = setNames(QP_results$solution, selected_methods)
    }
  } else {
    w_ens_vec = 1; names(w_ens_vec) = selected_methods
  }
  return(w_ens_vec)
}

rm_comma = function(x) { gsub(":.*", "", x) }

AUC = function(y, prob) { 
  if (TRUE %in% is.na(y) | TRUE %in% is.na(prob)) {
    out = -Inf
  } else if (length(unique(y)) == 1) {
    out = -Inf
  } else {
    out = pROC::roc(as.vector(y), as.vector(prob), levels = c(0,1), direction = '<')$auc
  }
  return(out[1])
}  

MSE = function(y, pred) { 
  if (TRUE %in% is.na(y) | TRUE %in% is.na(pred)) {
    out = Inf
  } else {
    out = mean((y-pred)^2)
  }
  return(out)
}  

D2K_gaussian = function(distmat_ls,
                        scale_ls)
{
  Kmat_ls = list()
  l_id = 1
  dist_names = names(distmat_ls)
  for (d_id in 1:length(distmat_ls)) {
    dist_mat = distmat_ls[[d_id]]
    scale_vec = scale_ls[[d_id]]
    for (scale_id in 1:length(scale_vec)) {
      scaling = scale_vec[scale_id]
      Kmat_ls[[l_id]] = exp(- dist_mat^2 / 2 / scaling^2 )
      names(Kmat_ls)[l_id] = paste0(dist_names[d_id], ':', scale_id)
      l_id = l_id + 1
    }
  }
  return(Kmat_ls)
}

ensemble = function(TR_data, # data to fit
                    TR_distmat_ls, # list of microbiome distances
                    n_CV = 5, # number of folds in cross-validataions
                    eta_vec = seq(0.1, 1, 0.1), # possible values of hyper-parameter for gaussian kernel
                    outcome = 'Y', # outcome name in TR_data
                    outcome_type = 'binary', # continuous or binary outcome
                    ...
                    )
{
  TR_data_y = TR_data[, outcome]
  TR_data_x = TR_data[, - which(colnames(TR_data) == outcome)]
  n_TR = nrow(TR_data)
  n_dist = length(TR_distmat_ls)
  dist_names = names(TR_distmat_ls)
  eta_ls = setNames(rep(list(eta_vec), n_dist), dist_names)
  while (T) {
    cv_ls = split(sample(1:n_TR), ceiling(seq(0, n_CV, length.out = n_TR+1)[-1]))
    if (outcome_type == 'continuous') {
      break
    } else if (min(sapply(cv_ls, function(x) {length(table(TR_data_y[x]))})) >= 2
               & min(sapply(cv_ls, function(x) {min(table(TR_data_y[x]))})) >= 2) {
      break
    }
  }
  opt_eta_mat_RF = matrix(NA, nrow = n_CV, ncol = n_dist)
  colnames(opt_eta_mat_RF) = dist_names
  opt_eta_mat_RF_id = opt_eta_mat_RF
  
  for (cv_id in 1:n_CV) {
    bs_ind = as.vector(unlist(cv_ls[-cv_id]))
    no_ind = as.vector(unlist(cv_ls[cv_id]))
    
    n_TR_sub = length(bs_ind)
    while (T) {
      sub_cv_ls = split(sample(bs_ind), ceiling(seq(0, n_CV, length.out = n_TR_sub+1)[-1]))
      if (outcome_type == 'continuous') {
        break
      } else if (min(sapply(sub_cv_ls, function(x) {length(table(TR_data_y[x]))})) >= 2) {
        break
      }
    }
    sub_AUC_mat_RF = NULL
    
    for (sub_cv_id in 1:n_CV) {
      sub_bs_ind = as.vector(unlist(sub_cv_ls[-sub_cv_id]))
      sub_no_ind = as.vector(unlist(sub_cv_ls[sub_cv_id]))
      sub_ind = c(sub_bs_ind, sub_no_ind)
      n_sub_bs = length(sub_bs_ind) 

      sub_TR_distmat_ls = lapply(TR_distmat_ls, function(x) x[sub_ind, sub_ind])
      sub_mean_vec = lapply(sub_TR_distmat_ls,
                            function(x) mean(x[1:n_sub_bs, 1:n_sub_bs][upper.tri(x[1:n_sub_bs, 1:n_sub_bs])]))
      sub_scale_ls = Map("*",
                         lapply(sub_mean_vec, function(x) rep(x, length(eta_vec))),
                         eta_ls)
      sub_Kmat_ls = D2K_gaussian(distmat_ls = sub_TR_distmat_ls,
                                 scale_ls = sub_scale_ls)
      sub_new_Kmat_ls = lapply(sub_Kmat_ls, function(x) x[-(1:n_sub_bs), 1:n_sub_bs])
      
      sub_RF = ranger(dependent.variable.name = outcome,
                      data = data.frame(TR_data[sub_bs_ind, ]),
                      keep.inbag = T,
                      probability = F,
                      classification = ifelse(outcome_type == 'binary', T, F),
                      ...
                      )
      sub_fitted_Re_RF = reweighted_RF(RF = sub_RF,
                                       data_x = TR_data_x[sub_bs_ind, ],
                                       data_y = TR_data_y[sub_bs_ind])
      sub_pred_mat_RF = predict_reweighted_RF(fit = sub_fitted_Re_RF,
                                              new_data_x = TR_data_x[sub_no_ind, ],
                                              new_Kmat_ls = sub_new_Kmat_ls)
      if (outcome_type == 'binary') {
        sub_AUC_vec_RF = apply(sub_pred_mat_RF$pred_mat, 2, function(x) AUC(TR_data_y[sub_no_ind], x))
      } else {
        sub_AUC_vec_RF = apply(sub_pred_mat_RF$pred_mat, 2, function(x) - MSE(TR_data_y[sub_no_ind], x))
      }
      sub_AUC_mat_RF = rbind(sub_AUC_mat_RF, sub_AUC_vec_RF)
    }
    sub_AUC_mean_RF = colMeans(sub_AUC_mat_RF)
    current_all_methods = names(sub_AUC_mean_RF)
    for (d_id in 1:n_dist) {
      current_dist = dist_names[d_id]
      temp_id = which(grepl(paste0(current_dist, ':'), current_all_methods))
      temp_AUCs = sub_AUC_mean_RF[temp_id]
      max_id = which(temp_AUCs == max(temp_AUCs))
      if (length(max_id) > 1) { max_id = sample(max_id)[1] }
      opt_eta_mat_RF[cv_id, current_dist] = eta_vec[max_id]
      opt_eta_mat_RF_id[cv_id, current_dist] = max_id
    }
  }
  
  AUC_mat_RF = NULL
  stacked_CV_pred_mat = NULL
  for (cv_id in 1:n_CV) {
    bs_ind = as.vector(unlist(cv_ls[-cv_id]))
    no_ind = as.vector(unlist(cv_ls[cv_id]))
    n_TR_sub = length(bs_ind)
    re_TR_distmat_ls = lapply(TR_distmat_ls, function(x) x[c(bs_ind, no_ind), c(bs_ind, no_ind)])
    
    CV_RF = ranger(dependent.variable.name = outcome,
                   data = data.frame(TR_data[bs_ind, ]),
                   keep.inbag = T,
                   probability = F,
                   classification = ifelse(outcome_type == 'binary', T, F),
                   ...
                   )
    CV_fitted_Re_RF = reweighted_RF(RF = CV_RF,
                                    data_x = TR_data_x[bs_ind, ],
                                    data_y = TR_data_y[bs_ind])
    CV_scale_ls = Map("*",
                      lapply(re_TR_distmat_ls,
                             function(x) mean(x[1:n_TR_sub, 1:n_TR_sub][upper.tri(x[1:n_TR_sub, 1:n_TR_sub])])),
                      eta_ls)
    
    CV_Kmat_ls_RF = D2K_gaussian(distmat_ls = re_TR_distmat_ls,
                                 scale_ls = CV_scale_ls)
    CV_new_Kmat_ls_RF = lapply(CV_Kmat_ls_RF, function(x) x[-c(1:n_TR_sub), 1:n_TR_sub])
    CV_pred_Re_RF = predict_reweighted_RF(fit = CV_fitted_Re_RF,
                                          new_data_x = TR_data_x[no_ind, ],
                                          new_Kmat_ls = CV_new_Kmat_ls_RF)
    CV_pred_mat_RF = CV_pred_Re_RF$pred_mat
    colnames(CV_pred_mat_RF) = paste0('RF_', colnames(CV_pred_mat_RF))
    
    CV_pred_mat = cbind(CV_pred_mat_RF, TR_data_y[no_ind])
    colnames(CV_pred_mat)[ncol(CV_pred_mat)] = outcome
    current_selected_methods = c('RF_Original',
                                 paste0('RF_', dist_names, ':', opt_eta_mat_RF_id[cv_id, ]),
                                 outcome)
    stacked_CV_pred_mat = rbind(stacked_CV_pred_mat,
                                CV_pred_mat[, current_selected_methods])
    
    if (outcome_type == 'binary') {
      AUC_vec_RF = apply(CV_pred_mat[, -ncol(CV_pred_mat)], 2, function(x) AUC(TR_data_y[no_ind], x))
    } else {
      AUC_vec_RF = apply(CV_pred_mat[, -ncol(CV_pred_mat)], 2, function(x) - MSE(TR_data_y[no_ind], x))
    }
    AUC_mat_RF = rbind(AUC_mat_RF, AUC_vec_RF)
  }
  colnames(stacked_CV_pred_mat) = rm_comma(colnames(stacked_CV_pred_mat))
  w_ens_vec = ensemble_weights(stacked_CV_pred_mat = stacked_CV_pred_mat,
                               outcome_type = outcome_type,
                               selected_methods = paste0('RF_', c('Original', dist_names)),
                               outcome = outcome)
  w_ens_vec_ls = list(Ensemble_RF = w_ens_vec)
  
  AUC_mean_RF = colMeans(AUC_mat_RF)
  final_opt_eta_id_vec = vector()
  for (d_id in 1:n_dist) {
    current_dist = dist_names[d_id]
    temp_id = which(grepl(paste0(current_dist, ':'), colnames(AUC_mat_RF)))
    temp_AUCs = AUC_mean_RF[temp_id]
    max_id = which(temp_AUCs == max(temp_AUCs))
    if (length(max_id) > 1) { max_id = sample(max_id)[1] }
    final_opt_eta_id_vec[ d_id ] = max_id
  }
  names(final_opt_eta_id_vec) = dist_names
  final_opt_eta_vec = setNames(eta_vec[ final_opt_eta_id_vec ],
                               dist_names)
  
  RF = ranger(dependent.variable.name = outcome,
              data = data.frame(TR_data),
              keep.inbag = T,
              probability = F,
              classification = ifelse(outcome_type == 'binary', T, F),
              ...
              )
  fitted_Re_RF = reweighted_RF(RF = RF,
                               data_x = TR_data_x,
                               data_y = TR_data_y)
  scale_ls = Map("*",
                 lapply(TR_distmat_ls,
                        function(x) mean(x[upper.tri(x)])),
                 as.list(final_opt_eta_vec))
  
  return(list(RF = RF,
              fitted_Re_RF = fitted_Re_RF,
              scale_ls = scale_ls,
              optimal_eta = final_opt_eta_vec,
              ens_weights = w_ens_vec_ls$Ensemble_RF
  ))
}

predict_ensemble = function(fit,
                            distmat_ls,
                            TR_ind,
                            TE_data
                            ) 
{
  scale_ls = fit$scale_ls
  fitted_Re_RF = fit$fitted_Re_RF
  ens_weights = fit$ens_weights
  
  Kmat_ls_RF = D2K_gaussian(distmat_ls = distmat_ls,
                            scale_ls = scale_ls)
  new_Kmat_ls_RF = lapply(Kmat_ls_RF, function(x) x[-TR_ind, TR_ind])
  pred_Re_RF = predict_reweighted_RF(fit = fitted_Re_RF,
                                     new_data_x = TE_data,
                                     new_Kmat_ls = new_Kmat_ls_RF)
  TE_pred_mat_RF = pred_Re_RF$pred_mat
  colnames(TE_pred_mat_RF) = paste0('RF_', rm_comma(colnames(TE_pred_mat_RF)))
  
  TE_pred_mat = cbind(TE_pred_mat_RF, 
                      TE_pred_mat_RF[, names(ens_weights), drop = F] %*% ens_weights)
  colnames(TE_pred_mat)[ncol(TE_pred_mat)] = 'Ensemble'
  
  return(TE_pred_mat)
}

