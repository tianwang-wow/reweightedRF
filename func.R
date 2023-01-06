Rcpp::sourceCpp('func_cpp.cpp')

reweighted_RF = function(RF, # a fitted random forest
                         data_x, # a matrix of input features: row: subjects; column: features
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


predict_reweighted_RF = function(fit, # a fitted reweighted_RF
                                 new_data_x, # a matrix of input features: row: new subjects; column: features
                                 new_Kmat_ls # a list of kernel values: row TE x column TR samples
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
  
  pred_mat = NULL
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
