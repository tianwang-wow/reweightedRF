reweighted_RF = function(RF, # a fitted random forest
                         data_x, # a matrix of input features: row: subjects; column: features
                         data_y, # a vector of responses: either continuous or binary class
                         Kmat_ls, # a list of kernel matrices (not needed if pred = F)
                         pred = F) # whether to find predictions on training samples
{
  if (is.factor(data_y)) { data_y = as.numeric(as.character(data_y)) }
  n_sam = RF$num.samples # number of all samples
  inbag_counts = matrix(unlist(RF$inbag.counts), nrow = n_sam, byrow = FALSE) # row: subject; column: tree;
  terminal_nodes = predict(object = RF, data = data.frame(data_x), type = 'terminalNodes', predict.all = TRUE)$predictions
  
  pred_mat = NULL
  if (pred) { # if we want to prediction on TR samples
    weight_counts = OOB_weight_cts(inbag_counts, terminal_nodes)
    Kmat_RF = weight_counts / rowSums(weight_counts) # kernel from RF
    pred_mat = Kmat_RF %*% data_y # predicted prob on TR set using RF
    colnames(pred_mat)[1] = 'Original' # record this prediction prob
    
    n_K = length(Kmat_ls) # number of distance matrices
    for (l_id in 1:n_K) { # for each kernel
      Kmat = Kmat_ls[[l_id]] # current Kmat
      temp_weight = Kmat * Kmat_RF # weight matrix for prediction
      temp_weight = temp_weight / rowSums(temp_weight) # scale the weight matrix
      temp_predictions = temp_weight %*% data_y # predicted prob in TR set
      pred_mat = cbind(pred_mat, temp_predictions) # cbind predcited prob
    }
    colnames(pred_mat)[-1] = names(Kmat_ls) # add method names
  }
  
  out = list(RF = RF,
             inbag_counts = inbag_counts, # row: subject; column: tree;
             terminal_nodes = terminal_nodes,
             data_x = data_x,
             data_y = data_y,
             pred_mat = pred_mat)
  return(out)
}


predict_reweighted_RF = function(fit, # a fitted reweighted_RF
                                 new_data_x, # a matrix of input features: row: subjects; column: features
                                 new_Kmat_ls) # a list of kernel values: row TE x column TR samples
{
  new_terminal_nodes = predict(object = fit$RF,
                               data = data.frame(new_data_x), 
                               type = 'terminalNodes', 
                               predict.all = TRUE)$predictions
  new_wt = new_weight_cts(fit$inbag_counts, 
                          fit$terminal_nodes, 
                          new_terminal_nodes)
  Kmat_RF = new_wt / rowSums(new_wt) # kernel from RF
  
  pred_mat = NULL
  pred_mat = Kmat_RF %*% fit$data_y # predicted prob on new data set using RF
  colnames(pred_mat)[1] = 'Original' # record this prediction prob
  n_K = length(new_Kmat_ls) # number of kernels
  for (l_id in 1:n_K) { # for each kernel
    new_Kmat = new_Kmat_ls[[l_id]] # current Kmat: row TE x column TR
    temp_weight = new_Kmat * Kmat_RF # weight matrix for prediction
    temp_weight = temp_weight / rowSums(temp_weight) # scale the weight matrix
    temp_predictions = temp_weight %*% fit$data_y # predicted prob in TR set
    pred_mat = cbind(pred_mat, temp_predictions) # cbind predcited prob
  }
  colnames(pred_mat)[-1] = names(new_Kmat_ls) # add method names
  
  out = list(pred_mat = pred_mat)
  
  return(out)
}
