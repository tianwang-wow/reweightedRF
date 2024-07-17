source('func.R')

library(ranger) # for random forest
library(castor) # for microbiome distance calculation
library(Rsolnp) # for finding ensemble weights
library(quadprog) # for finding ensemble weights
library(pROC) # for AUC calculation
library(Rcpp)
library(RcppArmadillo)

### Load "dataset" (600 samples of 856 OTU counts + outcome) and "tree" (a phylogenetic tree)
load('./Data/Simulated_data.RData')

OTUtab = dataset[, 1:856] # OTU count table
dataset[,1:856] = dataset[,1:856] / rowSums(dataset[,1:856]) # convert OTU counts to abundances
TR_ind = 1:300 # training samples
TR_data = dataset[TR_ind, ] # training data
TE_data_withY = dataset[- TR_ind, ] # testing data with outcome
TE_data = TE_data_withY[, 1:856] # testing data without outcome

### Calculate different microbiome distances of both training and testing samples
distmat_ls = UniFrac_dist(OTUtab = OTUtab, # OTU count table
                          type = c('BC', 'Unweighted', 0, 0.5, 'Weighted'), # microbiome distances to calculate
                          tree = tree)
TR_distmat_ls = lapply(distmat_ls, function(x) x[TR_ind, TR_ind]) # distances for training samples

### Fit ensemble estimate
set.seed(123)
fit = ensemble(TR_data = TR_data, # data to fit
               TR_distmat_ls = TR_distmat_ls, # list of microbiome distances
               n_CV = 5, # number of folds in cross-validataions
               eta_vec = seq(0.1, 1, 0.1), # possible values of hyper-parameter for gaussian kernel
               outcome = 'Y', # outcome name in TR_data
               outcome_type = 'binary', # continuous or binary outcome
               num.trees = 1000, # number of trees in RF
               min.node.size = 10 # min.node.size in RF
               )
               
fit$ens_weights # ensemble weights
fit$optimal_eta # optimal etas

### Do prediction on testing samples
pred_ensemble = predict_ensemble(fit = fit, # ensemble fit
                                 distmat_ls = distmat_ls, # list of microbiome distances of both training and testing samples
                                 TR_ind = TR_ind, # ids of training samples in the above list of distance matrices
                                 TE_data = TE_data # testing data
                                 ) 

### Classification error
apply(pred_ensemble, 2, function(x) mean(ifelse(x >= 0.5, 1, 0) != TE_data_withY[, 'Y']))
