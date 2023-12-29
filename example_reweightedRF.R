library(ranger) # for random forest
library(Rcpp)
library(RcppArmadillo)

source('func.R')

### Generate data
set.seed(123)
n = 1000
p = 50
data_x = matrix(nrow = n, rnorm(n * 100))
prob = 1/(1 + 0.5 * exp(- data_x[, c(1, 100)]))
data_y = rbinom(n, 1, prob)
dataset = cbind(data_x[, 1:p], data_y)
colnames(dataset) = c(paste0('X', 1:p), 'Y')
TR_ind = 1:500 # training samples

### Fit a regular random forest using training samples
RF = ranger(dependent.variable.name = "Y",
            data = data.frame(dataset[TR_ind, ]),
            keep.inbag = T,
            probability = F, # must be FALSE for reweighted RF
            num.trees = 1000,
            classification = T)

### Fit a reweightedRF
fit_Re_RF = reweighted_RF(RF = RF,
                          data_x = dataset[TR_ind, 1:p],
                          data_y = dataset[TR_ind, 'Y'])
### Get some kernel values between testing and training samples. You may use your preferred kernels.
eta = c(0.1, 0.5, 1) # 3 hyper-parameters
distmat = as.matrix(dist(data_x[, c(1, 100)]))
mean_dist = mean(distmat[TR_ind, TR_ind][upper.tri(distmat[TR_ind, TR_ind])])
Kmat_ls = list()
for (l_id in 1:length(eta)) { Kmat_ls[[l_id]] = exp(- distmat^2 / 2 / (mean_dist * eta[l_id])^2) }
names(Kmat_ls) = paste0('eta = ', eta)
new_Kmat_ls_RF = lapply(Kmat_ls, function(x) x[-TR_ind, TR_ind]) # a list of 3 kernel matrices: row: testing samples x column: training samples

### Find predictions on testing samples using reweightedRF
pred_Re_RF = predict_reweighted_RF(fit = fit_Re_RF,
                                   new_data_x = dataset[- TR_ind, 1:p],
                                   new_Kmat_ls = new_Kmat_ls_RF)
pred_Re_RF$pred_mat[1:5, ] # 1st column includes predictions for testing samples using the originial RF; the other 3 columns are predictions using reweightedRF

### Classification error
apply(pred_Re_RF$pred_mat, 2, function(x) mean(ifelse(x >= 0.5, 1, 0) != dataset[-TR_ind, 'Y']))
