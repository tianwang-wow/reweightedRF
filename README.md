 
================

## File list

The zip file in the supplement includes the following files:

- **`example_reweightedRF.R`**: an example of fitting reweighted RFs.
- **`example_ensemble.R`**: an example of fitting ensemble for
  microbiome data.
- **`func.R & func_cpp.cpp`**: two R files of functions for calculating
  microbiome distances and fitting reweighted RFs and ensemble.
- **`Data`**: a folder containing an R dataset
  **`Simulated_data.RData`** for running the example code.

## Main functions

<!-- UniFrac_dist -->
<!-- ensemble -->
<!-- predict_ensemble -->
<!-- reweighted_RF -->
<!-- predict_reweighted_RF -->

**`UniFrac_dist`**: Microbiome distance calculation.

- Argument
  - `OTUtab`: OTU count table (row: samples; column: OTUs).
  - `type`: A vector of distance names. Default = c(‘BC’, ‘Unweighted’,
    0, 0.5, ‘Weighted’).
  - `tree`: Phylogenetic tree.
- Value
  - A list of distance matrices.

**`ensemble`**: Fit reweighted RFs and ensemble estimate.

- Argument
  - `TR_data`: Training data (row: samples; column: features + outcome).
    It would be OTU abundance data + outcome if using microbiome data.
  - `TR_distmat_ls`: A list of microbiome distances for samples in
    `TR_data`.
  - `n_CV`: Number of folds in cross-validataions. Default = 5.
  - `eta_vec`: A vector of possible values of hyper-parameter for
    gaussian kernel. Default = seq(0.1, 1, 0.1).
  - `outcome`: Outcome name in `TR_data`.
  - `outcome_type`: ‘binary’ or ‘continuous’ outcome.
- Value
  - `RF`: Fitted original RF by `ranger`.
  - `fitted_Re_RF`: Results of fitted reweighted RFs for prediction.
  - `scale_ls`: Scaling factors in gaussian kernel to be used in
    prediction.
  - `optimal_eta`: A vector of optimal eta’s for all versions of
    reweighted RFs.
  - `ens_weights`: A vector of ensemble weights for all RF learners.

**`pred_ensemble`**: Predict test samples using ensemble estimate.

- Argument
  - `fit`: Fitted ensemble estimate by `ensemble`.
  - `distmat_ls`: A list of microbiome distances of both training and
    test samples.
  - `TR_ind`: Indices of training samples in `distmat_ls`.
  - `TE_data`: Test data (row: samples; column: features).
- Value
  - A matrix of predicted values for test samples (row: test samples;
    column: RF learners + Ensemble).

## Example for fitting ensemble

Below is a step-by-step example of fitting ensemble based on the R file
**`example_ensemble.R`**.

Load packages and functions:

``` r
library(ranger) # for random forest
library(castor) # for microbiome distance calculation
library(Rsolnp) # for finding ensemble weights
library(quadprog) # for finding ensemble weights
library(pROC) # for AUC calculation
library(Rcpp)
library(RcppArmadillo)
source('func.R')
```

Load “dataset” (600 simulated samples of 856 OTU counts + outcome) and
“tree” (a phylogenetic tree):

``` r
load('./Data/Simulated_data.RData')

OTUtab = dataset[, 1:856] # OTU count table
dataset[, 1:856] = dataset[, 1:856]/rowSums(dataset[, 1:856]) # convert counts to abundances
TR_ind = 1:300 # training sample indices
TR_data = dataset[TR_ind, ] # training data
TE_data_withY = dataset[- TR_ind, ] # test data with outcome
TE_data = TE_data_withY[, 1:856] # test data without outcome
```

Calculate microbiome distances of both training and test samples:

``` r
distmat_ls = UniFrac_dist(
  OTUtab = OTUtab, # OTU count table
  type = c('BC', 'Unweighted', 0, 0.5, 'Weighted'), # microbiome distances to calculate
  tree = tree)
TR_distmat_ls = lapply(distmat_ls,
                       function(x) x[TR_ind, TR_ind]) # distances for training samples
```

Fit ensemble estimate:

``` r
set.seed(123)
fit = ensemble(
  TR_data = TR_data, # training data to fit
  TR_distmat_ls = TR_distmat_ls, # list of microbiome distances for training samples
  n_CV = 5, # number of folds in cross-validataions
  eta_vec = seq(0.1, 1, 0.1), # possible values of hyper-parameter for gaussian kernel
  outcome = 'Y', # outcome name in TR_data
  outcome_type = 'binary', # continuous or binary outcome
  # below are arguments for RF
  num.trees = 1000, # number of trees in RF
  min.node.size = 10 # min.node.size in RF
  )
```

Check ensemble weights:

``` r
prmatrix(t(t(fit$ens_weights)), collab = 'Ensemble weight')
```

    ##               Ensemble weight
    ## RF_Original      5.017736e-09
    ## RF_BC            9.110720e-09
    ## RF_Unweighted    5.348847e-09
    ## RF_0             1.568756e-08
    ## RF_0.5           8.645484e-08
    ## RF_Weighted      9.999999e-01

``` r
barplot(fit$ens_weights, 
        names.arg = gsub('RF_', '', names(fit$ens_weights)),
        ylab = 'Ensemble weight',
        ylim = c(0, 1))
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

We can see that reweighted-RF using the weighted UniFrac distance has a
large weight.

Do prediction on test samples:

``` r
pred_ensemble = predict_ensemble(
  fit = fit, # ensemble fit
  distmat_ls = distmat_ls, # list of microbiome distances of training and test samples
  TR_ind = TR_ind, # indices of training samples in the above list of distance matrices
  TE_data = TE_data # test data
  )
print(round(pred_ensemble[1:5, ], 5))
```

    ##     RF_Original   RF_BC RF_Unweighted    RF_0  RF_0.5 RF_Weighted Ensemble
    ## 301     0.39728 0.40552       0.39664 0.39296 0.39980     0.39467  0.39467
    ## 302     0.45393 0.46690       0.45604 0.48763 0.47756     0.45451  0.45451
    ## 303     0.58228 0.61282       0.58588 0.61556 0.66280     0.62557  0.62557
    ## 304     0.62876 0.63973       0.63110 0.71348 0.73579     0.64858  0.64858
    ## 305     0.52879 0.56134       0.52841 0.57475 0.62102     0.58941  0.58941

Here “pred_ensemble” is matrix of predicted values for test samples
(row) using the 6 RF learners + ensemble (column).

Check classification errors:

``` r
class_error = apply(pred_ensemble, 2,
                    function(x) mean(ifelse(x >= 0.5, 1, 0) != TE_data_withY[, 'Y']))
prmatrix(t(t(round(class_error, 5))), collab = 'Classification error')
```

    ##               Classification error
    ## RF_Original                0.31667
    ## RF_BC                      0.30667
    ## RF_Unweighted              0.31667
    ## RF_0                       0.32000
    ## RF_0.5                     0.27333
    ## RF_Weighted                0.24667
    ## Ensemble                   0.24667

Note that the R file **`example_reweightedRF.R`** also includes an
example of fitting reweightRF with a given list of kernel matrices. You
can refer to that file if you are interested.
