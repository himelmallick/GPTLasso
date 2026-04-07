# GPTLasso

This repository houses the R package for multistudy multimodal transfer learning using **Global Pretraining and LASSO (`GPTLasso`)**.

It extends the `ptLasso` framework to support global multiview learning across multiple studies through cooperative-learning-based pretraining and transfer learning, with the current workflow centered on `gptLasso()` and `cv.gptLasso()`.

![](figures/workflow.png)

## Background

Modern biomedical prediction problems often involve multiple data modalities collected across several related cohorts, while each individual study may still be too small for stable model fitting. This project addresses that setting by combining multiview modeling with transfer learning so that information can be shared across both studies and modalities while still allowing study-specific refinement.

The framework uses a pretrained multiview model to learn shared structure, then fine-tunes study-level fits for improved local prediction. In simulations and motivating multi-omics applications, the goal is to improve predictive accuracy, estimation quality, and cross-study generalization relative to single-study or single-view approaches.

**Keywords:** Transfer learning, Cooperative learning, Multistudy analysis, Multimodal Integration, LASSO, Pretraining

## Installation

You can install the development version directly from GitHub:

``` r
install.packages("devtools")
devtools::install_github("himelmallick/GPTLasso")
library(GPTLasso)
```

## Get Started

### Bioconductor-style Input

GPTLasso requires a `x` as a named list containing `feature_table`, `sample_metadata`, and `feature_metadata`.

![](figures/Input.png)

**Note on `study` Labels:**

`gptLasso` detects the study levels directly from `x$sample_metadata$study`. The unique values in that column, in their observed order, define the groups used to split the training data, constract study-specific models, and organize prediction summaries.

In this repository, `study` is the motivating example for handling different groups of data. In practice, users are welcome to define their own grouping variable in the same column, such as `area`, site, cohort, or any other project-tailored grouping that serves the scientific or operational goal of the analysis.

### Tool

``` r
fit <- gptLasso(
    x,                                                      # Input
    alpha_ptlasso = 0.5,                                    # Transfer-learning level in `[0, 1]`
    family = c("gaussian", "binomial"),                     # Response family
    type.measure = c("default", "mse", "auc", "deviance"),  # Cross-validation metric used inside the multiview fits
    rho = c(0, 0.1, 0.25, 0.5, 1, 5, 10),                   # Multiview cooperative learning confusion stage parameter, e.g., rho=0 -> early fusion
    overall.lambda = c("lambda.1se", "lambda.min"),         # Lambda rule used for the stage-one overall model
    ind.lambda = c("lambda.1se", "lambda.min"),             # Lambda rule used for the individual models
    pre.lambda = c("lambda.1se", "lambda.min"),             # Lambda rule used for the pretrained models
    nfolds = 10,                                            # Number of folds used when `foldid` is not supplied
    alpha_glmnet = 1                                        # Elastic-net mixing parameter, default is 1 indicating Lasso regression
)
```

## Tutorial

This Gaussian example tutorial walks through simulated data generation, model fitting, transfer-level tuning, and prediction on held-out data.

### 1. Simulate Gaussian multistudy multiview data

``` r
set.seed(1234)

gaussian_sim <- sim.gaussian.data()
x_train <- gaussian_sim$x_train
x_test <- gaussian_sim$x_test
y_test <- split(x_test$sample_metadata$Y, x_test$sample_metadata$study)
```

Inspect the container layout before fitting:

``` r
str(x_train, max.level = 1)
colnames(x_train$sample_metadata)
table(x_train$sample_metadata$study)
```

Example output:

``` text
==============================
Gaussian x_train 
==============================
List of 3
 $ feature_table   : num [1:658, 1:840] 0.0261 0.0215 0.7405 0.0733 0.2195 ...
  ..- attr(*, "dimnames")=List of 2
 $ sample_metadata :'data.frame':       840 obs. of  5 variables:
 $ feature_metadata:'data.frame':       658 obs. of  2 variables:

Sample metadata columns:
[1] "subjectID" "Y"         "Xbeta"     "study"     "sample_id"

Study counts:

Study_1 Study_2 Study_3 
    350     280     210 
```

### 2. Fit a multistudy multiview transfer-learning model with `cv.gptLasso()`

``` r
cv_fit <- cv.gptLasso(
  x = x_train,                                  # Training input
  family = "gaussian",                          # Response family
  type.measure = "mse",                         # Cross-validation metric used inside the multiview fits
  rho = c(0, 0.5, 1),                           # Multiview cooperative learning fusion stage parameter
  alpha_ptlasso_list = seq(0, 1, length = 11),  # Numeric vector of transfer-learning values to compare
  nfolds = 10,                                  # Cross-validation fold
  verbose = TRUE                                # Track model fitting
)
```

Inspect the selected alpha and the performance grid:

``` r
cv_fit$alpha_ptlasso_hat          # selected fixed transfer-learning value
cv_fit$varying.alpha_ptlasso_hat  # study-specific transfer-learning values chosen from the same grid
cv_fit$errpre                     # performance of pretrained model for each candidate alpha
cv_fit$fitpre.rho                 # study-specific fusion stage value of pretrained model
```

Example output:

``` text
Gaussian cv.gptLasso() alpha grid summary:
[1] 0.8

Study_1 Study_2 Study_3 
    1.0     0.6     0.6 

      alpha_ptlasso  overall     mean group_Study_1 group_Study_2 group_Study_3
 [1,]           0.0 13.92094 14.14564     12.675854      14.38881      15.37226
 [2,]           0.1 12.75010 13.05121     11.217638      13.10503      14.83095
 [3,]           0.2 12.25312 12.40972     11.222819      12.90437      13.10197
 [4,]           0.3 11.40037 11.59173     10.283390      11.91201      12.57980
 [5,]           0.4 10.97536 11.18962      9.780630      11.43652      12.35171
 [6,]           0.5 10.52049 10.69164      9.374715      11.27168      11.42853
 [7,]           0.6 10.47295 10.60504      9.702109      10.82588      11.28712
 [8,]           0.7 10.66208 10.86011      9.551828      11.10027      11.92824
 [9,]           0.8 10.32312 10.57686      8.858216      10.96920      11.90317
[10,]           0.9 10.90699 11.23096      8.999471      11.80635      12.88705
[11,]           1.0 11.16569 11.56169      8.703100      12.52686      13.45512
```

### 4. Predict on held-out data

``` r
pred <- predict(
  cv_fit,
  xtest = x_test,
  ytest = y_test,
  alpha_ptlasso =  NULL,          # Optional user-specified transfer-learning choice. May be one value or one per study
  alpha_ptlasso_type = "varying", # Either `"fixed"` or `"varying", when `alpha_ptlasso` is not supplied
  type = "response"
)
```

Inspect the prediction object and compare held-out performance across all models from `metrics`:

``` r
names(pred)
pred$metrics$MSE
pred$metrics$r2
```

Example output:

``` text
 [1] "call"              "alpha_ptlasso"     "yhatoverall"      
 [4] "yhatind"           "yhatpre"           "supoverall"       
 [7] "supind"            "suppre.common"     "suppre.individual"
[10] "type.measure"      "metrics"           "erroverall"       
[13] "errind"            "errpre"            "fit"

         overall ind_group_Study_1 ind_group_Study_2 ind_group_Study_3 
        16.572974          8.939839         13.491641         11.968390 
   ind_group_mean pre_group_Study_1 pre_group_Study_2 pre_group_Study_3 
        11.466623          8.939839         12.432750          9.877647 
   pre_group_mean 
        10.416746 
        
                  overall ind_group_Study_1 ind_group_Study_2 ind_group_Study_3 
        0.5358854         0.7405784         0.6517837         0.6101880 
   ind_group_mean pre_group_Study_1 pre_group_Study_2 pre_group_Study_3 
        0.6675167         0.7405784         0.6791134         0.6782837 
   pre_group_mean 
        0.6993252 
```

## Citation

If you use this repository, please cite it as:

``` text
Gao C, Mallick H (2026). Multistudy Multimodal Pretraining and Transfer Learning.
Research abstract and open-source software for multistudy, multimodal transfer learning.
```

## Issues

For bugs, questions, or feature requests: contact information to be added.
