# GPTLasso

This repository houses the R package for multistudy multimodal transfer learning using **Global Pretraining and LASSO (`GPTLasso`)**.

It extends the `ptLasso` framework to support global multiview learning across multiple studies through cooperative-learning-based pretraining and transfer learning, with the current workflow centered on `gptLasso()` and `cv.gptLasso()`.

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

### Input

### Tool

``` r
fit <- gptLasso(
    x,                                                      # Input
    alpha_ptlasso = 0.5,                                    # Transfer-learning level in `[0, 1]`
    family = c("gaussian", "binomial"),                     # Response family
    type.measure = c("default", "mse", "auc", "deviance"),  # Cross-validation metric used inside the multiview fits
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
 $ feature_table   : num [1:658, 1:420] 0.0308 0.0189 0.4211 0.1075 0.33 ...
  ..- attr(*, "dimnames")=List of 2
 $ sample_metadata :'data.frame':       420 obs. of  5 variables:
 $ feature_metadata:'data.frame':       658 obs. of  2 variables:

Sample metadata columns:
[1] "subjectID" "Y"         "Xbeta"     "study"     "sample_id"

Study counts:

Study_1 Study_2 Study_3 
    140     140     140 
```

### 2. Fit a multistudy multiview transfer-learning model with `cv.gptLasso()`

``` r
cv_fit <- cv.gptLasso(
  x = x_train,                                  # Training input
  family = "gaussian",                          # Response family
  type.measure = "mse",                         # # Cross-validation metric used inside the multiview fits
  alpha_ptlasso_list = seq(0, 1, length = 11),  # Numeric vector of transfer-learning values to compare
  nfolds = 10,                                  # Cross-validation fold
  verbose = TRUE                                # Track model fitting
)
```

Inspect the selected alpha and the performance grid:

``` r
cv_fit$alpha_ptlasso_hat
cv_fit$varying.alpha_ptlasso_hat
cv_fit$errpre
```

Example output:

``` text
Gaussian cv.gptLasso() alpha grid summary:
[1] 0.4

Study_1 Study_2 Study_3 
    0.4     0.4     0.3 

      alpha_ptlasso  overall     mean group_Study_1 group_Study_2 group_Study_3
 [1,]           0.0 17.37729 17.37729      16.89311      19.25987      15.97890
 [2,]           0.1 14.14973 14.14973      16.00852      14.15400      12.28666
 [3,]           0.2 12.44223 12.44223      14.36584      12.20872      10.75211
 [4,]           0.3 12.47641 12.47641      14.33569      12.41417      10.67936
 [5,]           0.4 11.75454 11.75454      13.41357      11.05694      10.79311
 [6,]           0.5 12.43800 12.43800      14.91034      11.56325      10.84039
 [7,]           0.6 12.10708 12.10708      13.86219      11.14541      11.31363
 [8,]           0.7 12.42077 12.42077      14.40414      11.88224      10.97594
 [9,]           0.8 12.34025 12.34025      14.25577      11.71032      11.05467
[10,]           0.9 13.25401 13.25401      14.14932      12.04975      13.56297
[11,]           1.0 12.70732 12.70732      14.93202      12.24143      10.94851
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

Inspect the prediction object and held-out performance summary:

``` r
names(pred)
pred$errpre
```

Example output:

``` text
[1] "call" "alpha_ptlasso" "yhatoverall" "yhatind" "yhatpre"
[6] "supoverall" "supind" "suppre.common" "suppre.individual"
[10] "type.measure" "erroverall" "errind" "errpre"

Pretrained performance summary:
    allGroups          mean group_Study_1 group_Study_2 group_Study_3 
   13.4524277    13.4524277    12.7809198    19.1572324     8.4191310 
          r^2 
    0.6839485 
```

## Citation

If you use this repository, please cite it as:

``` text
Gao C, Mallick H (2026). Multistudy Multimodal Pretraining and Transfer Learning.
Research abstract and open-source software for multistudy, multimodal transfer learning.
```

## Issues

For bugs, questions, or feature requests: contact information to be added.
