# GPTLasso

This repository houses the `GPTLasso` R package for multistudy multimodal transfer learning using pretraining and LASSO.

It extends the `ptLasso` framework to support global multiview learning across multiple studies through cooperative-learning-based pretraining and transfer learning, with the current workflow centered on `gptLasso()` and `cv.gptLasso()`.

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

## Quick Example

The shortest workflow is to simulate a multistudy multiview dataset, fit a transfer-learning model with `gptLasso()`, and then tune the transfer level with `cv.gptLasso()`.

``` r
set.seed(1234)

sim_dat <- sim.gaussian.data()
x_train <- sim_dat$x_train
x_test <- sim_dat$x_test

fit <- gptLasso(
  x = x_train,
  alpha_ptlasso = 0.5,
  family = "gaussian",
  type.measure = "mse",
  nfolds = 3
)

cv_fit <- cv.gptLasso(
  x = x_train,
  family = "gaussian",
  type.measure = "mse",
  alpha_ptlasso_list = c(0, 0.5, 1),
  nfolds = 3
)
```

## Tutorial

This Gaussian example tutorial walks through simulation, model fitting, transfer-level tuning, and prediction on held-out data.

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
List of 3
 $ feature_table   : num [1:497, 1:420] ...
 $ sample_metadata : 'data.frame': 420 obs. of ...
 $ feature_metadata: 'data.frame': 497 obs. of ...

[1] "subjectID" "Y" "Xbeta" "study" "sample_id"

Study_1 Study_2 Study_3
    140     140     140
```

### 2. Fit a multistudy multiview transfer-learning model

``` r
fit <- gptLasso(
  x = x_train,
  alpha_ptlasso = 0.5,
  family = "gaussian",
  type.measure = "mse",
  nfolds = 3
)
```

Inspect the main returned components:

``` r
names(fit)
fit$study_names
fit$view_names
fit$n_by_study
```

Example output:

``` text
[1] "call" "k" "N_all" "n_by_study" "study_names" "group.levels"
[7] "view_names" "n_views" "p_by_view" "features_all" "alpha_ptlasso" ...

[1] "Study_1" "Study_2" "Study_3"
[1] "methyl" "expr" "protein"

Study_1 Study_2 Study_3
    140     140     140
```

### 3. Tune the transfer-learning level with `cv.gptLasso()`

``` r
cv_fit <- cv.gptLasso(
  x = x_train,
  family = "gaussian",
  type.measure = "mse",
  alpha_ptlasso_list = c(0, 0.5, 1),
  nfolds = 3
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
[1] 0.5

 Study_1 Study_2 Study_3
     0.5     0.5     0.5

     alpha_ptlasso   overall      mean group_Study_1 group_Study_2 group_Study_3
[1,]            0  ...
[2,]          0.5  ...
[3,]            1  ...
```

### 4. Predict on held-out data

``` r
pred <- predict(
  fit,
  xtest = x_test,
  ytest = y_test,
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

 allGroups       mean group_Study_1 group_Study_2 group_Study_3         r^2
      ...        ...          ...          ...          ...          ...
```

## Citation

If you use this repository, please cite it as:

``` text
Gao C, Mallick H (2026). Multistudy Multimodal Pretraining and Transfer Learning.
Research abstract and open-source software for multistudy, multimodal transfer learning.
```

## Issues

For bugs, questions, or feature requests: contact information to be added.
