test_that("single-view multistudy fit keeps transfer and fixes rho", {
  mae <- make_test_mae(2, 1)
  fit <- gptLasso(mae, alpha.ptlasso = 0.5, rho = c(0, 1), nfolds = 3)
  expect_s3_class(fit, "gptLasso")
  expect_true(fit$transfer.applicable)
  expect_false(fit$fusion.applicable)
  expect_equal(fit$alpha.ptlasso, 0.5)
  expect_equal(fit$rho, 0)
  expect_equal(fit$rho.requested, c(0, 1))
  expect_false(requireNamespace("ptLasso", quietly = TRUE))
  prediction <- predict(fit, mae)
  expect_s3_class(prediction, "predict.gptLasso")
  expect_named(prediction$yhatpre, fit$study.names)
})

test_that("single-study multiview fit fixes alpha and tunes rho", {
  mae <- make_test_mae(1, 2, n_each = 18)
  fit <- gptLasso(mae, alpha.ptlasso = 0.25, rho = c(0, 1), nfolds = 3)
  expect_s3_class(fit, "gptLasso")
  expect_false(fit$transfer.applicable)
  expect_true(fit$fusion.applicable)
  expect_equal(fit$alpha.ptlasso, 1)
  expect_equal(fit$alpha.ptlasso.requested, 0.25)
  expect_identical(fit$fitind[[1]], fit$fitoverall)
  expect_identical(fit$fitpre[[1]], fit$fitoverall)
})

test_that("single-study single-view returns stable fit and prediction classes", {
  mae <- make_test_mae(1, 1, n_each = 18)
  fit <- gptLasso(mae, alpha.ptlasso = 0.2, rho = c(0, 1), nfolds = 3)
  expect_s3_class(fit, "gptLasso")
  expect_equal(fit$alpha.ptlasso, 1)
  expect_equal(fit$rho, 0)
  expect_s3_class(predict(fit, mae), "predict.gptLasso")
})

test_that("multistudy multiview Gaussian fit and target prediction work", {
  mae <- make_test_mae(2, 2)
  fit <- gptLasso(mae, alpha.ptlasso = 0.5, rho = c(0, 1),
                  target = "study2", nfolds = 3)
  expect_s3_class(fit, "gptLasso")
  expect_identical(fit$target.study.names, "study2")
  expect_named(fit$fitpre, "study2")
  expect_error(predict(fit, mae), "non-target")
  pred <- predict(fit, mae[, SummarizedExperiment::colData(mae)$study == "study2"])
  expect_true("MSE" %in% names(pred$metrics))
})

test_that("single-study cross-validation bypasses alpha search", {
  mae <- make_test_mae(1, 2, n_each = 18)
  fit <- cv.gptLasso(mae, alpha.ptlasso.list = c(0, 0.5, 1),
                     rho = c(0, 1), nfolds = 3)
  expect_s3_class(fit, "cv.gptLasso")
  expect_equal(fit$alpha.ptlasso.list, 1)
  expect_equal(fit$alpha.ptlasso.requested, c(0, 0.5, 1))
  expect_equal(fit$alpha.ptlasso.hat, 1)
  expect_equal(unname(fit$varying.alpha.ptlasso.hat), 1)
  expect_s3_class(predict(fit, mae), "predict.cv.gptLasso")
})

test_that("multistudy CV supports fixed and study-varying prediction", {
  mae <- make_test_mae(2, 1)
  fit <- cv.gptLasso(
    mae, alpha.ptlasso.list = c(0, 1), rho = c(0, 1), nfolds = 3
  )
  expect_true(fit$transfer.applicable)
  expect_equal(nrow(fit$errpre), 2)
  fixed <- predict(fit, mae, alpha.ptlasso.type = "fixed")
  varying <- predict(
    fit, mae, alpha.ptlasso = c(study1 = 0, study2 = 1),
    alpha.ptlasso.type = "varying", return.link = TRUE
  )
  expect_s3_class(fixed, "predict.cv.gptLasso")
  expect_s3_class(varying, "predict.cv.gptLasso")
  expect_named(varying$yhatpre, c("study1", "study2"))
  expect_named(varying$linkpre, c("study1", "study2"))
  expect_true("MSE" %in% names(varying$metrics))
})

test_that("binomial and Cox families retain metrics and prediction shapes", {
  binary <- make_test_mae(2, 1, family = "binomial", n_each = 15)
  binary_fit <- gptLasso(binary, family = "binomial", type.measure = "deviance",
                         rho = 0, nfolds = 3)
  binary_pred <- predict(binary_fit, binary, type = "response")
  expect_true("deviance" %in% names(binary_pred$metrics))

  cox <- make_test_mae(2, 1, family = "cox", n_each = 15)
  cox_fit <- suppressWarnings(gptLasso(
    cox, family = "cox", type.measure = "deviance", rho = 0, nfolds = 3
  ))
  cox_pred <- predict(cox_fit, cox, type = "link")
  expect_true("deviance" %in% names(cox_pred$metrics))
  expect_equal(length(cox_pred$yhatpre), 2)
})
