test_that("the packaged Gaussian README workflow runs end to end", {
  gaussian <- load_packaged_dataset("gaussianData")
  set.seed(2026)
  cv_fit <- cv.gptLasso(
    x = gaussian$x_train,
    target = NULL,
    family = "gaussian",
    type.measure = "mse",
    rho = c(0, 0.5, 1),
    alpha.ptlasso.list = seq(0, 1, 0.25),
    overall.lambda = "lambda.min",
    ind.lambda = "lambda.1se",
    pre.lambda = "lambda.1se",
    nfolds = 5,
    verbose = FALSE
  )
  prediction <- predict(
    cv_fit,
    xtest = gaussian$x_test,
    target = NULL,
    alpha.ptlasso = NULL,
    alpha.ptlasso.type = "varying",
    overall.lambda = "lambda.min",
    ind.lambda = "lambda.1se",
    pre.lambda = "lambda.1se",
    type = "response"
  )

  expect_s3_class(cv_fit, "cv.gptLasso")
  expect_true(all(vapply(cv_fit$fit, inherits, logical(1), "gptLasso")))
  expect_s3_class(prediction, "predict.cv.gptLasso")
  expect_named(prediction$yhatpre, c("Study_1", "Study_2", "Study_3"))
  expect_true(all(c("MSE", "r2") %in% names(prediction$metrics)))
  expect_identical(dim(prediction$metrics$MSE), c(3L, 4L))
  expect_identical(
    colnames(prediction$metrics$MSE),
    c("mean", "Study_1", "Study_2", "Study_3")
  )
  expect_identical(
    colnames(prediction$metrics$r2),
    c("mean", "Study_1", "Study_2", "Study_3")
  )
  expect_identical(
    colnames(cv_fit$errpre),
    c("alpha.ptlasso", "mean", "Study_1", "Study_2", "Study_3")
  )
  expect_false("pooled" %in% colnames(prediction$metrics$MSE))
  expect_false("pooled" %in% colnames(prediction$metrics$r2))
  expect_false("pooled" %in% colnames(cv_fit$errpre))
  for (name in c("erroverall", "errind", "errpre")) {
    expect_false(any(grepl("pooled", names(prediction[[name]]), fixed = TRUE)))
  }
  expect_false(any(grepl("pooled", names(cv_fit$erroverall), fixed = TRUE)))
  expect_false(any(grepl("pooled", names(cv_fit$errind), fixed = TRUE)))
  expect_true(all(is.finite(prediction$metrics$MSE)))
  expect_true(all(prediction$metrics$r2 > 0))
  expect_identical(
    rownames(prediction$metrics$MSE)[which.min(prediction$metrics$MSE[, "mean"])],
    "pre"
  )
})
