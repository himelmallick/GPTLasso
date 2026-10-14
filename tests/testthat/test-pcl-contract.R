test_that("all four PCL topologies normalize", {
  for (n_study in 1:2) {
    for (n_view in 1:2) {
      pcl <- make_test_pcl(n_study, n_view)
      normalized <- GPTLasso:::ptmv_normalize_container(pcl)
      expect_identical(normalized$container_type, "pcl")
      expect_equal(normalized$k, n_study)
      expect_equal(normalized$n_views, n_view)
      expect_equal(sum(normalized$n_by_study), n_study * 12L)
      expect_null(normalized$assay_names)
    }
  }
})

test_that("equivalent PCL and MAE normalize identically", {
  pcl <- make_test_pcl()
  mae <- pcl_to_mae(pcl)
  pcl_norm <- GPTLasso:::ptmv_normalize_container(pcl, family = "gaussian")
  mae_norm <- GPTLasso:::ptmv_normalize_container(mae, family = "gaussian")
  for (field in c("x", "y", "study_names", "view_names", "n_by_study",
                  "p_by_view", "feature_ids", "features_by_view")) {
    expect_equal(pcl_norm[[field]], mae_norm[[field]])
  }
})

test_that("malformed PCL inputs fail clearly", {
  pcl <- make_test_pcl()
  expect_error(GPTLasso:::ptmv_normalize_container(list()), "non-empty list")

  missing_y <- pcl
  missing_y$sample_metadata$Y <- NULL
  expect_error(GPTLasso:::ptmv_normalize_container(missing_y), "must contain")

  duplicate_sample <- pcl
  duplicate_sample$sample_metadata$sample_id[2] <- duplicate_sample$sample_metadata$sample_id[1]
  expect_error(GPTLasso:::ptmv_normalize_container(duplicate_sample), "must be unique")

  empty_sample <- pcl
  empty_sample$sample_metadata$sample_id[1] <- ""
  expect_error(GPTLasso:::ptmv_normalize_container(empty_sample), "non-missing and non-empty")

  duplicate_feature <- pcl
  duplicate_feature$feature_metadata$featureID[2] <- duplicate_feature$feature_metadata$featureID[1]
  expect_error(GPTLasso:::ptmv_normalize_container(duplicate_feature), "must be unique")

  empty_feature <- pcl
  empty_feature$feature_metadata$featureID[1] <- ""
  expect_error(GPTLasso:::ptmv_normalize_container(empty_feature), "non-missing and non-empty")

  bad_sample_order <- pcl
  bad_sample_order$sample_metadata <- bad_sample_order$sample_metadata[nrow(bad_sample_order$sample_metadata):1, ]
  expect_error(GPTLasso:::ptmv_normalize_container(bad_sample_order), "Sample order mismatch")

  bad_feature_order <- pcl
  bad_feature_order$feature_metadata <- bad_feature_order$feature_metadata[nrow(bad_feature_order$feature_metadata):1, ]
  expect_error(GPTLasso:::ptmv_normalize_container(bad_feature_order), "Feature order mismatch")

  missing_study <- pcl
  missing_study$sample_metadata$study[1] <- NA_character_
  expect_error(GPTLasso:::ptmv_normalize_container(missing_study), "study must be non-missing")

  missing_view <- pcl
  missing_view$feature_metadata$featureType[1] <- ""
  expect_error(GPTLasso:::ptmv_normalize_container(missing_view), "non-missing and non-empty")

  expect_error(
    GPTLasso:::ptmv_normalize_container(pcl, assay.names = "value"),
    "applies only to MultiAssayExperiment"
  )
})

test_that("PCL preserves Cox responses", {
  pcl <- make_test_pcl(family = "cox")
  normalized <- GPTLasso:::ptmv_normalize_container(pcl, family = "cox")
  expect_s3_class(normalized$y_all, "Surv")
  expect_true(all(vapply(normalized$y, inherits, logical(1), "Surv")))

  invalid <- pcl
  invalid$sample_metadata$Y <- seq_len(nrow(invalid$sample_metadata))
  expect_error(
    GPTLasso:::ptmv_normalize_container(invalid, family = "cox"),
    "Surv"
  )
})

test_that("fit and prediction allow all PCL and MAE combinations", {
  pcl <- make_test_pcl(n_study = 1, n_view = 1, n_each = 18)
  mae <- pcl_to_mae(pcl)
  foldid <- rep(seq_len(3), length.out = 18)
  set.seed(991)
  fit_pcl <- gptLasso(pcl, rho = 0, nfolds = 3, foldid = foldid)
  set.seed(991)
  fit_mae <- gptLasso(mae, rho = 0, nfolds = 3, foldid = foldid)

  predictions <- list(
    pcl_pcl = predict(fit_pcl, pcl),
    pcl_mae = predict(fit_pcl, mae),
    mae_pcl = predict(fit_mae, pcl),
    mae_mae = predict(fit_mae, mae)
  )
  expect_true(all(vapply(predictions, inherits, logical(1), "predict.gptLasso")))
  expect_equal(predictions$pcl_pcl$yhatoverall, predictions$pcl_mae$yhatoverall)
  expect_equal(predictions$mae_pcl$yhatoverall, predictions$mae_mae$yhatoverall)
  expect_equal(predictions$pcl_pcl$yhatoverall, predictions$mae_mae$yhatoverall)
})
