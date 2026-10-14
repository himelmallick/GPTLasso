test_that("packaged simulations contain aligned train and test PCL objects", {
  for (name in c("gaussianData", "binaryData", "survivalData")) {
    object <- load_packaged_dataset(name)
    expect_named(object, c("x_train", "x_test", "true_betas", "settings"))
    for (pcl in list(object$x_train, object$x_test)) {
      expect_named(pcl, c("feature_table", "sample_metadata", "feature_metadata"))
      expect_identical(colnames(pcl$feature_table), pcl$sample_metadata$sample_id)
      expect_identical(rownames(pcl$feature_table), pcl$feature_metadata$featureID)
      expect_identical(unique(as.character(pcl$feature_metadata$featureType)),
                       c("methyl", "expr", "protein"))
      expect_identical(unique(as.character(pcl$sample_metadata$study)),
                       c("Study_1", "Study_2", "Study_3"))
    }
    expect_identical(object$x_train$feature_metadata,
                     object$x_test$feature_metadata)
  }
})

test_that("packaged outcomes satisfy their family contracts", {
  gaussian <- load_packaged_dataset("gaussianData")
  expect_type(gaussian$x_train$sample_metadata$Y, "double")
  expect_type(gaussian$x_test$sample_metadata$Y, "double")

  binary <- load_packaged_dataset("binaryData")
  for (pcl in list(binary$x_train, binary$x_test)) {
    metadata <- pcl$sample_metadata
    by_study <- split(metadata$Y, metadata$study)
    expect_true(all(vapply(by_study, function(y) {
      identical(sort(unique(y)), c(0L, 1L))
    }, logical(1))))
  }

  survival <- load_packaged_dataset("survivalData")
  for (pcl in list(survival$x_train, survival$x_test)) {
    metadata <- pcl$sample_metadata
    expect_s3_class(metadata$Y, "Surv")
    event_counts <- tapply(metadata$Y[, 2], metadata$study, sum)
    expect_true(all(event_counts >= 3L))
  }
})

test_that("packaged simulation settings record deterministic seeds", {
  expected <- list(
    gaussianData = list(
      seed = 1234, nsample = c(300, 200, 100), rho.beta = 0.75,
      train = c(Study_1 = 210L, Study_2 = 140L, Study_3 = 70L),
      test = c(Study_1 = 90L, Study_2 = 60L, Study_3 = 30L)
    ),
    binaryData = list(
      seed = 5678, nsample = c(60, 50, 40), rho.beta = 0.5,
      train = c(Study_1 = 42L, Study_2 = 35L, Study_3 = 28L),
      test = c(Study_1 = 18L, Study_2 = 15L, Study_3 = 12L)
    ),
    survivalData = list(
      seed = 9012, nsample = c(60, 50, 40), rho.beta = 0.5,
      train = c(Study_1 = 42L, Study_2 = 35L, Study_3 = 28L),
      test = c(Study_1 = 18L, Study_2 = 15L, Study_3 = 12L)
    )
  )

  for (name in names(expected)) {
    object <- load_packaged_dataset(name)
    spec <- expected[[name]]
    expect_equal(object$settings$seed, spec$seed)
    expect_equal(object$settings$nsample, spec$nsample)
    expect_equal(object$settings$rho.beta, spec$rho.beta)
    train_counts <- table(object$x_train$sample_metadata$study)
    test_counts <- table(object$x_test$sample_metadata$study)
    expect_identical(names(train_counts), names(spec$train))
    expect_identical(names(test_counts), names(spec$test))
    expect_equal(as.integer(train_counts), unname(spec$train))
    expect_equal(as.integer(test_counts), unname(spec$test))
    expect_equal(as.integer(train_counts + test_counts), spec$nsample)
  }
})
