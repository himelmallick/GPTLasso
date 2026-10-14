test_that("all four MAE topologies normalize", {
  for (n_study in 1:2) {
    for (n_view in 1:2) {
      mae <- make_test_mae(n_study, n_view)
      normalized <- GPTLasso:::ptmv_normalize_container(mae)
      expect_equal(normalized$k, n_study)
      expect_equal(normalized$n_views, n_view)
      expect_equal(sum(normalized$n_by_study), n_study * 12L)
    }
  }
})

test_that("sampleMap aligns distinct assay and primary identifiers", {
  mae <- make_test_mae()
  normalized <- GPTLasso:::ptmv_normalize_container(mae)
  expect_identical(rownames(normalized$x$study1$view1),
                   paste0("primary", 1:12))
  expect_identical(colnames(normalized$x$study1$view1),
                   paste0("feature", 1:4))
})

test_that("assay selection is explicit and checked", {
  mae <- make_test_mae(assays_per_view = 2L)
  expect_error(GPTLasso:::ptmv_normalize_container(mae),
               "exactly one assay")
  selected <- GPTLasso:::ptmv_normalize_container(mae, assay.names = "value")
  expect_identical(unname(selected$assay_names), rep("value", 2))
  selected2 <- GPTLasso:::ptmv_normalize_container(
    mae, assay.names = c(view2 = "alternative", view1 = "value")
  )
  expect_identical(names(selected2$assay_names), c("view1", "view2"))
  expect_error(GPTLasso:::ptmv_normalize_container(mae, assay.names = "absent"),
               "absent")
  expect_error(
    GPTLasso:::ptmv_normalize_container(
      mae, assay.names = c(view1 = "value", wrong_view = "value")
    ),
    "names exactly match"
  )
})

test_that("invalid MAE layouts fail clearly", {
  expect_error(GPTLasso:::ptmv_normalize_container(list()),
               "non-empty list")
  expect_error(GPTLasso:::ptmv_normalize_container(
    MultiAssayExperiment::MultiAssayExperiment()), "at least one experiment")

  mae <- make_test_mae()
  no_study <- mae
  SummarizedExperiment::colData(no_study)$study <- NULL
  expect_error(GPTLasso:::ptmv_normalize_container(no_study), "must contain")
  no_y <- mae
  SummarizedExperiment::colData(no_y)$Y <- NULL
  expect_error(GPTLasso:::ptmv_normalize_container(no_y), "must contain")

  partial_experiments <- MultiAssayExperiment::experiments(mae)
  partial_experiments[[1]] <- partial_experiments[[1]][, -1]
  partial_map <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  partial_map <- partial_map[!(partial_map$assay == "view1" &
                               partial_map$primary == "primary1"), ]
  partial <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = partial_experiments,
    colData = SummarizedExperiment::colData(mae),
    sampleMap = S4Vectors::DataFrame(partial_map)
  )
  expect_error(GPTLasso:::ptmv_normalize_container(partial),
               "complete one-to-one")

  empty_feature <- mae
  empty_experiments <- MultiAssayExperiment::experiments(empty_feature)
  rownames(empty_experiments[[1]])[1] <- ""
  MultiAssayExperiment::experiments(empty_feature) <- empty_experiments
  expect_error(GPTLasso:::ptmv_normalize_container(empty_feature),
               "unique, non-empty feature")

  duplicated_feature <- mae
  duplicated_experiments <- MultiAssayExperiment::experiments(duplicated_feature)
  rownames(duplicated_experiments[[1]])[2] <- rownames(duplicated_experiments[[1]])[1]
  MultiAssayExperiment::experiments(duplicated_feature) <- duplicated_experiments
  expect_error(GPTLasso:::ptmv_normalize_container(duplicated_feature),
               "unique, non-empty feature")
})

test_that("new data must preserve within-view feature layouts", {
  mae <- make_test_mae()
  template <- GPTLasso:::ptmv_normalize_container(mae, family = "gaussian")

  one_view_map <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  one_view_map <- one_view_map[one_view_map$assay == "view1", ]
  one_view_map$assay <- factor(one_view_map$assay, levels = "view1")
  missing_view <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::experiments(mae)[1],
    colData = SummarizedExperiment::colData(mae),
    sampleMap = S4Vectors::DataFrame(one_view_map)
  )
  expect_error(GPTLasso:::ptmv_normalize_newdata(missing_view, template),
               "same views")

  changed <- mae
  rownames(MultiAssayExperiment::experiments(changed)[[1]])[1] <- "changed"
  expect_error(GPTLasso:::ptmv_normalize_newdata(changed, template),
               "feature identifiers")

  study <- SummarizedExperiment::colData(mae)$study
  reordered <- mae[, c(which(study == "study2"), which(study == "study1"))]
  expect_error(
    GPTLasso:::ptmv_normalize_newdata(reordered, template),
    "xtest studies.*ordering"
  )

  unknown_study <- mae
  SummarizedExperiment::colData(unknown_study)$study[
    SummarizedExperiment::colData(unknown_study)$study == "study2"
  ] <- "unknown"
  expect_error(
    GPTLasso:::ptmv_normalize_newdata(unknown_study, template),
    "xtest study labels.*subset"
  )
})

test_that("replicated primary mappings are rejected", {
  mae <- make_test_mae()
  map <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
  duplicate <- map[map$assay == "view1" & map$primary == "primary1", ]
  duplicate$colname <- "view1_column2"
  replicated <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::experiments(mae),
    colData = SummarizedExperiment::colData(mae),
    sampleMap = S4Vectors::DataFrame(rbind(map, duplicate))
  )
  expect_error(GPTLasso:::ptmv_normalize_container(replicated),
               "complete one-to-one")
})

test_that("unsupported containers fail with the dual-format contract", {
  expect_error(gptLasso(1, nfolds = 2), "MultiAssayExperiment or a PCL")
  expect_error(gptLasso(list(), nfolds = 2), "non-empty list")
})
