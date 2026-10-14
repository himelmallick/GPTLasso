test_that("simulation output is PCL and RNG state is scoped", {
  set.seed(9)
  before <- .Random.seed
  simulated <- sim.gaussian.data(nsample = 12, nstudy = 1, nrep = 1,
                                 seed = 300, de.prob = rep(0.5, 3))
  after <- .Random.seed
  expect_identical(after, before)
  expect_named(simulated$x_train,
               c("feature_table", "sample_metadata", "feature_metadata"))
  expect_named(simulated$x_test,
               c("feature_table", "sample_metadata", "feature_metadata"))
  expect_true(all(c("sample_id", "study", "Y") %in%
                  colnames(simulated$x_train$sample_metadata)))
  expect_identical(colnames(simulated$x_train$feature_table),
                   simulated$x_train$sample_metadata$sample_id)
  expect_identical(rownames(simulated$x_train$feature_table),
                   simulated$x_train$feature_metadata$featureID)
})

test_that("raw generator is explicit and PCL wrappers are ready to fit", {
  raw <- gen_simmba_multistudy(
    nsample = c(12, 10), nstudy = 2, nrep = 1,
    de.prob = rep(0.5, 3), seed = 301
  )
  expect_identical(raw$nsample, c(12L, 10L))
  expect_type(raw$trainDat, "list")
  expect_type(raw$testDat, "list")
  expect_false(inherits(raw$trainDat, "MultiAssayExperiment"))

  gaussian <- sim.gaussian.data(
    nsample = c(12, 10), nstudy = 2, nrep = 1,
    de.prob = rep(0.5, 3), seed = 301
  )
  binary <- sim.binary.data(
    nsample = c(12, 10), nstudy = 2, nrep = 1,
    de.prob = rep(0.5, 3), seed = 302
  )
  for (simulated in list(gaussian, binary)) {
    expect_type(simulated$x_train, "list")
    expect_type(simulated$x_test, "list")
    expect_named(simulated$x_train,
                 c("feature_table", "sample_metadata", "feature_metadata"))
  }
})
