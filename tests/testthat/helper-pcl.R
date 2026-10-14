make_test_pcl <- function(n_study = 2L, n_view = 2L,
                          family = "gaussian", n_each = 12L) {
  set.seed(411)
  n <- n_study * n_each
  sample_id <- paste0("sample", seq_len(n))
  study <- rep(paste0("study", seq_len(n_study)), each = n_each)
  view <- rep(paste0("view", seq_len(n_view)), each = 4L)
  feature_id <- paste0(view, "_feature", rep(seq_len(4L), n_view))
  feature_table <- matrix(
    rnorm(length(feature_id) * n), nrow = length(feature_id),
    dimnames = list(feature_id, sample_id)
  )
  signal <- colSums(feature_table[seq_len(2L), , drop = FALSE])
  response <- switch(
    family,
    gaussian = signal + rnorm(n, sd = 0.5),
    binomial = as.integer(signal > stats::median(signal)),
    cox = I(survival::Surv(
      time = seq_len(n) / exp(as.numeric(scale(signal))),
      event = rep(c(1, 1, 0), length.out = n)
    ))
  )
  sample_metadata <- data.frame(
    sample_id = sample_id, study = study, stringsAsFactors = FALSE,
    row.names = sample_id
  )
  sample_metadata$Y <- response
  feature_metadata <- data.frame(
    featureID = feature_id, featureType = view, stringsAsFactors = FALSE,
    row.names = feature_id
  )
  list(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  )
}

pcl_to_mae <- function(pcl, assay_name = "value") {
  views <- unique(as.character(pcl$feature_metadata$featureType))
  experiments <- lapply(views, function(view) {
    keep <- pcl$feature_metadata$featureType == view
    SummarizedExperiment::SummarizedExperiment(
      assays = stats::setNames(list(pcl$feature_table[keep, , drop = FALSE]), assay_name)
    )
  })
  names(experiments) <- views
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiments),
    colData = S4Vectors::DataFrame(pcl$sample_metadata)
  )
}
