make_test_mae <- function(n_study = 2L, n_view = 2L,
                          family = "gaussian", n_each = 12L,
                          assays_per_view = 1L) {
  set.seed(410)
  n <- n_study * n_each
  primary <- paste0("primary", seq_len(n))
  studies <- rep(paste0("study", seq_len(n_study)), each = n_each)
  experiment_list <- lapply(seq_len(n_view), function(v) {
    matrix_value <- matrix(rnorm(4L * n), nrow = 4L,
                           dimnames = list(paste0("feature", seq_len(4L)),
                                           paste0("view", v, "_column", seq_len(n))))
    assay_list <- list(value = matrix_value)
    if (assays_per_view > 1L) assay_list$alternative <- matrix_value + 1
    SummarizedExperiment::SummarizedExperiment(assays = assay_list)
  })
  names(experiment_list) <- paste0("view", seq_len(n_view))
  signal <- colSums(SummarizedExperiment::assay(experiment_list[[1]], "value")[1:2, , drop = FALSE])
  response <- switch(
    family,
    gaussian = signal + rnorm(n, sd = 0.5),
    binomial = as.integer(signal > stats::median(signal)),
    cox = I(survival::Surv(
      time = seq_len(n) / exp(as.numeric(scale(signal))),
      event = rep(c(1, 1, 0), length.out = n)
    ))
  )
  map <- do.call(rbind, lapply(names(experiment_list), function(view) {
    S4Vectors::DataFrame(
      assay = factor(view, levels = names(experiment_list)),
      primary = primary,
      colname = colnames(experiment_list[[view]])
    )
  }))
  MultiAssayExperiment::MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(experiment_list),
    colData = S4Vectors::DataFrame(
      study = studies, Y = response, row.names = primary
    ),
    sampleMap = map
  )
}

load_packaged_dataset <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "GPTLasso", envir = env)
  env[[name]]
}
