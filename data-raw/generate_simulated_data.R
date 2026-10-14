# Rebuild the packaged Gaussian, binary, and survival example datasets.
# Run from the package root with:
# Rscript data-raw/generate_simulated_data.R

r_files <- list.files("R", pattern = "[.][Rr]$", full.names = TRUE)
invisible(lapply(r_files, source))

dir.create("data", showWarnings = FALSE)

common_settings <- list(
  nstudy = 3,
  snr = 5,
  p.train = 0.7,
  sigma.alpha = 0.1,
  tau.snr = 0.1,
  nrep = 1,
  ygen.mode = "LM",
  de.prob = rep(0.1, 3),
  de.downProb = rep(0.5, 3),
  de.facLoc = rep(1, 3),
  de.facScale = rep(0.4, 3)
)

gaussian_settings <- c(
  list(nsample = c(300, 200, 100), rho.beta = 0.75),
  common_settings
)

small_example_settings <- c(
  list(nsample = c(60, 50, 40), rho.beta = 0.5),
  common_settings
)

build_dataset <- function(outcome.type, seed, settings) {
  raw <- do.call(
    gen_simmba_multistudy,
    c(settings, list(outcome.type = outcome.type, seed = seed))
  )
  list(
    x_train = ptmv_build_sim_container(raw, "Rep_1", "trainDat"),
    x_test = ptmv_build_sim_container(raw, "Rep_1", "testDat"),
    true_betas = raw$true_betas[["Rep_1"]],
    settings = raw[setdiff(
      names(raw), c("trainDat", "testDat", "true_betas")
    )]
  )
}

dataset_specs <- list(
  gaussianData = list(
    outcome.type = "continuous", seed = 1234,
    settings = gaussian_settings
  ),
  binaryData = list(
    outcome.type = "binary", seed = 5678,
    settings = small_example_settings
  ),
  survivalData = list(
    outcome.type = "survival", seed = 9012,
    settings = small_example_settings
  )
)

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) {
  requested <- names(dataset_specs)
}
unknown <- setdiff(requested, names(dataset_specs))
if (length(unknown) > 0L) {
  stop("Unknown dataset name(s): ", paste(unknown, collapse = ", "))
}

for (dataset_name in requested) {
  spec <- dataset_specs[[dataset_name]]
  dataset <- build_dataset(spec$outcome.type, spec$seed, spec$settings)
  assign(dataset_name, dataset)
  save(
    list = dataset_name,
    file = file.path("data", paste0(dataset_name, ".rda")),
    compress = "xz",
    version = 3
  )
}
