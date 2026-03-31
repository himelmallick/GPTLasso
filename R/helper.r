# Helper: simulation base function
.trigger_InterSIM_safe <- function(n) {
  sim.data <- InterSIM::InterSIM(
    n.sample = n,
    cluster.sample.prop = c(0.51, 0.49),
    delta.methyl = 0,
    delta.expr = 0,
    delta.protein = 0
  )
  
  data_names <- c("methyl", "expr", "protein")
  list_features <- vector("list", length(data_names))
  names(list_features) <- data_names
  list_features[[1]] <- t(sim.data$dat.methyl)
  list_features[[2]] <- t(sim.data$dat.expr)
  list_features[[3]] <- t(sim.data$dat.protein)
  
  feature_table <- as.data.frame(Reduce(rbind, list_features))
  feature_ids <- make.unique(rownames(feature_table))
  rownames(feature_table) <- feature_ids
  colnames(feature_table) <- stringr::str_to_title(colnames(feature_table))
  
  sample_metadata <- sim.data$clustering.assignment
  colnames(sample_metadata) <- c("subjectID", "Y")
  sample_metadata$subjectID <- stringr::str_to_title(sample_metadata$subjectID)
  rownames(sample_metadata) <- sample_metadata$subjectID
  
  row_id <- rep(data_names, vapply(list_features, nrow, integer(1)))
  feature_metadata <- cbind.data.frame(
    featureID = feature_ids,
    featureType = row_id
  )
  rownames(feature_metadata) <- feature_metadata$featureID
  
  list(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  )
}

#' Generate simulated multi-study, multiview data for benchmarking
#'
#' Generates multi-study simulated multiview datasets with specified parameters,
#' including per-study sample size, signal-to-noise ratio, differential
#' expression probabilities, and response generation mode. The output is split
#' into training and testing sets within each study.
#'
#' The model uses a shared and study-specific coefficient structure:
#' for study `s`,
#' `beta_s = sqrt(rho.beta) * beta_shared + sqrt(1 - rho.beta) * beta_indiv_s`,
#' where `beta_shared` is common across studies and `beta_indiv_s` is
#' study-specific.
#'
#' Special cases:
#' - `nstudy = 1`, `rho.beta = 1`, `sigma.alpha = 0`, `tau.snr = 0` recovers
#'   the single-study generator behavior.
#' - `rho.beta = 0` yields independent study-specific coefficients.
#'
#' @param nsample Number of samples per study.
#' @param nstudy Number of studies.
#' @param snr Signal-to-noise ratio for continuous outcomes.
#' @param p.train Train-test split ratio.
#' @param de.prob Differential-expression probability across all modalities.
#' @param de.downProb Down-regulation probability across modalities.
#' @param de.facLoc Differential-expression factor location across modalities.
#' @param de.facScale Differential-expression factor scale across modalities.
#' @param ygen.mode Outcome-generation mode: `"LM"`, `"Friedman"`, or `"Friedman2"`.
#' @param outcome.type Outcome type: `"continuous"`, `"binary"`, or `"survival"`.
#' @param surv.hscale Multiplicative scale on hazard for survival outcomes.
#' @param cens.lower Lower bound for uniform censoring time.
#' @param cens.upper Upper bound for uniform censoring time.
#' @param rho.beta Shared-versus-study-specific mixing level in `[0, 1]`.
#' @param sigma.alpha Standard deviation of study-specific intercepts.
#' @param tau.snr Standard deviation of the log signal-to-noise ratio across studies.
#' @param nrep Number of repetitions.
#' @param seed Random seed.
#'
#' @return A list containing:
#' \itemize{
#'   \item `trainDat`: list of simulated training datasets indexed by repetition, then study.
#'   \item `testDat`: list of simulated testing datasets indexed by repetition, then study.
#'   \item `true_betas`: shared and study-specific coefficient vectors used to generate the outcomes.
#'   \item Simulation settings echoed back, including `snr`, `p.train`, `de.prob`,
#'   `de.downProb`, `de.facLoc`, `de.facScale`, `nrep`, `seed`, `ygen.mode`,
#'   `outcome.type`, `nstudy`, `rho.beta`, `sigma.alpha`, and `tau.snr`.
#' }
#' @export
gen_simmba_multistudy <- function(nsample, # Number of samples per study.
                                  nstudy = 1, # Number of studies.
                                  snr = 1, # Signal-to-noise ratio for continuous outcomes.
                                  p.train = 0.7, # Train-test split ratio.
                                  de.prob = rep(0.1, 3), # Differential-expression probability across all modalities.
                                  de.downProb = rep(0.5, 3), # Down-regulation probability across modalities.
                                  de.facLoc = rep(1, 3), # Differential-expression factor location across modalities.
                                  de.facScale = rep(0.4, 3), # Differential-expression factor scale across modalities.
                                  ygen.mode = c("LM", "Friedman", "Friedman2"), # Outcome-generation mode
                                  outcome.type = c("continuous", "binary", "survival"), # Outcome type
                                  surv.hscale = 1, # Multiplicative scale on hazard for survival outcomes.
                                  cens.lower = 1, # Lower bound for uniform censoring time.
                                  cens.upper = 3, # Upper bound for uniform censoring time.
                                  rho.beta = 1, # Shared-versus-study-specific mixing level in `[0, 1]`.
                                  sigma.alpha = 0, # Standard deviation of study-specific intercepts.
                                  tau.snr = 0, # Standard deviation of the log signal-to-noise ratio across studies.
                                  nrep = 100, # Number of repetitions.
                                  seed = 1234) # Random seed.
{
  set.seed(seed)
  
  ygen.mode <- match.arg(ygen.mode)
  outcome.type <- match.arg(outcome.type)
  
  trainDat <- testDat <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- paste("Rep", seq_len(nrep), sep = "_")
  
  true_betas <- vector("list", nrep)
  names(true_betas) <- paste0("Rep_", seq_len(nrep))
  
  if (nstudy == 1 && rho.beta == 1 && sigma.alpha == 0 && tau.snr == 0) {
    for (k in seq_len(nrep)) {
      pcl <- .trigger_InterSIM_safe(n = nsample)
      X <- as.matrix(t(pcl$feature_table))
      
      nfeature <- table(pcl$feature_metadata$featureType)
      de.facs <- vector("list", 3)
      for (i in seq_len(3)) {
        de.facs[[i]] <- splatter:::getLNormFactors(
          n.facs = nfeature[i],
          sel.prob = de.prob[i],
          neg.prob = de.downProb[i],
          fac.loc = de.facLoc[i],
          fac.scale = de.facScale[i]
        )
      }
      
      beta0 <- log2(unlist(de.facs))
      true_betas[[k]] <- list(
        beta_shared = beta0,
        beta_indiv = list(Study_1 = rep(0, length(beta0))),
        beta_s = list(Study_1 = beta0)
      )
      
      eta_lin <- as.numeric(X %*% beta0)
      
      if (ygen.mode %in% c("Friedman", "Friedman2")) {
        nonzero_index <- which(beta0 != 0)
        if (length(nonzero_index) < 5) {
          stop("Not enough non-zero coefficients to select 5 Friedman features.")
        }
        friedman_index <- sample(nonzero_index, 5)
        X.friedman <- X[, friedman_index, drop = FALSE]
        Xbeta.friedman <- 10 * sin(pi * X.friedman[, 1] * X.friedman[, 2]) +
          20 * (X.friedman[, 3] - 0.5)^2 +
          10 * X.friedman[, 4] +
          5 * X.friedman[, 5]
      }
      
      if (ygen.mode == "LM") {
        eta <- eta_lin
      } else if (ygen.mode == "Friedman") {
        eta <- Xbeta.friedman
      } else {
        eta <- eta_lin + Xbeta.friedman
      }
      
      pcl$sample_metadata$Xbeta <- eta
      
      if (outcome.type == "continuous") {
        sigma2 <- as.vector(stats::var(eta) / snr)
        pcl$sample_metadata$Y <- as.vector(eta + stats::rnorm(nsample) * sqrt(sigma2))
      } else if (outcome.type == "binary") {
        p <- stats::plogis(eta)
        pcl$sample_metadata$Y <- stats::rbinom(nsample, size = 1, prob = p)
      } else {
        h <- as.vector(surv.hscale * exp(eta))
        X0 <- stats::rexp(nsample, rate = h)
        C <- stats::runif(nsample, cens.lower, cens.upper)
        pcl$sample_metadata$time <- ifelse(C >= X0, X0, C)
        pcl$sample_metadata$status <- ifelse(C >= X0, 1L, 0L)
      }
      
      train <- test <- pcl
      tr.row <- sample.int(nsample, size = round(nsample * p.train), replace = FALSE)
      train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
      test$sample_metadata <- pcl$sample_metadata[-tr.row, , drop = FALSE]
      train$feature_table <- pcl$feature_table[, tr.row, drop = FALSE]
      test$feature_table <- pcl$feature_table[, -tr.row, drop = FALSE]
      trainDat[[k]] <- list(Study_1 = train)
      testDat[[k]] <- list(Study_1 = test)
    }
    
    return(list(
      trainDat = trainDat,
      testDat = testDat,
      true_betas = true_betas,
      snr = snr,
      p.train = p.train,
      de.prob = de.prob,
      de.downProb = de.downProb,
      de.facLoc = de.facLoc,
      de.facScale = de.facScale,
      nrep = nrep,
      seed = seed,
      ygen.mode = ygen.mode,
      outcome.type = outcome.type,
      surv.hscale = surv.hscale,
      cens.lower = cens.lower,
      cens.upper = cens.upper,
      nstudy = nstudy,
      rho.beta = rho.beta,
      sigma.alpha = sigma.alpha,
      tau.snr = tau.snr
    ))
  }
  
  for (k in seq_len(nrep)) {
    pcl_template <- .trigger_InterSIM_safe(n = nsample)
    X_template <- as.matrix(t(pcl_template$feature_table))
    nfeature <- table(pcl_template$feature_metadata$featureType)
    p <- nrow(pcl_template$feature_table)
    
    de.facs.shared <- vector("list", 3)
    for (i in seq_len(3)) {
      de.facs.shared[[i]] <- splatter:::getLNormFactors(
        n.facs = nfeature[i],
        sel.prob = de.prob[i],
        neg.prob = de.downProb[i],
        fac.loc = de.facLoc[i],
        fac.scale = de.facScale[i]
      )
    }
    beta_shared <- log2(unlist(de.facs.shared))
    
    train_list <- vector("list", nstudy)
    test_list <- vector("list", nstudy)
    names(train_list) <- names(test_list) <- paste0("Study_", seq_len(nstudy))
    
    true_betas[[k]] <- list(
      beta_shared = beta_shared,
      beta_indiv = vector("list", nstudy),
      beta_s = vector("list", nstudy)
    )
    names(true_betas[[k]]$beta_indiv) <- paste0("Study_", seq_len(nstudy))
    names(true_betas[[k]]$beta_s) <- paste0("Study_", seq_len(nstudy))
    
    for (s in seq_len(nstudy)) {
      if (rho.beta < 1) {
        de.facs.indiv <- vector("list", 3)
        for (i in seq_len(3)) {
          de.facs.indiv[[i]] <- splatter:::getLNormFactors(
            n.facs = nfeature[i],
            sel.prob = de.prob[i],
            neg.prob = de.downProb[i],
            fac.loc = de.facLoc[i],
            fac.scale = de.facScale[i]
          )
        }
        beta_indiv <- log2(unlist(de.facs.indiv))
      } else {
        beta_indiv <- rep(0, length(beta_shared))
      }
      
      beta_s <- sqrt(rho.beta) * beta_shared + sqrt(1 - rho.beta) * beta_indiv
      true_betas[[k]]$beta_indiv[[s]] <- beta_indiv
      true_betas[[k]]$beta_s[[s]] <- beta_s
      
      alpha_s <- if (sigma.alpha > 0) stats::rnorm(1, 0, sigma.alpha) else 0
      log_snr_s <- if (tau.snr > 0) log(snr) + stats::rnorm(1, 0, tau.snr) else log(snr)
      snr_s <- exp(log_snr_s)
      
      pcl <- if (s == 1) pcl_template else .trigger_InterSIM_safe(n = nsample)
      X <- as.matrix(t(pcl$feature_table))
      n_s <- nrow(X)
      eta_lin <- as.numeric(X %*% beta_s) + alpha_s
      
      if (ygen.mode %in% c("Friedman", "Friedman2")) {
        nonzero_index <- which(beta_s != 0)
        if (length(nonzero_index) < 5) {
          stop("Not enough non-zero coefficients to select 5 Friedman features.")
        }
        friedman_index <- sample(nonzero_index, 5)
        X.friedman <- X[, friedman_index, drop = FALSE]
        Xbeta.friedman <- 10 * sin(pi * X.friedman[, 1] * X.friedman[, 2]) +
          20 * (X.friedman[, 3] - 0.5)^2 +
          10 * X.friedman[, 4] +
          5 * X.friedman[, 5]
      }
      
      if (ygen.mode == "LM") {
        eta <- eta_lin
      } else if (ygen.mode == "Friedman") {
        eta <- Xbeta.friedman
      } else {
        eta <- eta_lin + Xbeta.friedman
      }
      
      pcl$sample_metadata$Xbeta <- eta
      pcl$sample_metadata$study <- paste0("Study_", s)
      
      if (outcome.type == "continuous") {
        sigma2 <- as.vector(stats::var(eta) / snr_s)
        pcl$sample_metadata$Y <- as.vector(eta + stats::rnorm(n_s) * sqrt(sigma2))
      } else if (outcome.type == "binary") {
        p_prob <- stats::plogis(eta)
        pcl$sample_metadata$Y <- stats::rbinom(n_s, size = 1, prob = p_prob)
      } else {
        h <- as.vector(surv.hscale * exp(eta))
        X0 <- stats::rexp(n_s, rate = h)
        C <- stats::runif(n_s, cens.lower, cens.upper)
        pcl$sample_metadata$time <- ifelse(C >= X0, X0, C)
        pcl$sample_metadata$status <- ifelse(C >= X0, 1L, 0L)
      }
      
      train <- test <- pcl
      tr.row <- sample.int(n_s, size = round(n_s * p.train), replace = FALSE)
      train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
      test$sample_metadata <- pcl$sample_metadata[-tr.row, , drop = FALSE]
      train$feature_table <- pcl$feature_table[, tr.row, drop = FALSE]
      test$feature_table <- pcl$feature_table[, -tr.row, drop = FALSE]
      train_list[[s]] <- train
      test_list[[s]] <- test
    }
    
    trainDat[[k]] <- train_list
    testDat[[k]] <- test_list
  }
  
  list(
    trainDat = trainDat,
    testDat = testDat,
    true_betas = true_betas,
    snr = snr,
    p.train = p.train,
    de.prob = de.prob,
    de.downProb = de.downProb,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,
    nrep = nrep,
    seed = seed,
    ygen.mode = ygen.mode,
    outcome.type = outcome.type,
    surv.hscale = surv.hscale,
    cens.lower = cens.lower,
    cens.upper = cens.upper,
    nstudy = nstudy,
    rho.beta = rho.beta,
    sigma.alpha = sigma.alpha,
    tau.snr = tau.snr
  )
}

ptmv_build_sim_container <- function(sim_obj, rep_id = "Rep_1", dataset = c("trainDat", "testDat")) {
  dataset <- match.arg(dataset)
  rep_obj <- sim_obj[[dataset]][[rep_id]]
  if (is.null(rep_obj)) {
    stop(sprintf("%s not found in sim_obj$%s.", rep_id, dataset))
  }
  
  study_names <- names(rep_obj)
  reference_feature_metadata <- rep_obj[[study_names[1]]]$feature_metadata
  feature_order <- as.character(reference_feature_metadata$featureID)
  
  feature_table <- do.call(cbind, lapply(study_names, function(study_name) {
    study_obj <- rep_obj[[study_name]]
    if (!identical(as.character(study_obj$feature_metadata$featureID), feature_order)) {
      stop("All studies must share the same feature ordering.")
    }
    ft <- as.matrix(study_obj$feature_table)
    if (!identical(rownames(ft), feature_order)) {
      stop("feature_table rownames must match feature_metadata$featureID.")
    }
    smd <- study_obj$sample_metadata
    sample_ids <- if ("sample_id" %in% colnames(smd)) {
      as.character(smd$sample_id)
    } else if ("subjectID" %in% colnames(smd)) {
      as.character(smd$subjectID)
    } else {
      colnames(ft)
    }
    colnames(ft) <- paste(study_name, sample_ids, sep = "__")
    ft
  }))
  
  sample_metadata <- do.call(rbind, lapply(study_names, function(study_name) {
    study_obj <- rep_obj[[study_name]]
    smd <- as.data.frame(study_obj$sample_metadata, stringsAsFactors = FALSE)
    if (!("sample_id" %in% colnames(smd))) {
      if ("subjectID" %in% colnames(smd)) {
        smd$sample_id <- as.character(smd$subjectID)
      } else {
        smd$sample_id <- colnames(study_obj$feature_table)
      }
    }
    smd$sample_id <- paste(study_name, as.character(smd$sample_id), sep = "__")
    if (!("study" %in% colnames(smd))) {
      smd$study <- study_name
    } else {
      smd$study <- as.character(smd$study)
      smd$study[smd$study == ""] <- study_name
    }
    rownames(smd) <- NULL
    smd
  }))
  
  sample_metadata <- sample_metadata[match(colnames(feature_table), sample_metadata$sample_id), , drop = FALSE]
  if (!identical(sample_metadata$sample_id, colnames(feature_table))) {
    stop("Failed to align sample_metadata to the concatenated feature table.")
  }
  rownames(sample_metadata) <- sample_metadata$sample_id
  
  feature_metadata <- as.data.frame(reference_feature_metadata, stringsAsFactors = FALSE)
  if (!("featureType" %in% colnames(feature_metadata)) && "view" %in% colnames(feature_metadata)) {
    feature_metadata$featureType <- feature_metadata$view
  }
  feature_metadata <- feature_metadata[match(feature_order, feature_metadata$featureID), , drop = FALSE]
  rownames(feature_metadata) <- feature_metadata$featureID
  
  list(
    feature_table = feature_table,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata
  )
}

#' Simulate Gaussian multistudy multiview data for documentation examples
#'
#' Generate one simulated training container and one simulated test container
#' aligned with the `x` input expected by `gptLasso()` and related functions.
#'
#' @inheritParams gen_simmba_multistudy
#'
#' @return A list with:
#' \itemize{
#'   \item `x_train`: a Bioconductor-style multistudy multiview training container.
#'   \item `x_test`: a matching held-out test container with the same feature and view layout.
#' }
#' @examples
#' sim_dat <- sim.gaussian.data()
#' names(sim_dat)
#' str(sim_dat$x_train$sample_metadata)
#' @export
sim.gaussian.data <- function(nsample = 200,
                              nstudy = 3,
                              snr = 5,
                              rho.beta = 0.5,
                              sigma.alpha = 0.1,
                              tau.snr = 0.1,
                              outcome.type = "continuous",
                              nrep = 1,
                              seed = 1234,
                              p.train = 0.7,
                              de.prob = rep(0.1, 3),
                              de.downProb = rep(0.5, 3),
                              de.facLoc = rep(1, 3),
                              de.facScale = rep(0.4, 3),
                              ygen.mode = "LM",
                              surv.hscale = 1,
                              cens.lower = 1,
                              cens.upper = 3) {
  sim_obj <- gen_simmba_multistudy(
    nsample = nsample,
    nstudy = nstudy,
    snr = snr,
    p.train = p.train,
    de.prob = de.prob,
    de.downProb = de.downProb,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,
    ygen.mode = ygen.mode,
    outcome.type = outcome.type,
    surv.hscale = surv.hscale,
    cens.lower = cens.lower,
    cens.upper = cens.upper,
    rho.beta = rho.beta,
    sigma.alpha = sigma.alpha,
    tau.snr = tau.snr,
    nrep = nrep,
    seed = seed
  )
  
  list(
    x_train = ptmv_build_sim_container(sim_obj, rep_id = "Rep_1", dataset = "trainDat"),
    x_test = ptmv_build_sim_container(sim_obj, rep_id = "Rep_1", dataset = "testDat")
  )
}

#' Simulate binary multistudy multiview data for documentation examples
#'
#' Generate one simulated training container and one simulated test container
#' aligned with the `x` input expected by `gptLasso()` and related functions.
#'
#' @inheritParams gen_simmba_multistudy
#'
#' @return A list with:
#' \itemize{
#'   \item `x_train`: a Bioconductor-style multistudy multiview training container.
#'   \item `x_test`: a matching held-out test container with the same feature and view layout.
#' }
#' @examples
#' sim_dat <- sim.binary.data()
#' names(sim_dat)
#' table(sim_dat$x_train$sample_metadata$Y)
#' @export
sim.binary.data <- function(nsample = 300,
                            nstudy = 3,
                            snr = 5,
                            rho.beta = 0.5,
                            sigma.alpha = 0.1,
                            tau.snr = 0.1,
                            outcome.type = "binary",
                            nrep = 1,
                            seed = 5678,
                            p.train = 0.7,
                            de.prob = rep(0.1, 3),
                            de.downProb = rep(0.5, 3),
                            de.facLoc = rep(1, 3),
                            de.facScale = rep(0.4, 3),
                            ygen.mode = "LM",
                            surv.hscale = 1,
                            cens.lower = 1,
                            cens.upper = 3) {
  sim_obj <- gen_simmba_multistudy(
    nsample = nsample,
    nstudy = nstudy,
    snr = snr,
    p.train = p.train,
    de.prob = de.prob,
    de.downProb = de.downProb,
    de.facLoc = de.facLoc,
    de.facScale = de.facScale,
    ygen.mode = ygen.mode,
    outcome.type = outcome.type,
    surv.hscale = surv.hscale,
    cens.lower = cens.lower,
    cens.upper = cens.upper,
    rho.beta = rho.beta,
    sigma.alpha = sigma.alpha,
    tau.snr = tau.snr,
    nrep = nrep,
    seed = seed
  )
  
  list(
    x_train = ptmv_build_sim_container(sim_obj, rep_id = "Rep_1", dataset = "trainDat"),
    x_test = ptmv_build_sim_container(sim_obj, rep_id = "Rep_1", dataset = "testDat")
  )
}


# Helper: choose whether a metric should be minimized or maximized.
# Used in `cv.gptLasso()` when selecting `alpha_ptlasso_hat`.
ptmv_match_metric <- function(type.measure) {
  if (type.measure %in% c("auc")) {
    list(best = which.max, aggregate = max)
  } else {
    list(best = which.min, aggregate = min)
  }
}

# Helper: resolve `s` into a numeric lambda value for a cv-style multiview fit.
# Used in `cv.gptLasso()`, `predict.gptLasso()`, and `predict.cv.gptLasso()`.
ptmv_resolve_s <- function(object, s = c("lambda.min", "lambda.1se")) {
  if (is.numeric(s)) {
    return(as.numeric(s)[1])
  }
  s <- match.arg(s)
  object[[s]]
}

# Helper: stable logit transform used for smoothed binomial baselines.
# Used in `gptLasso()` and prediction helpers.
ptmv_safe_logit <- function(p) log(p / (1 - p))

# Helper: compute one baseline value per study for future prediction.
# Used in `gptLasso()` and `predict.gptLasso()`.
ptmv_compute_group_baseline <- function(y, family) {
  if (family == "gaussian") {
    return(vapply(y, mean, numeric(1)))
  }

  if (family == "binomial") {
    a <- 0.5
    b <- 0.5
    return(vapply(y, function(yi) {
      yi <- as.numeric(yi)
      n1 <- sum(yi == 1)
      n0 <- sum(yi == 0)
      ptmv_safe_logit((n1 + a) / (n1 + n0 + a + b))
    }, numeric(1)))
  }

  stop("group baseline only implemented for gaussian and binomial.")
}

# Helper: compute fold-wise study-specific offsets for stage-one training.
# Used in `gptLasso()`.
ptmv_compute_baseline_offset <- function(y_all, groups_all, foldid_all, family_obj) {
  off <- numeric(length(y_all))
  folds <- sort(unique(foldid_all))
  grps <- levels(groups_all)

  if (family_obj$family == "gaussian") {
    for (f in folds) {
      train <- foldid_all != f
      for (g in grps) {
        idx <- which(groups_all == g & !train)
        if (length(idx) == 0L) {
          next
        }
        off[idx] <- mean(y_all[train & groups_all == g])
      }
    }
    return(off)
  }

  if (family_obj$family == "binomial") {
    a <- 0.5
    b <- 0.5
    for (f in folds) {
      train <- foldid_all != f
      for (g in grps) {
        idx <- which(groups_all == g & !train)
        if (length(idx) == 0L) {
          next
        }
        y_tr_g <- as.numeric(y_all[train & groups_all == g])
        n1 <- sum(y_tr_g == 1)
        n0 <- sum(y_tr_g == 0)
        off[idx] <- ptmv_safe_logit((n1 + a) / (n1 + n0 + a + b))
      }
    }
    return(off)
  }

  stop("Baseline offset only implemented for gaussian and binomial families.")
}

# Helper: remap fold ids to consecutive integers.
# Used in `gptLasso()` and fallback calls.
ptmv_renumber_foldid <- function(fid) {
  used <- sort(unique(as.integer(fid)))
  fold_map <- setNames(seq_along(used), used)
  as.integer(fold_map[as.character(as.integer(fid))])
}

# Helper: create study-level folds with stratification for binomial outcomes.
# Used in `gptLasso()`.
ptmv_make_foldid <- function(n, nfolds, family, y = NULL) {
  nfolds_use <- min(nfolds, n)
  if (family == "binomial") {
    y01 <- as.integer(as.numeric(y) > 0)
    idx0 <- which(y01 == 0L)
    idx1 <- which(y01 == 1L)
    foldid <- integer(n)
    foldid[idx0] <- sample(rep(seq_len(nfolds_use), length.out = length(idx0)))
    foldid[idx1] <- sample(rep(seq_len(nfolds_use), length.out = length(idx1)))
    return(foldid)
  }
  sample(rep(seq_len(nfolds_use), length.out = n))
}

# Helper: switch between serial and study-level parallel lapply.
# Used in `gptLasso()`.
ptmv_maybe_parallel_lapply <- function(X, FUN, parallel = FALSE, ncores = 1L,
                                       verbose = FALSE, ...) {
  if (!isTRUE(parallel) || length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }

  ncores_use <- max(1L, as.integer(ncores))
  if (.Platform$OS.type == "windows" || ncores_use == 1L) {
    if (verbose) {
      message("parallel = TRUE requested, but using serial lapply for this platform/configuration.")
    }
    return(lapply(X, FUN, ...))
  }

  parallel::mclapply(X, FUN, ..., mc.cores = ncores_use)
}

ptmv_require_named_list <- function(x, required_names, context) {
  if (!is.list(x) || length(x) == 0L) {
    stop(sprintf("%s must be a non-empty list.", context))
  }
  missing_names <- setdiff(required_names, names(x))
  if (length(missing_names) > 0L) {
    stop(sprintf(
      "%s must contain the items: %s.",
      context,
      paste(required_names, collapse = ", ")
    ))
  }
}

ptmv_detect_view_column <- function(feature_metadata, context) {
  candidates <- c("featureType", "view", "feature_view", "view_name")
  hit <- intersect(candidates, colnames(feature_metadata))
  if (length(hit) > 0L) {
    return(hit[1])
  }

  if (ncol(feature_metadata) == 2L) {
    return(colnames(feature_metadata)[2])
  }

  stop(sprintf(
    "%s must contain a view column such as 'featureType' or 'view'.",
    context
  ))
}

ptmv_validate_matrix <- function(x, context) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop(sprintf("%s must be a matrix or data.frame coercible to a matrix.", context))
  }
  x
}

# Helper: normalize and validate Bioconductor-style multistudy input.
# Used in `gptLasso()`.
ptmv_normalize_container <- function(x, require_y = TRUE, context = "x") {
  ptmv_require_named_list(x, c("feature_table", "sample_metadata", "feature_metadata"), context)

  feature_table <- ptmv_validate_matrix(x$feature_table, sprintf("%s$feature_table", context))
  sample_metadata <- as.data.frame(x$sample_metadata, stringsAsFactors = FALSE)
  feature_metadata <- as.data.frame(x$feature_metadata, stringsAsFactors = FALSE)

  sample_required <- c("sample_id", "study")
  if (require_y) {
    sample_required <- c(sample_required, "Y")
  }
  sample_missing <- setdiff(sample_required, colnames(sample_metadata))
  if (length(sample_missing) > 0L) {
    stop(sprintf(
      "%s$sample_metadata must contain the columns: %s.",
      context,
      paste(sample_required, collapse = ", ")
    ))
  }

  view_col <- ptmv_detect_view_column(feature_metadata, sprintf("%s$feature_metadata", context))
  feature_required <- c("featureID", view_col)
  feature_missing <- setdiff(feature_required, colnames(feature_metadata))
  if (length(feature_missing) > 0L) {
    stop(sprintf(
      "%s$feature_metadata must contain the columns: %s.",
      context,
      paste(feature_required, collapse = ", ")
    ))
  }

  if (is.null(colnames(feature_table))) {
    stop(sprintf("%s$feature_table must have sample IDs as column names.", context))
  }
  if (is.null(rownames(feature_table))) {
    stop(sprintf("%s$feature_table must have feature IDs as row names.", context))
  }

  sample_id <- as.character(sample_metadata$sample_id)
  feature_id <- as.character(feature_metadata$featureID)
  if (anyNA(sample_id) || any(sample_id == "")) {
    stop(sprintf("%s$sample_metadata$sample_id must be non-missing and non-empty.", context))
  }
  if (anyDuplicated(sample_id)) {
    stop(sprintf("%s$sample_metadata$sample_id must be unique.", context))
  }
  if (anyNA(feature_id) || any(feature_id == "")) {
    stop(sprintf("%s$feature_metadata$featureID must be non-missing and non-empty.", context))
  }
  if (anyDuplicated(feature_id)) {
    stop(sprintf("%s$feature_metadata$featureID must be unique.", context))
  }

  if (!identical(colnames(feature_table), sample_id)) {
    stop(sprintf(
      "Sample order mismatch: %s$sample_metadata$sample_id must exactly match colnames(%s$feature_table).",
      context, context
    ))
  }
  if (!identical(rownames(feature_table), feature_id)) {
    stop(sprintf(
      "Feature order mismatch: %s$feature_metadata$featureID must exactly match rownames(%s$feature_table).",
      context, context
    ))
  }

  study <- as.character(sample_metadata$study)
  if (anyNA(study) || any(study == "")) {
    stop(sprintf("%s$sample_metadata$study must be non-missing and non-empty.", context))
  }
  if (require_y && anyNA(sample_metadata$Y)) {
    stop(sprintf("%s$sample_metadata$Y must be non-missing for training.", context))
  }

  view <- as.character(feature_metadata[[view_col]])
  if (anyNA(view) || any(view == "")) {
    stop(sprintf("%s$feature_metadata[[%s]] must be non-missing and non-empty.", context, view_col))
  }

  study_names <- unique(study)
  view_names <- unique(view)
  if (length(study_names) == 0L) {
    stop(sprintf("%s must contain at least one study.", context))
  }
  if (length(view_names) == 0L) {
    stop(sprintf("%s must contain at least one view.", context))
  }

  features_by_view <- split(feature_id, view)
  features_by_view <- features_by_view[view_names]
  p_by_view <- vapply(features_by_view, length, integer(1))
  if (any(p_by_view == 0L)) {
    empty_views <- names(p_by_view)[p_by_view == 0L]
    stop(sprintf(
      "%s contains view(s) with zero features: %s.",
      context,
      paste(empty_views, collapse = ", ")
    ))
  }

  samples_by_study <- split(sample_id, study)
  samples_by_study <- samples_by_study[study_names]
  n_by_study <- vapply(samples_by_study, length, integer(1))
  if (any(n_by_study == 0L)) {
    empty_studies <- names(n_by_study)[n_by_study == 0L]
    stop(sprintf(
      "%s contains study/studies with zero samples: %s.",
      context,
      paste(empty_studies, collapse = ", ")
    ))
  }

  x_norm <- lapply(study_names, function(study_name) {
    study_samples <- samples_by_study[[study_name]]
    lapply(view_names, function(view_name) {
      feature_ids <- features_by_view[[view_name]]
      mat <- t(feature_table[feature_ids, study_samples, drop = FALSE])
      colnames(mat) <- feature_ids
      rownames(mat) <- study_samples
      mat
    })
  })
  names(x_norm) <- study_names
  for (i in seq_along(x_norm)) {
    names(x_norm[[i]]) <- view_names
  }

  y_norm <- NULL
  y_all <- NULL
  if (require_y) {
    y_norm <- lapply(study_names, function(study_name) {
      idx <- sample_metadata$study == study_name
      drop(sample_metadata$Y[idx])
    })
    names(y_norm) <- study_names
    y_all <- drop(sample_metadata$Y)
  }

  x_list_all <- lapply(view_names, function(view_name) {
    mats <- lapply(x_norm, function(study_obj) study_obj[[view_name]])
    out <- do.call(rbind, mats)
    colnames(out) <- colnames(mats[[1]])
    out
  })
  names(x_list_all) <- view_names

  groups_all <- factor(rep(study_names, times = n_by_study), levels = study_names)

  list(
    x = x_norm,
    y = y_norm,
    k = length(study_names),
    study_names = study_names,
    view_names = view_names,
    n_views = length(view_names),
    n_by_study = n_by_study,
    N_all = sum(n_by_study),
    p_by_view = p_by_view,
    p = sum(p_by_view),
    x_list_all = x_list_all,
    y_all = y_all,
    groups_all = groups_all,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    feature_table = feature_table,
    view_col = view_col,
    sample_ids = sample_id,
    feature_ids = feature_id,
    features_by_view = features_by_view,
    samples_by_study = samples_by_study
  )
}

# Helper: normalize and validate Bioconductor-style test input against training layout.
# Used in `predict.gptLasso()` and `predict.cv.gptLasso()`.
ptmv_normalize_newdata <- function(xtest, template) {
  new_input <- ptmv_normalize_container(xtest, require_y = FALSE, context = "xtest")

  if (!all(new_input$study_names %in% template$study_names)) {
    stop("xtest$sample_metadata$study must be a subset of the training study names.")
  }
  template_subset <- template$study_names[template$study_names %in% new_input$study_names]
  if (!identical(new_input$study_names, template_subset)) {
    stop("xtest$sample_metadata$study must follow the training study ordering for the studies it contains.")
  }
  if (!identical(new_input$view_names, template$view_names)) {
    stop("xtest must use the same views and view ordering as the training data.")
  }
  if (!identical(new_input$feature_ids, template$feature_ids)) {
    stop("xtest$feature_metadata$featureID must exactly match the training feature ordering.")
  }
  if (!identical(new_input$feature_metadata[[new_input$view_col]], template$feature_metadata[[template$view_col]])) {
    stop("xtest$feature_metadata view assignments must exactly match the training data.")
  }

  new_input
}

# Helper: normalize and validate study-list test outcomes.
# Used in `predict.gptLasso()` and `predict.cv.gptLasso()`.
ptmv_prepare_ytest <- function(ytest, study_names, study_sizes = NULL) {
  if (is.null(ytest)) {
    return(NULL)
  }
  if (!is.list(ytest) || length(ytest) != length(study_names)) {
    stop("ytest must be a list with one response vector per study.")
  }
  if (is.null(names(ytest)) || any(names(ytest) == "")) {
    names(ytest) <- study_names
  }
  if (!setequal(names(ytest), study_names)) {
    stop("ytest study names must match the xtest study names exactly.")
  }
  ytest <- ytest[study_names]
  out <- lapply(seq_along(ytest), function(i) {
    yi <- drop(ytest[[i]])
    if (length(yi) == 0L) {
      stop(sprintf("ytest[['%s']] is empty.", study_names[i]))
    }
    if (!is.null(study_sizes) && length(yi) != study_sizes[i]) {
      stop(sprintf(
        "ytest[['%s']] has length %d but xtest has %d samples for that study.",
        study_names[i], length(yi), study_sizes[i]
      ))
    }
    yi
  })
  names(out) <- study_names
  out
}

# Helper: expand one study-level baseline per study into one offset per sample.
# Used in `predict.gptLasso()`.
ptmv_study_offsets <- function(study_sizes, group_baseline) {
  unlist(Map(function(n, b) rep(b, n), study_sizes, as.list(group_baseline)), use.names = FALSE)
}

# Helper: stack study data into one multiview object view by view.
# Used in `predict.gptLasso()`.
ptmv_stack_by_view <- function(study_list, view_names) {
  out <- lapply(view_names, function(v) {
    do.call(rbind, lapply(study_list, function(study) study[[v]]))
  })
  names(out) <- view_names
  out
}

# Helper: split a stacked prediction vector back into studies.
# Used in `gptLasso()` and `predict.gptLasso()`.
ptmv_split_vector_by_study <- function(x, study_sizes, study_names) {
  split(x, rep(study_names, times = study_sizes))[study_names]
}

# Helper: predict one list of study-specific models and return study-wise outputs.
# Used in `predict.gptLasso()`.
ptmv_predict_by_study <- function(model_list, xtest, s, type, offsets = NULL) {
  study_names <- names(xtest)
  out <- vector("list", length(study_names))
  names(out) <- study_names
  for (study_name in study_names) {
    model <- model_list[[study_name]]
    newoffset <- if (is.null(offsets)) NULL else offsets[[study_name]]
    pred <- predict(model, newx = xtest[[study_name]], s = s, type = type, newoffset = newoffset)
    out[[study_name]] <- as.numeric(pred)
  }
  out
}

# Helper: extract nonzero support from one cvar/cv.multiview object at `s`.
# Used in `predict.gptLasso()`.
ptmv_get_support <- function(model, s) {
  lambda <- ptmv_resolve_s(model, s)
  beta <- as.numeric(coef(model$multiview.fit, s = lambda))
  which(beta[-1] != 0)
}

# Helper: union nonzero support across multiple multiview fits.
# Used in `predict.gptLasso()`.
ptmv_get_union_support <- function(models, s) {
  sort(unique(unlist(lapply(models, ptmv_get_support, s = s))))
}

# Helper: compute one scalar performance metric for one study or stacked data.
# Used in `cv.gptLasso()`, `predict.gptLasso()`, and `predict.cv.gptLasso()`.
ptmv_metric_value <- function(y, pred, family, type.measure) {
  y <- as.numeric(y)
  pred <- as.numeric(pred)

  if (family == "gaussian") {
    if (type.measure %in% c("mse", "deviance")) {
      return(mean((y - pred)^2))
    }
    stop("Unsupported type.measure for gaussian.")
  }

  if (family == "binomial") {
    if (type.measure == "auc") {
      if (length(unique(stats::na.omit(y))) < 2L) {
        return(NA_real_)
      }
      return(as.numeric(pROC::auc(
        pROC::roc(response = y, predictor = pred, levels = c(0, 1), direction = "<", quiet = TRUE)
      )))
    }
    if (type.measure == "deviance") {
      eps <- 1e-8
      pp <- pmin(pmax(pred, eps), 1 - eps)
      return(-2 * mean(y * log(pp) + (1 - y) * log(1 - pp)))
    }
    stop("Unsupported type.measure for binomial.")
  }

  stop("Unsupported family.")
}

# Helper: summarize study-wise predictions into overall and per-study metrics.
# Used in `cv.gptLasso()`, `predict.gptLasso()`, and `predict.cv.gptLasso()`.
ptmv_summarize_metric <- function(preds, y, family, type.measure, add_r2 = FALSE) {
  study_names <- names(y)
  study_err <- vapply(study_names, function(study_name) {
    ptmv_metric_value(y[[study_name]], preds[[study_name]], family, type.measure)
  }, numeric(1))

  all_y <- unlist(y, use.names = FALSE)
  all_pred <- unlist(preds, use.names = FALSE)
  out <- c(
    allGroups = ptmv_metric_value(all_y, all_pred, family, type.measure),
    mean = mean(study_err, na.rm = TRUE),
    setNames(study_err, paste0("group_", study_names))
  )

  if (add_r2 && family == "gaussian") {
    out <- c(out, "r^2" = 1 - sum((all_y - all_pred)^2) / sum((all_y - mean(all_y))^2))
  }
  out
}

# Helper: assemble a common prediction object for direct-fit and cv-fit methods.
# Used in `predict.gptLasso()` and `predict.cv.gptLasso()`.
ptmv_build_prediction_object <- function(call, alpha_ptlasso, type.measure,
                                         yhatoverall, yhatind, yhatpre,
                                         supoverall, supind, suppre.common, suppre.individual,
                                         linkoverall = NULL, linkind = NULL, linkpre = NULL,
                                         erroverall = NULL, errind = NULL, errpre = NULL,
                                         fit = NULL, class_name = "predict.gptLasso") {
  out <- list(
    call = call,
    alpha_ptlasso = alpha_ptlasso,
    yhatoverall = yhatoverall,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = suppre.individual,
    type.measure = type.measure
  )
  if (!is.null(linkoverall)) out$linkoverall <- linkoverall
  if (!is.null(linkind)) out$linkind <- linkind
  if (!is.null(linkpre)) out$linkpre <- linkpre
  if (!is.null(erroverall)) out$erroverall <- erroverall
  if (!is.null(errind)) out$errind <- errind
  if (!is.null(errpre)) out$errpre <- errpre
  if (!is.null(fit)) out$fit <- fit
  if (identical(class_name, "predict.cv.gptLasso")) {
    class(out) <- c("predict.cv.gptLasso", "predict.gptLasso")
  } else {
    class(out) <- "predict.gptLasso"
  }
  out
}
