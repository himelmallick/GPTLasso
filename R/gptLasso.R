#' Fit the multistudy multiview gptLasso pipeline
#'
#' Fit the overall, study-specific individual, and pretrained transfer models
#' for a Bioconductor-style multistudy multiview input container.
#'
#' @param x A named list with entries `feature_table`, `sample_metadata`, and
#'   `feature_metadata`. The feature table must be feature-by-sample, the sample
#'   metadata must contain `sample_id`, `study`, and `Y`, and the feature
#'   metadata must contain `featureID` plus a view-mapping column such as
#'   `featureType`.
#' @param alpha.ptlasso Transfer-learning level in `[0, 1]`.
#' @param family Response family. Currently `"gaussian"` and `"binomial"` are supported.
#' @param type.measure Cross-validation metric optimized inside the multiview fits.
#' @param rho Multiview cooperative-learning fusion parameter, supplied as one value or a tuning grid.
#' @param overall.lambda Lambda rule used for the stage-one overall model.
#' @param ind.lambda Lambda rule used for the individual models.
#' @param pre.lambda Lambda rule used for the pretrained models.
#' @param target Optional target study/studies. Default `NULL` reproduces current
#'   behavior and fits study-specific models for all studies. When provided,
#'   must be a character vector of study names; `fitind`/`fitpre` are restricted
#'   to these targets (in training-study order), while `fitoverall` is still fit
#'   on all studies.
#' @param foldid Optional stacked fold assignment across all studies.
#' @param nfolds Number of folds used when `foldid` is not supplied.
#' @param verbose Should progress messages be printed?
#' @param fitoverall Optional pre-fit overall multiview model to reuse.
#' @param fitind Optional pre-fit list of study-specific multiview models to reuse.
#' @param penalty.factor Optional penalty factors across concatenated views.
#' @param group.intercepts Should study-specific stage-one baselines be used.
#' @param alpha.glmnet Elastic-net mixing parameter passed to the multiview base learner, default is 1 indicating Lasso regression.
#' @param parallel Logical; if `TRUE`, allow study-level parallel fits where available.
#' @param ncores Number of worker cores for study-level parallel fits.
#' @param ... Additional arguments forwarded to the multiview base fitter.
#'
#' @return A `gptLasso` object, stored as a list. Important components include:
#' \itemize{
#'   \item `study.names`: study labels detected from `x$sample_metadata$study`.
#'   \item `view.names`: view labels detected from `x$feature_metadata`.
#'   \item `n.by.study`: number of samples in each study.
#'   \item `alpha.ptlasso`: the transfer-learning level used for the fit.
#'   \item `fitoverall`: the pooled multiview stage-one fit.
#'   \item `fitind`: a named list of study-specific individual fits.
#'   \item `fitpre`: a named list of study-specific pretrained transfer fits.
#'   \item `target`: the user-supplied target specification.
#'   \item `target.study.names`: normalized target studies used for study-specific fitting.
#'   \item `support.vars`: feature indices selected by the overall stage-one fit.
#'   \item `group.baseline`: study-level baseline offsets used when `group.intercepts = TRUE`.
#'   \item `preval.offset`: study-wise stage-one linear predictors reused as offsets in pretraining.
#'   \item `training.layout`: the normalized training container split by study and view.
#' }
#'
#' @examples
#' # Gaussian example
#' set.seed(1234)
#' sim_dat <- sim.gaussian.data()
#' x_train <- sim_dat$x_train
#' str(x_train)
#'
#' fit_gaussian <- gptLasso(
#'   x = x_train,
#'   alpha.ptlasso = 0.5,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   nfolds = 3
#' )
#'
#' fit_gaussian_target <- gptLasso(
#'   x = x_train,
#'   alpha.ptlasso = 0.5,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   target = "Study_2",
#'   nfolds = 3
#' )
#'
#' names(fit_gaussian)
#' fit_gaussian$study.names
#' fit_gaussian$view.names
#' fit_gaussian$n.by.study
#' fit_gaussian$support.vars
#'
#' # Binomial example
#' set.seed(5678)
#' sim_dat_bin <- sim.binary.data()
#' x_train_bin <- sim_dat_bin$x_train
#'
#' fit_binomial <- gptLasso(
#'   x = x_train_bin,
#'   alpha.ptlasso = 0.5,
#'   family = "binomial",
#'   type.measure = "auc",
#'   nfolds = 3
#' )
#'
#' names(fit_binomial)
#' fit_binomial$group.baseline
#' @export
gptLasso <- function(
    x,
    alpha.ptlasso = 0.5,
    family = c("gaussian", "binomial"),
    type.measure = c("default", "mse", "auc", "deviance"),
    rho = seq(0, 1, length = 11),
    overall.lambda = c("lambda.1se", "lambda.min"),
    ind.lambda = c("lambda.1se", "lambda.min"),
    pre.lambda = c("lambda.1se", "lambda.min"),
    target = NULL,
    foldid = NULL,
    nfolds = 10,
    verbose = FALSE,
    fitoverall = NULL,
    fitind = NULL,
    penalty.factor = NULL,
    group.intercepts = TRUE,
    alpha.glmnet = 1,
    parallel = FALSE,
    ncores = 1L,
    ...
) {
  this.call <- match.call()
  dot.args <- list(...)
  legacy.names <- c("alpha_ptlasso", "alpha_glmnet")
  bad.args <- intersect(names(dot.args), legacy.names)
  if (length(bad.args) > 0L) {
    stop(
      sprintf(
        "Legacy underscore argument(s) not supported in gptLasso(): %s. Use dot-style names instead.",
        paste(bad.args, collapse = ", ")
      )
    )
  }
  
  family <- match.arg(family)
  type.measure <- match.arg(type.measure)
  if (type.measure == "default") {
    type.measure <- if (family == "gaussian") "mse" else "deviance"
  }
  
  overall.lambda <- match.arg(overall.lambda, c("lambda.1se", "lambda.min"))
  ind.lambda <- match.arg(ind.lambda, c("lambda.1se", "lambda.min"))
  pre.lambda <- match.arg(pre.lambda, c("lambda.1se", "lambda.min"))
  
  if (!is.numeric(alpha.ptlasso) || length(alpha.ptlasso) != 1L ||
      alpha.ptlasso < 0 || alpha.ptlasso > 1) {
    stop("alpha.ptlasso must be a single number between 0 and 1.")
  }
  if (!is.numeric(alpha.glmnet) || length(alpha.glmnet) != 1L ||
      alpha.glmnet < 0 || alpha.glmnet > 1) {
    stop("alpha.glmnet must be a single number between 0 and 1.")
  }
  rho <- unique(as.numeric(rho))
  if (length(rho) < 1L || any(is.na(rho))) {
    stop("rho must contain at least one numeric value.")
  }
  if (!is.numeric(nfolds) || length(nfolds) != 1L || nfolds < 2) {
    stop("nfolds must be a single integer greater than or equal to 2.")
  }
  
  family_fn <- switch(family, gaussian = gaussian, binomial = binomial)
  input <- ptmv_normalize_container(x, require_y = TRUE, context = "x")
  x <- input$x
  y <- input$y
  k <- input$k
  study_names <- input$study_names
  view_names <- input$view_names
  n_views <- input$n_views
  n_by_study <- input$n_by_study
  N_all <- input$N_all
  p_by_view <- input$p_by_view
  p <- input$p
  x_list_all <- input$x_list_all
  y_all <- input$y_all
  groups_all <- input$groups_all
  target.study.names <- ptmv_resolve_target_studies(
    target = target,
    study_names = study_names,
    context = "gptLasso(target)"
  )
  k_target <- length(target.study.names)
  
  if (k == 1L) {
    message("Single-study input detected; falling back to cvar.multiview().")
    return(cvar.multiview(
      x.list = x[[1]],
      y = y[[1]],
      family = family_fn(),
      alpha = alpha.glmnet,
      rho = rho,
      s = overall.lambda,
      nfolds = min(nfolds, length(y[[1]])),
      foldid = if (is.null(foldid)) NULL else ptmv_renumber_foldid(foldid),
      penalty.factor = penalty.factor,
      type.measure = type.measure,
      keep = TRUE,
      ...
    ))
  }
  
  if (n_views == 1L) {
    message("Single-view multi-study input detected; falling back to ptLasso().")
    x_single <- x_list_all[[1]]
    if (is.null(penalty.factor)) {
      penalty.factor <- rep(1, ncol(x_single))
    }
    return(ptLasso::ptLasso(
      x = x_single,
      y = y_all,
      groups = groups_all,
      alpha = alpha.ptlasso,
      family = family,
      type.measure = type.measure,
      use.case = "inputGroups",
      overall.lambda = overall.lambda,
      foldid = foldid,
      nfolds = nfolds,
      verbose = verbose,
      penalty.factor = penalty.factor,
      fitoverall = fitoverall,
      fitind = fitind,
      en.alpha = alpha.glmnet,
      group.intercepts = group.intercepts,
      parallel = parallel,
      ...
    ))
  }
  
  if (is.null(penalty.factor)) {
    penalty.factor <- rep(1, p)
  }
  if (length(penalty.factor) != p) {
    stop(sprintf("penalty.factor must have length %d, the total number of features across views.", p))
  }
  
  if (!is.null(fitoverall)) {
    valid_overall <- inherits(fitoverall, "cv.multiview") ||
      inherits(fitoverall, "cvar.multiview") ||
      inherits(fitoverall, "cv.multiview.revised")
    if (!valid_overall) {
      stop("fitoverall must be a cv.multiview/cvar.multiview/cv.multiview.revised object.")
    }
  }
  
  if (!is.null(fitind)) {
    if (length(fitind) != k_target) {
      stop("fitind must contain one model per target study.")
    }
    valid_ind <- vapply(
      fitind,
      function(z) inherits(z, "cv.multiview") ||
        inherits(z, "cvar.multiview"),
      logical(1)
    )
    if (!all(valid_ind)) {
      stop("All elements of fitind must be cv.multiview/cvar.multiview objects.")
    }
    if (is.null(names(fitind)) || any(names(fitind) == "")) {
      names(fitind) <- target.study.names
    }
    if (!setequal(names(fitind), target.study.names)) {
      stop("fitind names must match normalized target study names.")
    }
    fitind <- fitind[target.study.names]
  }
  
  if (is.null(foldid)) {
    foldid_all <- integer(N_all)
    start <- 1L
    for (g in seq_len(k)) {
      idx <- start:(start + n_by_study[g] - 1L)
      foldid_all[idx] <- ptmv_make_foldid(
        n = n_by_study[g],
        nfolds = nfolds,
        family = family,
        y = y[[g]]
      )
      start <- start + n_by_study[g]
    }
    foldid_all <- ptmv_renumber_foldid(foldid_all)
  } else {
    if (length(foldid) != N_all) {
      stop("Provided foldid must have length equal to total N across studies.")
    }
    foldid_all <- ptmv_renumber_foldid(foldid)
  }
  nfolds_all <- length(unique(foldid_all))
  
  foldid_within <- split(foldid_all, groups_all)
  foldid_within <- lapply(foldid_within, ptmv_renumber_foldid)
  foldid_within <- foldid_within[target.study.names]
  
  group_baseline <- ptmv_compute_group_baseline(y, family)
  baseline_offset_all <- NULL
  if (isTRUE(group.intercepts)) {
    baseline_offset_all <- ptmv_compute_baseline_offset(
      y_all = y_all,
      groups_all = groups_all,
      foldid_all = foldid_all,
      family_obj = family_fn()
    )
  }
  
  if (is.null(fitoverall)) {
    if (verbose) {
      message("Fitting overall multiview model.")
    }
    fitoverall <- cvar.multiview(
      x.list = x_list_all,
      y = y_all,
      family = family_fn(),
      alpha = alpha.glmnet,
      rho = rho,
      s = overall.lambda,
      type.measure = type.measure,
      foldid = foldid_all,
      nfolds = nfolds_all,
      offset = baseline_offset_all,
      penalty.factor = penalty.factor,
      keep = TRUE,
      ...
    )
  }
  
  fitoverall_fit <- fitoverall$multiview.fit
  lamhat <- fitoverall$lambda.choice
  rho_overall <- fitoverall$rho.choice
  r_idx <- which(fitoverall$rho == rho_overall)[1]
  rho_obj <- fitoverall$cv_by_rho[[r_idx]]
  lam_idx <- which.min(abs(rho_obj$lambda - lamhat))
  preval_all <- as.numeric(rho_obj$fit.preval[, lam_idx])
  preval.offset <- ptmv_split_vector_by_study(preval_all, n_by_study, study_names)
  preval.offset.target <- preval.offset[target.study.names]
  
  coef_vec <- as.numeric(coef(fitoverall_fit, s = lamhat))
  supall <- which(coef_vec[-1] != 0)
  x_target <- x[target.study.names]
  y_target <- y[target.study.names]
  
  if (is.null(fitind)) {
    if (verbose) {
      message("Fitting individual multiview models.")
    }
    fitind <- ptmv_maybe_parallel_lapply(seq_len(k_target), function(kk) {
      if (verbose) {
        message(sprintf("  Study %d / %d", kk, k_target))
      }
      cvar.multiview(
        x.list = x_target[[kk]],
        y = y_target[[kk]],
        family = family_fn(),
        alpha = alpha.glmnet,
        rho = rho,
        foldid = foldid_within[[kk]],
        s = ind.lambda,
        type.measure = type.measure,
        keep = TRUE,
        ...
      )
    }, parallel = parallel, ncores = ncores, verbose = verbose)
    names(fitind) <- target.study.names
  }
  
  if (verbose) {
    message("Fitting pretrained multiview models.")
  }
  
  if (alpha.ptlasso == 1) {
    fitpre <- fitind
  } else {
    fitpre <- ptmv_maybe_parallel_lapply(seq_len(k_target), function(kk) {
      if (verbose) {
        message(sprintf("  Pretrained model %d / %d", kk, k_target))
      }
      alpha_eff <- max(alpha.ptlasso, 1e-9)
      fac <- rep(1 / alpha_eff, p)
      fac[supall] <- 1
      pf <- penalty.factor * fac
      if (alpha.ptlasso == 0 && length(supall) == 0L) {
        pf <- penalty.factor * rep(1e9, p)
      }
      cvar.multiview(
        x.list = x_target[[kk]],
        y = y_target[[kk]],
        family = family_fn(),
        alpha = alpha.glmnet,
        rho = rho,
        foldid = foldid_within[[kk]],
        s = pre.lambda,
        type.measure = type.measure,
        offset = (1 - alpha.ptlasso) * preval.offset.target[[kk]],
        penalty.factor = pf,
        keep = TRUE,
        ...
      )
    }, parallel = parallel, ncores = ncores, verbose = verbose)
    names(fitpre) <- target.study.names
  }
  
  out <- list(
    call = this.call,
    k = k,
    N.all = N_all,
    n.by.study = n_by_study,
    study.names = study_names,
    group.levels = study_names,
    view.names = view_names,
    n.views = n_views,
    p.by.view = p_by_view,
    features.all = p,
    alpha.ptlasso = alpha.ptlasso,
    alpha.glmnet = alpha.glmnet,
    target = target,
    target.study.names = target.study.names,
    rho = rho,
    family = family,
    type.measure = type.measure,
    overall.lambda = overall.lambda,
    ind.lambda = ind.lambda,
    pre.lambda = pre.lambda,
    fitoverall.lambda = lamhat,
    fitoverall.rho = rho_overall,
    fitind.rho = stats::setNames(vapply(fitind, function(model) {
      if (is.null(model$rho.choice)) NA_real_ else model$rho.choice
    }, numeric(1)), target.study.names),
    fitpre.rho = stats::setNames(vapply(fitpre, function(model) {
      if (is.null(model$rho.choice)) NA_real_ else model$rho.choice
    }, numeric(1)), target.study.names),
    group.intercepts = group.intercepts,
    group.baseline = group_baseline,
    foldid = foldid_all,
    foldid.within = foldid_within,
    support.vars = supall,
    penalty.factor = penalty.factor,
    fitoverall = fitoverall,
    fitind = fitind,
    fitpre = fitpre,
    baseline.offset = baseline_offset_all,
    preval.offset = preval.offset,
    training.layout = input,
    parallel = parallel,
    ncores = max(1L, as.integer(ncores))
  )
  
  class(out) <- "gptLasso"
  out
}
