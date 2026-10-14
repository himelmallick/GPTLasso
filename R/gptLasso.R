#' Fit the multistudy multiview gptLasso pipeline
#'
#' Fit the overall, study-specific individual, and pretrained transfer models
#' for a Bioconductor-style multistudy multiview input container.
#'
#' @param x Either a [MultiAssayExperiment::MultiAssayExperiment] with one
#'   named [SummarizedExperiment::SummarizedExperiment] per view, or a PCL-style
#'   list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#'   PCL sample metadata must contain `sample_id`, `study`, and `Y`; feature
#'   metadata must contain `featureID` and a recognized view column.
#' @param assay.names Assays to use for `MultiAssayExperiment` input. Must be
#'   `NULL` for PCL input. For MAE input, `NULL` requires one assay per experiment;
#'   a scalar selects the same assay name in every experiment; a named character
#'   vector selects one assay per experiment.
#' @param alpha.ptlasso Transfer-learning level in `[0, 1]`.
#' @param family Response family: `"gaussian"`, `"binomial"`, or `"cox"`.
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
#' @return A `gptLasso` object, stored as a list with these components:
#' \itemize{
#'   \item `call`: the matched fitting call.
#'   \item `k`, `N.all`, `n.by.study`: numbers of studies, total samples, and
#'   samples per study.
#'   \item `study.names`, `group.levels`: study labels detected from
#'   the input sample metadata, in fitting order.
#'   \item `view.names`, `n.views`, `p.by.view`, `features.all`: input view
#'   names and feature-count summaries.
#'   \item `alpha.ptlasso`, `alpha.ptlasso.requested`,
#'   `transfer.applicable`: effective/requested transfer values and whether
#'   cross-study transfer applies.
#'   \item `rho`, `rho.requested`, `fusion.applicable`: effective/requested
#'   multiview fusion grids and whether fusion tuning applies.
#'   \item `assay.names`: the selected assay name for each experiment.
#'   \item `alpha.glmnet`, `family`, `type.measure`: base-learner and response
#'   configuration.
#'   \item `target`, `target.study.names`: supplied and normalized targets.
#'   \item `overall.lambda`, `ind.lambda`, `pre.lambda`, `fitoverall.lambda`:
#'   requested lambda rules and the selected overall-model lambda.
#'   \item `fitoverall.rho`, `fitind.rho`, `fitpre.rho`: selected rho values
#'   for the overall, individual, and pretrained fits.
#'   \item `group.intercepts`, `group.intercepts.requested`, `group.baseline`,
#'   `cox.study.effect`: effective/requested study-baseline configuration.
#'   \item `foldid`, `foldid.within`: stacked and within-study fold assignments.
#'   \item `support.vars`, `penalty.factor`: selected overall-model feature
#'   indices and applied penalty factors.
#'   \item `fitoverall`: the overall all-study multiview fit.
#'   \item `fitind`: a named list of study-specific individual fits.
#'   \item `fitpre`: a named list of study-specific pretrained transfer fits.
#'   \item `baseline.offset`, `preval.offset`: study-baseline and stage-one
#'   linear-predictor offsets.
#'   \item `training.layout`: the internal study-by-view matrix layout derived
#'   from the input container.
#'   \item `parallel`, `ncores`: effective study-level parallel settings.
#' }
#'
#' @examples
#' # Gaussian example
#' set.seed(1234)
#' sim_dat <- sim.gaussian.data(nsample = c(15, 15), nstudy = 2,
#'                              de.prob = rep(0.5, 3))
#' x_train <- sim_dat$x_train
#' names(x_train)
#' colnames(x_train$sample_metadata)
#' table(x_train$sample_metadata$study)
#'
#' fit_gaussian <- gptLasso(
#'   x = x_train,
#'   alpha.ptlasso = 0.5,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   rho = 0,
#'   nfolds = 3
#' )
#'
#' fit_gaussian_target <- gptLasso(
#'   x = x_train,
#'   alpha.ptlasso = 0.5,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   rho = 0,
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
#' @export
gptLasso <- function(
    x,
    assay.names = NULL,
    alpha.ptlasso = 0.5,
    family = c("gaussian", "binomial", "cox"),
    type.measure = c("default", "mse", "auc", "deviance", "C"),
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
  allowed.metrics <- switch(
    family,
    gaussian = c("mse", "deviance"),
    binomial = c("auc", "deviance"),
    cox = c("deviance", "C")
  )
  if (!type.measure %in% allowed.metrics) {
    stop(sprintf(
      "For %s family, type.measure must be one of: %s.",
      family, paste(allowed.metrics, collapse = ", ")
    ))
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
  min.folds <- if (family == "cox") 3L else 2L
  if (!is.numeric(nfolds) || length(nfolds) != 1L || nfolds < min.folds) {
    stop(sprintf("nfolds must be a single integer greater than or equal to %d.", min.folds))
  }
  
  family_arg <- switch(family, gaussian = stats::gaussian(),
                       binomial = stats::binomial(), cox = "cox")
  input <- ptmv_normalize_container(
    x, require_y = TRUE, context = "x", family = family,
    assay.names = assay.names
  )
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
  alpha.ptlasso.requested <- alpha.ptlasso
  rho.requested <- rho
  group.intercepts.requested <- group.intercepts
  transfer.applicable <- k > 1L
  fusion.applicable <- n_views > 1L
  if (!transfer.applicable) {
    alpha.ptlasso <- 1
    group.intercepts <- FALSE
    fitind <- NULL
  }
  if (!fusion.applicable) rho <- 0
  target.study.names <- ptmv_resolve_target_studies(
    target = target,
    study_names = study_names,
    context = "gptLasso(target)"
  )
  k_target <- length(target.study.names)
  
  if (family == "cox") {
    event.counts <- vapply(y, function(yi) sum(as.numeric(yi[, 2]) == 1), integer(1))
    if (any(event.counts < 3L)) {
      bad <- names(event.counts)[event.counts < 3L]
      stop(sprintf(
        "Cox fitting requires at least 3 observed events in every study; insufficient events in: %s.",
        paste(bad, collapse = ", ")
      ))
    }
    if (is.null(foldid)) {
      nfolds <- min(nfolds, min(event.counts), min(n_by_study))
    } else {
      if (length(foldid) != N_all) {
        stop("Provided foldid must have length equal to total N across studies.")
      }
      ptmv_validate_cox_folds(
        y_all, groups_all, ptmv_renumber_foldid(foldid), context = "foldid"
      )
    }
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
  if (family == "cox") {
    ptmv_validate_cox_folds(y_all, groups_all, foldid_all, context = "foldid")
  }
  
  foldid_within <- split(foldid_all, groups_all)
  foldid_within <- lapply(foldid_within, ptmv_renumber_foldid)
  foldid_within <- foldid_within[target.study.names]
  
  group_baseline <- if (family == "cox") {
    stats::setNames(rep(0, k), study_names)
  } else {
    ptmv_compute_group_baseline(y, family)
  }
  baseline_offset_all <- NULL
  if (isTRUE(group.intercepts) && family != "cox") {
    baseline_offset_all <- ptmv_compute_baseline_offset(
      y_all = y_all,
      groups_all = groups_all,
      foldid_all = foldid_all,
      family_obj = family_arg
    )
  }
  
  if (is.null(fitoverall)) {
    if (verbose) {
      message("Fitting overall multiview model.")
    }
    fitoverall <- cvar.multiview(
      x.list = x_list_all,
      y = y_all,
      family = family_arg,
      alpha = alpha.glmnet,
      rho = rho,
      s = overall.lambda,
      type.measure = type.measure,
      foldid = foldid_all,
      nfolds = nfolds_all,
      offset = baseline_offset_all,
      cox.groups = if (family == "cox" && isTRUE(group.intercepts)) groups_all else NULL,
      penalty.factor = penalty.factor,
      keep = TRUE,
      ...
    )
  }
  
  if (!identical(ptmv_family_name(fitoverall$family), family)) {
    stop("fitoverall family is not compatible with the requested family.")
  }
  
  if (family == "cox" && isTRUE(group.intercepts)) {
    if (!identical(fitoverall$cox.study.effect, "fixed_shift") ||
        is.null(fitoverall$cox.group.shift)) {
      stop("A supplied Cox fitoverall with group.intercepts = TRUE must contain compatible fixed study shifts.")
    }
    if (!setequal(names(fitoverall$cox.group.shift), study_names)) {
      stop("Stored Cox study shifts must be named with the training study names.")
    }
    group_baseline <- fitoverall$cox.group.shift[study_names]
    baseline_offset_all <- fitoverall$offset.full
  } else if (family == "cox" && !is.null(fitoverall$cox.study.effect)) {
    stop("A supplied Cox fitoverall with group.intercepts = FALSE must not contain study shifts.")
  }
  
  lamhat <- fitoverall$lambda.choice
  rho_overall <- fitoverall$rho.choice
  r_idx <- which(fitoverall$rho == rho_overall)[1]
  rho_obj <- fitoverall$cv_by_rho[[r_idx]]
  lam_idx <- which.min(abs(rho_obj$lambda - lamhat))
  preval_all <- as.numeric(rho_obj$fit.preval[, lam_idx])
  preval.offset <- ptmv_split_vector_by_study(preval_all, n_by_study, study_names)
  preval.offset.target <- preval.offset[target.study.names]
  
  supall <- ptmv_get_support(fitoverall, lamhat)
  x_target <- x[target.study.names]
  y_target <- y[target.study.names]
  
  if (!transfer.applicable) {
    fitind <- stats::setNames(list(fitoverall), study_names)
  } else if (is.null(fitind)) {
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
        family = family_arg,
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
  
  if (!transfer.applicable || alpha.ptlasso == 1) {
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
        family = family_arg,
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
    alpha.ptlasso.requested = alpha.ptlasso.requested,
    transfer.applicable = transfer.applicable,
    alpha.glmnet = alpha.glmnet,
    target = target,
    target.study.names = target.study.names,
    rho = rho,
    rho.requested = rho.requested,
    fusion.applicable = fusion.applicable,
    assay.names = input$assay_names,
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
    group.intercepts.requested = group.intercepts.requested,
    group.baseline = group_baseline,
    cox.study.effect = if (family == "cox" && isTRUE(group.intercepts)) "fixed_shift" else NULL,
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