#' Cross-validate the multistudy multiview gptLasso pipeline
#'
#' Refit `gptLasso()` over a grid of transfer-learning values and summarize the
#' pretrained performance profile used to select the final `alpha.ptlasso`.
#'
#' @param x A named list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#' @param alpha.ptlasso.list Numeric vector of transfer-learning values to compare.
#' @param family Response family. Currently `"gaussian"` and `"binomial"` are supported.
#' @param type.measure Cross-validation metric used to compare transfer-learning levels.
#' @param rho Multiview cooperative-learning fusion parameter, supplied as one value or a tuning grid.
#' @param nfolds Number of folds used inside each `gptLasso()` fit.
#' @param foldid Optional stacked fold assignment across all studies.
#' @param overall.lambda Lambda rule used for the overall model when summarizing CV performance.
#' @param ind.lambda Lambda rule used for the individual models when summarizing CV performance.
#' @param pre.lambda Lambda rule used for the pretrained models when summarizing CV performance.
#' @param verbose Should progress messages be printed?
#' @param fitoverall Optional pre-fit overall multiview model reused across alphas.
#' @param fitind Optional pre-fit list of study-specific multiview models reused across alphas.
#' @param group.intercepts Should study-specific stage-one baselines be used?
#' @param parallel Logical; if `TRUE`, allow study-level parallel fits where available.
#' @param ncores Number of worker cores for study-level parallel fits.
#' @param ... Additional arguments forwarded to `gptLasso()`.
#'
#' @return A `cv.gptLasso` object, stored as a list. Important components include:
#' \itemize{
#'   \item `alpha.ptlasso.hat`: the selected fixed transfer-learning value.
#'   \item `varying.alpha.ptlasso.hat`: study-specific transfer-learning values chosen from the same grid.
#'   \item `alpha.ptlasso.list`: the candidate transfer-learning grid that was evaluated.
#'   \item `errpre`: a matrix summarizing pretrained performance for each candidate alpha, with pooled and study-specific columns.
#'   \item `errind`: performance summary for the individual study fits on the training layout.
#'   \item `erroverall`: performance summary for the pooled stage-one fit on the training layout.
#'   \item `fitoverall`: the shared overall fit reused across the alpha grid.
#'   \item `fitind`: the shared list of individual study fits reused across the alpha grid.
#'   \item `fit`: the list of `gptLasso` objects corresponding to `alpha.ptlasso.list`.
#'   \item `overall.lambda`, `ind.lambda`, `pre.lambda`: lambda rules used when summarizing candidate performance.
#' }
#'
#' @examples
#' # Gaussian cross-validation example
#' set.seed(1234)
#' sim_dat <- sim.gaussian.data()
#' x_train <- sim_dat$x_train
#'
#' cv_fit_gaussian <- cv.gptLasso(
#'   x = x_train,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   alpha.ptlasso.list = c(0, 0.5, 1),
#'   nfolds = 3
#' )
#'
#' cv_fit_gaussian$alpha.ptlasso.hat
#' cv_fit_gaussian$varying.alpha.ptlasso.hat
#' cv_fit_gaussian$errpre
#'
#' # Binomial cross-validation example
#' set.seed(5678)
#' sim_dat_bin <- sim.binary.data()
#' x_train_bin <- sim_dat_bin$x_train
#'
#' cv_fit_binomial <- cv.gptLasso(
#'   x = x_train_bin,
#'   family = "binomial",
#'   type.measure = "auc",
#'   alpha.ptlasso.list = c(0, 0.5, 1),
#'   nfolds = 3
#' )
#'
#' cv_fit_binomial$alpha.ptlasso.hat
#' head(cv_fit_binomial$errpre)
#' @export
cv.gptLasso <- function(
    x,
    alpha.ptlasso.list = seq(0, 1, length = 11),
    family = c("gaussian", "binomial"),
    type.measure = c("default", "mse", "auc", "deviance"),
    rho = seq(0, 1, length = 11),
    nfolds = 10,
    foldid = NULL,
    overall.lambda = c("lambda.min", "lambda.1se"),
    ind.lambda = c("lambda.min", "lambda.1se"),
    pre.lambda = c("lambda.min", "lambda.1se"),
    verbose = FALSE,
    fitoverall = NULL,
    fitind = NULL,
    group.intercepts = TRUE,
    parallel = FALSE,
    ncores = 1L,
    ...
) {
  this.call <- match.call()
  dot.args <- list(...)
  legacy.names <- c("alpha_ptlasso_list", "alpha_ptlasso_hat.choice", "s")
  bad.args <- intersect(names(dot.args), legacy.names)
  if (length(bad.args) > 0L) {
    stop(
      sprintf(
        "Legacy underscore argument(s) not supported in cv.gptLasso(): %s. Use dot-style names instead.",
        paste(bad.args, collapse = ", ")
      )
    )
  }
  
  family <- match.arg(family)
  type.measure <- match.arg(type.measure)
  if (type.measure == "default") {
    type.measure <- if (family == "gaussian") "mse" else "deviance"
  }
  overall.lambda <- match.arg(overall.lambda, c("lambda.min", "lambda.1se"))
  ind.lambda <- match.arg(ind.lambda, c("lambda.min", "lambda.1se"))
  pre.lambda <- match.arg(pre.lambda, c("lambda.min", "lambda.1se"))
  
  if (family == "binomial" && !(type.measure %in% c("auc", "deviance"))) {
    stop("For binomial family, type.measure must be 'auc' or 'deviance'.")
  }
  if (family == "gaussian" && !(type.measure %in% c("mse", "deviance"))) {
    stop("For gaussian family, type.measure must be 'mse' or 'deviance'.")
  }
  
  alpha.ptlasso.list <- sort(unique(as.numeric(alpha.ptlasso.list)))
  if (length(alpha.ptlasso.list) < 2L ||
      any(is.na(alpha.ptlasso.list)) ||
      any(alpha.ptlasso.list < 0 | alpha.ptlasso.list > 1)) {
    stop("alpha.ptlasso.list must contain at least two values between 0 and 1.")
  }
  rho <- unique(as.numeric(rho))
  if (length(rho) < 1L || any(is.na(rho))) {
    stop("rho must contain at least one numeric value.")
  }
  
  metric.rule <- ptmv_match_metric(type.measure)
  fit <- vector("list", length(alpha.ptlasso.list))
  err.rows <- vector("list", length(alpha.ptlasso.list))
  
  for (ii in seq_along(alpha.ptlasso.list)) {
    alpha.ptlasso <- alpha.ptlasso.list[ii]
    if (verbose) {
      message(sprintf("alpha.ptlasso = %s", format(alpha.ptlasso)))
    }
    
    fit[[ii]] <- gptLasso(
      x = x,
      alpha.ptlasso = alpha.ptlasso,
      family = family,
      type.measure = type.measure,
      rho = rho,
      overall.lambda = overall.lambda,
      ind.lambda = ind.lambda,
      pre.lambda = pre.lambda,
      nfolds = nfolds,
      foldid = foldid,
      verbose = verbose,
      fitoverall = fitoverall,
      fitind = fitind,
      group.intercepts = group.intercepts,
      parallel = parallel,
      ncores = ncores,
      ...
    )
    
    if (is.null(fitoverall)) fitoverall <- fit[[ii]]$fitoverall
    if (is.null(fitind)) fitind <- fit[[ii]]$fitind
    
    pred.pre <- lapply(seq_along(fit[[ii]]$fitpre), function(kk) {
      model <- fit[[ii]]$fitpre[[kk]]
      lambda <- ptmv_resolve_s(model, pre.lambda)
      lam.idx <- which.min(abs(model$lambda - lambda))
      as.numeric(model$fit.preval[, lam.idx])
    })
    names(pred.pre) <- fit[[ii]]$study.names
    
    study.err <- vapply(fit[[ii]]$study.names, function(study.name) {
      ptmv_metric_value(
        fit[[ii]]$training.layout$y[[study.name]],
        pred.pre[[study.name]],
        family,
        type.measure
      )
    }, numeric(1))
    
    err.rows[[ii]] <- c(
      overall = ptmv_metric_value(
        unlist(fit[[ii]]$training.layout$y, use.names = FALSE),
        unlist(pred.pre, use.names = FALSE),
        family = family,
        type.measure = type.measure
      ),
      stats::setNames(study.err, fit[[ii]]$study.names)
    )
  }
  
  errpre <- cbind(alpha.ptlasso = alpha.ptlasso.list, do.call(rbind, err.rows))
  rownames(errpre) <- NULL
  
  base.fit <- fit[[1]]
  overall.pred <- lapply(seq_along(base.fit$study.names), function(kk) {
    study.name <- base.fit$study.names[kk]
    lambda <- ptmv_resolve_s(base.fit$fitoverall, overall.lambda)
    preds <- predict(
      base.fit$fitoverall,
      newx = base.fit$training.layout$x[[study.name]],
      s = lambda,
      type = "response",
      newoffset = if (base.fit$group.intercepts) {
        rep(base.fit$group.baseline[study.name], base.fit$n.by.study[kk])
      } else {
        NULL
      }
    )
    as.numeric(preds)
  })
  names(overall.pred) <- base.fit$study.names
  
  ind.pred <- lapply(base.fit$study.names, function(study.name) {
    as.numeric(
      predict(
        base.fit$fitind[[study.name]],
        newx = base.fit$training.layout$x[[study.name]],
        s = ind.lambda,
        type = "response"
      )
    )
  })
  names(ind.pred) <- base.fit$study.names
  
  erroverall <- ptmv_summarize_metric(
    overall.pred,
    base.fit$training.layout$y,
    family,
    type.measure,
    add_r2 = family == "gaussian"
  )
  errind <- ptmv_summarize_metric(
    ind.pred,
    base.fit$training.layout$y,
    family,
    type.measure,
    add_r2 = family == "gaussian"
  )
  
  overall.idx <- metric.rule$best(errpre[, "overall"])
  alpha.ptlasso.hat <- alpha.ptlasso.list[overall.idx]
  selected.fit <- fit[[overall.idx]]
  varying.alpha.ptlasso.hat <- vapply(base.fit$study.names, function(study.name) {
    alpha.ptlasso.list[metric.rule$best(errpre[, study.name])]
  }, numeric(1))
  
  out <- list(
    call = this.call,
    alpha.ptlasso.hat = alpha.ptlasso.hat,
    varying.alpha.ptlasso.hat = varying.alpha.ptlasso.hat,
    alpha.ptlasso.list = alpha.ptlasso.list,
    rho = rho,
    errpre = errpre,
    errind = errind,
    erroverall = erroverall,
    fitoverall.rho = selected.fit$fitoverall.rho,
    fitind.rho = selected.fit$fitind.rho,
    fitpre.rho = selected.fit$fitpre.rho,
    fitoverall = fitoverall,
    fitind = fitind,
    fit = fit,
    family = family,
    type.measure = type.measure,
    overall.lambda = overall.lambda,
    ind.lambda = ind.lambda,
    pre.lambda = pre.lambda
  )
  class(out) <- "cv.gptLasso"
  out
}
