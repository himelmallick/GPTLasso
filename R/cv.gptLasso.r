#' Cross-validation for multistudy multiview gptLasso
#'
#' Run `gptLasso()` across a grid of transfer-learning `alpha_ptlasso` values
#' and select the best value using held-out performance.
#'
#' @param x A named list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#' @param alpha_ptlasso_list Numeric vector of transfer-learning values to compare.
#' @param family Response family. Currently `"gaussian"` and `"binomial"` are supported.
#' @param type.measure Cross-validation metric used to compare transfer-learning levels.
#' @param nfolds Number of folds used inside each `gptLasso()` fit.
#' @param foldid Optional stacked fold assignment across all studies.
#' @param s Lambda rule used when summarizing CV performance.
#' @param alpha_ptlasso_hat.choice Whether to choose the final alpha using overall or mean performance.
#' @param verbose Should progress messages be printed?
#' @param fitoverall Optional pre-fit overall multiview model reused across alphas.
#' @param fitind Optional pre-fit list of individual multiview models reused across alphas.
#' @param group.intercepts Should study-specific stage-one baselines be used?
#' @param parallel Logical; if `TRUE`, allow study-level parallel fits where available.
#' @param ncores Number of worker cores for study-level parallel fits.
#' @param ... Additional arguments forwarded to `gptLasso()`.
#'
#' @return A `cv.gptLasso` object containing the selected transfer-learning
#'   values, performance summaries, and the fitted model list.
#' @export
cv.gptLasso <- function(
    x,
    alpha_ptlasso_list = seq(0, 1, length = 11),
    family = c("gaussian", "binomial"),
    type.measure = c("default", "mse", "auc", "deviance"),
    nfolds = 10,
    foldid = NULL,
    s = c("lambda.min", "lambda.1se"),
    alpha_ptlasso_hat.choice = c("overall", "mean"),
    verbose = FALSE,
    fitoverall = NULL,
    fitind = NULL,
    group.intercepts = TRUE,
    parallel = FALSE,
    ncores = 1L,
    ...
) {
  this.call <- match.call()
  family <- match.arg(family)
  type.measure <- match.arg(type.measure)
  if (type.measure == "default") {
    type.measure <- if (family == "gaussian") "mse" else "deviance"
  }
  s <- match.arg(s)
  alpha_ptlasso_hat.choice <- match.arg(alpha_ptlasso_hat.choice)

  if (family == "binomial" && !(type.measure %in% c("auc", "deviance"))) {
    stop("For binomial family, type.measure must be 'auc' or 'deviance'.")
  }
  if (family == "gaussian" && !(type.measure %in% c("mse", "deviance"))) {
    stop("For gaussian family, type.measure must be 'mse' or 'deviance'.")
  }

  alpha_ptlasso_list <- sort(unique(as.numeric(alpha_ptlasso_list)))
  if (length(alpha_ptlasso_list) < 2L ||
      any(is.na(alpha_ptlasso_list)) ||
      any(alpha_ptlasso_list < 0 | alpha_ptlasso_list > 1)) {
    stop("alpha_ptlasso_list must contain at least two values between 0 and 1.")
  }

  metric_rule <- ptmv_match_metric(type.measure)
  fit <- vector("list", length(alpha_ptlasso_list))
  err_rows <- vector("list", length(alpha_ptlasso_list))

  for (ii in seq_along(alpha_ptlasso_list)) {
    alpha_ptlasso <- alpha_ptlasso_list[ii]
    if (verbose) {
      message(sprintf("alpha_ptlasso = %s", format(alpha_ptlasso)))
    }

    fit[[ii]] <- gptLasso(
      x = x,
      alpha_ptlasso = alpha_ptlasso,
      family = family,
      type.measure = type.measure,
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

    pred_pre <- lapply(seq_along(fit[[ii]]$fitpre), function(kk) {
      model <- fit[[ii]]$fitpre[[kk]]
      lambda <- ptmv_resolve_s(model, s)
      lam_idx <- which.min(abs(model$lambda - lambda))
      as.numeric(model$fit.preval[, lam_idx])
    })
    names(pred_pre) <- fit[[ii]]$study_names

    err_row <- c(
      overall = ptmv_metric_value(
        unlist(fit[[ii]]$training_layout$y, use.names = FALSE),
        unlist(pred_pre, use.names = FALSE),
        family = family,
        type.measure = type.measure
      ),
      mean = mean(vapply(fit[[ii]]$study_names, function(study_name) {
        ptmv_metric_value(fit[[ii]]$training_layout$y[[study_name]], pred_pre[[study_name]], family, type.measure)
      }, numeric(1)), na.rm = TRUE),
      setNames(
        vapply(fit[[ii]]$study_names, function(study_name) {
          ptmv_metric_value(fit[[ii]]$training_layout$y[[study_name]], pred_pre[[study_name]], family, type.measure)
        }, numeric(1)),
        paste0("group_", fit[[ii]]$study_names)
      )
    )
    err_rows[[ii]] <- err_row
  }

  errpre <- cbind(alpha_ptlasso = alpha_ptlasso_list, do.call(rbind, err_rows))
  rownames(errpre) <- NULL

  base_fit <- fit[[1]]
  overall_pred <- lapply(seq_along(base_fit$study_names), function(kk) {
    study_name <- base_fit$study_names[kk]
    lambda <- ptmv_resolve_s(base_fit$fitoverall, s)
    preds <- predict(
      base_fit$fitoverall,
      newx = base_fit$training_layout$x[[study_name]],
      s = lambda,
      type = "response",
      newoffset = if (base_fit$group.intercepts) rep(base_fit$group_baseline[study_name], base_fit$n_by_study[kk]) else NULL
    )
    as.numeric(preds)
  })
  names(overall_pred) <- base_fit$study_names

  ind_pred <- lapply(base_fit$study_names, function(study_name) {
    as.numeric(predict(base_fit$fitind[[study_name]], newx = base_fit$training_layout$x[[study_name]], s = s, type = "response"))
  })
  names(ind_pred) <- base_fit$study_names

  erroverall <- ptmv_summarize_metric(overall_pred, base_fit$training_layout$y, family, type.measure, add_r2 = family == "gaussian")
  errind <- ptmv_summarize_metric(ind_pred, base_fit$training_layout$y, family, type.measure, add_r2 = family == "gaussian")

  overall_idx <- if (alpha_ptlasso_hat.choice == "mean") {
    metric_rule$best(errpre[, "mean"])
  } else {
    metric_rule$best(errpre[, "overall"])
  }
  alpha_ptlasso_hat <- alpha_ptlasso_list[overall_idx]
  varying.alpha_ptlasso_hat <- vapply(base_fit$study_names, function(study_name) {
    alpha_ptlasso_list[metric_rule$best(errpre[, paste0("group_", study_name)])]
  }, numeric(1))

  out <- list(
    call = this.call,
    alpha_ptlasso_hat = alpha_ptlasso_hat,
    varying.alpha_ptlasso_hat = varying.alpha_ptlasso_hat,
    alpha_ptlasso_list = alpha_ptlasso_list,
    errpre = errpre,
    errind = errind,
    erroverall = erroverall,
    fitoverall = fitoverall,
    fitind = fitind,
    fit = fit,
    family = family,
    type.measure = type.measure,
    s = s,
    alpha_ptlasso_hat.choice = alpha_ptlasso_hat.choice
  )
  class(out) <- "cv.gptLasso"
  out
}
