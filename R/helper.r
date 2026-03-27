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
