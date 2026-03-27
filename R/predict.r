#' Predict from a multistudy multiview gptLasso fit
#'
#' Generate overall, individual, and pretrained predictions from a `gptLasso`
#' object on Bioconductor-style multistudy multiview test data.
#'
#' @param object A fitted `gptLasso` object.
#' @param xtest A named list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#' @param ytest Optional list of response vectors used to compute performance summaries.
#' @param type Prediction scale: `"link"`, `"response"`, or `"class"` for binomial fits.
#' @param s Lambda rule used for prediction.
#' @param return.link Should link-scale predictions also be returned?
#' @param ... Reserved for future extensions.
#'
#' @return A `predict.gptLasso` object containing study-wise predictions
#'   for the overall, individual, and pretrained models, support summaries,
#'   and optional error summaries.
#' @export
predict.gptLasso <- function(object, xtest, ytest = NULL,
                             type = c("link", "response", "class"),
                             s = c("lambda.min", "lambda.1se"),
                             return.link = FALSE, ...) {
  if (missing(xtest)) {
    stop("Please supply xtest.")
  }
  if (!inherits(object, "gptLasso")) {
    stop("object must be a gptLasso fit.")
  }

  type <- match.arg(type)
  s <- match.arg(s)
  this.call <- match.call()

  if (object$family != "binomial" && type == "class") {
    stop("type = 'class' is only supported for binomial models.")
  }

  xtest_norm <- ptmv_normalize_newdata(xtest, object$training_layout)
  xtest_split <- xtest_norm$x
  test_study_names <- xtest_norm$study_names
  test_sizes <- xtest_norm$n_by_study
  ytest <- ptmv_prepare_ytest(ytest, test_study_names, study_sizes = test_sizes)

  x_list_test <- ptmv_stack_by_view(xtest_split, object$view_names)
  overall_baseline_offset <- NULL
  if (isTRUE(object$group.intercepts)) {
    overall_baseline_offset <- ptmv_study_offsets(test_sizes, object$group_baseline[test_study_names])
  }

  overall_link_stacked <- as.numeric(predict(
    object$fitoverall,
    newx = x_list_test,
    s = s,
    type = "link",
    newoffset = overall_baseline_offset
  ))
  overall_resp_stacked <- if (type == "class") {
    as.numeric(predict(
      object$fitoverall,
      newx = x_list_test,
      s = s,
      type = "response",
      newoffset = overall_baseline_offset
    ))
  } else {
    as.numeric(predict(
      object$fitoverall,
      newx = x_list_test,
      s = s,
      type = type,
      newoffset = overall_baseline_offset
    ))
  }

  overall_link <- ptmv_split_vector_by_study(overall_link_stacked, test_sizes, test_study_names)
  overall_resp <- ptmv_split_vector_by_study(overall_resp_stacked, test_sizes, test_study_names)
  overall_pred <- if (type == "class") {
    lapply(overall_resp, function(x) ifelse(x >= 0.5, 1, 0))
  } else {
    overall_resp
  }

  stage1_link <- ptmv_split_vector_by_study(
    as.numeric(predict(
      object$fitoverall,
      newx = x_list_test,
      s = object$fitoverall.lambda,
      type = "link",
      newoffset = overall_baseline_offset
    )),
    test_sizes,
    test_study_names
  )

  pre_offsets <- lapply(stage1_link, function(x) (1 - object$alpha_ptlasso) * x)
  linkpre <- ptmv_predict_by_study(object$fitpre, xtest_split, s = s, type = "link", offsets = pre_offsets)
  yhatpre <- if (type == "class") {
    ptmv_predict_by_study(object$fitpre, xtest_split, s = s, type = "response", offsets = pre_offsets)
  } else {
    ptmv_predict_by_study(object$fitpre, xtest_split, s = s, type = type, offsets = pre_offsets)
  }
  if (type == "class") {
    yhatpre <- lapply(yhatpre, function(x) ifelse(x >= 0.5, 1, 0))
  }

  linkind <- ptmv_predict_by_study(object$fitind, xtest_split, s = s, type = "link")
  yhatind <- if (type == "class") {
    ptmv_predict_by_study(object$fitind, xtest_split, s = s, type = "response")
  } else {
    ptmv_predict_by_study(object$fitind, xtest_split, s = s, type = type)
  }
  if (type == "class") {
    yhatind <- lapply(yhatind, function(x) ifelse(x >= 0.5, 1, 0))
  }

  supoverall <- ptmv_get_support(object$fitoverall, s)
  supind <- ptmv_get_union_support(object$fitind, s)
  suppre.common <- ptmv_get_support(object$fitoverall, object$fitoverall.lambda)
  suppre.individual <- setdiff(ptmv_get_union_support(object$fitpre, s), suppre.common)

  erroverall <- errind <- errpre <- NULL
  if (!is.null(ytest)) {
    erroverall <- ptmv_summarize_metric(overall_resp, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    errind <- ptmv_summarize_metric(yhatind, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    errpre <- ptmv_summarize_metric(yhatpre, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
  }

  ptmv_build_prediction_object(
    call = this.call,
    alpha_ptlasso = object$alpha_ptlasso,
    type.measure = object$type.measure,
    yhatoverall = overall_pred,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = suppre.individual,
    linkoverall = if (return.link) overall_link else NULL,
    linkind = if (return.link) linkind else NULL,
    linkpre = if (return.link) linkpre else NULL,
    erroverall = erroverall,
    errind = errind,
    errpre = errpre,
    class_name = "predict.gptLasso"
  )
}

#' Predict from a cross-validated multistudy multiview gptLasso fit
#'
#' Resolve one fixed or varying transfer-learning choice from a `cv.gptLasso`
#' object and generate predictions on Bioconductor-style multistudy multiview
#' test data.
#'
#' @param object A fitted `cv.gptLasso` object.
#' @param xtest A named list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#' @param ytest Optional list of response vectors used to compute performance summaries.
#' @param alpha_ptlasso Optional user-specified transfer-learning choice. May be one value or one per study.
#' @param alpha_ptlasso_type Either `"fixed"` or `"varying"` when `alpha_ptlasso` is not supplied.
#' @param type Prediction scale: `"link"`, `"response"`, or `"class"` for binomial fits.
#' @param s Lambda rule used for prediction.
#' @param return.link Should link-scale predictions also be returned?
#' @param ... Reserved for future extensions.
#'
#' @return A `predict.cv.gptLasso` object containing the chosen transfer-learning
#'   setting, study-wise predictions, support summaries, and optional error summaries.
#' @export
predict.cv.gptLasso <- function(object, xtest, ytest = NULL,
                                alpha_ptlasso = NULL,
                                alpha_ptlasso_type = c("fixed", "varying"),
                                type = c("link", "response", "class"),
                                s = c("lambda.min", "lambda.1se"),
                                return.link = FALSE, ...) {
  if (missing(xtest)) {
    stop("Please supply xtest.")
  }
  if (!inherits(object, "cv.gptLasso")) {
    stop("object must be a cv.gptLasso fit.")
  }

  this.call <- match.call()
  alpha_ptlasso_type <- match.arg(alpha_ptlasso_type)
  type <- match.arg(type)
  s <- match.arg(s)

  close.enough <- 1e-6
  if (is.null(alpha_ptlasso)) {
    alpha_ptlasso <- if (alpha_ptlasso_type == "fixed") {
      object$alpha_ptlasso_hat
    } else {
      object$varying.alpha_ptlasso_hat
    }
  }

  if (length(alpha_ptlasso) == 1L) {
    model_idx <- which(abs(object$alpha_ptlasso_list - alpha_ptlasso) < close.enough)
    if (length(model_idx) == 0L) {
      stop("Not a valid choice of alpha_ptlasso. Please choose from object$alpha_ptlasso_list.")
    }
    fit <- object$fit[[model_idx[1]]]
    out <- predict.gptLasso(
      fit,
      xtest = xtest,
      ytest = ytest,
      type = type,
      s = s,
      return.link = return.link,
      ...
    )
    out$call <- this.call
    out$fit <- object
    class(out) <- c("predict.cv.gptLasso", "predict.gptLasso")
    return(out)
  }

  if (is.null(names(alpha_ptlasso))) {
    if (length(alpha_ptlasso) != length(object$fit[[1]]$study_names)) {
      stop("Must have one alpha_ptlasso for each study.")
    }
    names(alpha_ptlasso) <- object$fit[[1]]$study_names
  }
  if (!setequal(names(alpha_ptlasso), object$fit[[1]]$study_names)) {
    stop("alpha_ptlasso vector names must match the training study names exactly.")
  }
  alpha_ptlasso <- alpha_ptlasso[object$fit[[1]]$study_names]

  if (!all(vapply(alpha_ptlasso, function(a) any(abs(object$alpha_ptlasso_list - a) < close.enough), logical(1)))) {
    stop("Includes at least one invalid alpha_ptlasso choice. Please choose from object$alpha_ptlasso_list.")
  }

  pred_by_alpha <- lapply(object$fit, function(fit) {
    predict.gptLasso(
      fit,
      xtest = xtest,
      ytest = ytest,
      type = type,
      s = s,
      return.link = return.link,
      ...
    )
  })
  names(pred_by_alpha) <- as.character(object$alpha_ptlasso_list)

  yhatpre <- yhatind <- yhatoverall <- vector("list", length(alpha_ptlasso))
  names(yhatpre) <- names(yhatind) <- names(yhatoverall) <- names(alpha_ptlasso)
  if (return.link) {
    linkpre <- linkind <- linkoverall <- vector("list", length(alpha_ptlasso))
    names(linkpre) <- names(linkind) <- names(linkoverall) <- names(alpha_ptlasso)
  } else {
    linkpre <- linkind <- linkoverall <- NULL
  }

  for (study_name in names(alpha_ptlasso)) {
    key <- as.character(alpha_ptlasso[[study_name]])
    pred <- pred_by_alpha[[key]]
    yhatpre[[study_name]] <- pred$yhatpre[[study_name]]
    yhatind[[study_name]] <- pred$yhatind[[study_name]]
    yhatoverall[[study_name]] <- pred$yhatoverall[[study_name]]
    if (return.link) {
      linkpre[[study_name]] <- pred$linkpre[[study_name]]
      linkind[[study_name]] <- pred$linkind[[study_name]]
      linkoverall[[study_name]] <- pred$linkoverall[[study_name]]
    }
  }

  supoverall <- pred_by_alpha[[as.character(object$alpha_ptlasso_hat)]]$supoverall
  supind <- sort(unique(unlist(lapply(names(alpha_ptlasso), function(study_name) {
    pred_by_alpha[[as.character(alpha_ptlasso[[study_name]])]]$supind
  }))))
  suppre.common <- pred_by_alpha[[as.character(object$alpha_ptlasso_hat)]]$suppre.common
  suppre.individual <- sort(unique(unlist(lapply(names(alpha_ptlasso), function(study_name) {
    pred_by_alpha[[as.character(alpha_ptlasso[[study_name]])]]$suppre.individual
  }))))

  erroverall <- errind <- errpre <- NULL
  if (!is.null(ytest)) {
    xtest_norm <- ptmv_normalize_newdata(xtest, object$fit[[1]]$training_layout)
    ytest_norm <- ptmv_prepare_ytest(ytest, names(alpha_ptlasso), study_sizes = xtest_norm$n_by_study)
    family <- object$family
    type.measure <- object$type.measure
    erroverall <- ptmv_summarize_metric(yhatoverall, ytest_norm, family, type.measure, add_r2 = family == "gaussian")
    errind <- ptmv_summarize_metric(yhatind, ytest_norm, family, type.measure, add_r2 = family == "gaussian")
    errpre <- ptmv_summarize_metric(yhatpre, ytest_norm, family, type.measure, add_r2 = family == "gaussian")
  }

  ptmv_build_prediction_object(
    call = this.call,
    alpha_ptlasso = alpha_ptlasso,
    type.measure = object$type.measure,
    yhatoverall = yhatoverall,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = setdiff(suppre.individual, suppre.common),
    linkoverall = linkoverall,
    linkind = linkind,
    linkpre = linkpre,
    erroverall = erroverall,
    errind = errind,
    errpre = errpre,
    fit = object,
    class_name = "predict.cv.gptLasso"
  )
}
