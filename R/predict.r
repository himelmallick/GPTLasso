#' Predict from a fitted multistudy multiview gptLasso model
#'
#' Generate held-out predictions from the overall, study-specific individual,
#' and pretrained transfer models stored in a `gptLasso` fit.
#'
#' @param object A fitted `gptLasso` object.
#' @param xtest A named list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#' @param ytest Optional list of response vectors used to compute performance summaries.
#' @param type Prediction scale: `"link"`, `"response"`, or `"class"` for binomial fits.
#' @param overall.lambda Lambda rule used for overall-model prediction.
#' @param ind.lambda Lambda rule used for individual-model prediction.
#' @param pre.lambda Lambda rule used for pretrained-model prediction.
#' @param return.link Should link-scale predictions also be returned?
#' @param ... Reserved for future extensions.
#'
#' @return A `predict.gptLasso` object, stored as a list. Important components include:
#' \itemize{
#'   \item `yhatoverall`: study-wise predictions from the pooled stage-one fit.
#'   \item `yhatind`: study-wise predictions from the individual study fits.
#'   \item `yhatpre`: study-wise predictions from the pretrained transfer fits.
#'   \item `supoverall`: selected feature indices from the pooled stage-one fit.
#'   \item `supind`: union of selected feature indices across the individual fits.
#'   \item `suppre.common`: stage-one support reused by the pretrained fits.
#'   \item `suppre.individual`: additional feature indices selected by the pretrained fits beyond `suppre.common`.
#'   \item `linkoverall`, `linkind`, `linkpre`: optional link-scale predictions when `return.link = TRUE`.
#'   \item `metrics`: optional held-out performance summaries when `ytest` is supplied. The primary metric is stored in `metrics$MSE`, `metrics$AUC`, or `metrics$deviance` as a `3 x (k + 1)` table with rows `overall`, `ind`, and `pre`, and columns `mean` plus one column per study. Gaussian fits also include `metrics$r2` with the same layout.
#'   \item `erroverall`, `errind`, `errpre`: backward-compatible performance summaries when `ytest` is supplied.
#' }
#'
#' @examples
#' # Gaussian prediction example
#' set.seed(1234)
#' sim_dat <- sim.gaussian.data()
#' x_train <- sim_dat$x_train
#' x_test <- sim_dat$x_test
#' ytest <- split(x_test$sample_metadata$Y, x_test$sample_metadata$study)
#'
#' fit_gaussian <- gptLasso(
#'   x = x_train,
#'   alpha.ptlasso = 0.5,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   nfolds = 3
#' )
#'
#' pred_gaussian <- predict(
#'   fit_gaussian,
#'   xtest = x_test,
#'   ytest = ytest,
#'   type = "response"
#' )
#'
#' names(pred_gaussian)
#' str(pred_gaussian$yhatpre)
#' pred_gaussian$errpre
#'
#' # Binomial prediction example
#' set.seed(5678)
#' sim_dat_bin <- sim.binary.data()
#' x_train_bin <- sim_dat_bin$x_train
#' x_test_bin <- sim_dat_bin$x_test
#' ytest_bin <- split(x_test_bin$sample_metadata$Y, x_test_bin$sample_metadata$study)
#'
#' fit_binomial <- gptLasso(
#'   x = x_train_bin,
#'   alpha.ptlasso = 0.5,
#'   family = "binomial",
#'   type.measure = "auc",
#'   nfolds = 3
#' )
#'
#' pred_binomial <- predict(
#'   fit_binomial,
#'   xtest = x_test_bin,
#'   ytest = ytest_bin,
#'   type = "response"
#' )
#'
#' pred_binomial$errpre
#' @export
predict.gptLasso <- function(object, xtest, ytest = NULL,
                             type = c("link", "response", "class"),
                             overall.lambda = c("lambda.min", "lambda.1se"),
                             ind.lambda = c("lambda.min", "lambda.1se"),
                             pre.lambda = c("lambda.min", "lambda.1se"),
                             return.link = FALSE, ...) {
  if (missing(xtest)) {
    stop("Please supply xtest.")
  }
  if (!inherits(object, "gptLasso")) {
    stop("object must be a gptLasso fit.")
  }
  
  dot.args <- list(...)
  legacy.names <- c("s")
  bad.args <- intersect(names(dot.args), legacy.names)
  if (length(bad.args) > 0L) {
    stop(
      sprintf(
        "Legacy underscore argument(s) not supported in predict.gptLasso(): %s. Use dot-style names instead.",
        paste(bad.args, collapse = ", ")
      )
    )
  }
  
  type <- match.arg(type)
  overall.lambda <- match.arg(overall.lambda, c("lambda.min", "lambda.1se"))
  ind.lambda <- match.arg(ind.lambda, c("lambda.min", "lambda.1se"))
  pre.lambda <- match.arg(pre.lambda, c("lambda.min", "lambda.1se"))
  this.call <- match.call()
  
  if (object$family != "binomial" && type == "class") {
    stop("type = 'class' is only supported for binomial models.")
  }
  
  xtest.norm <- ptmv_normalize_newdata(xtest, object$training.layout)
  xtest.split <- xtest.norm$x
  test.study.names <- xtest.norm$study_names
  test.sizes <- xtest.norm$n_by_study
  ytest <- ptmv_prepare_ytest(ytest, test.study.names, study_sizes = test.sizes)
  
  x.list.test <- ptmv_stack_by_view(xtest.split, object$view.names)
  overall.baseline.offset <- NULL
  if (isTRUE(object$group.intercepts)) {
    overall.baseline.offset <- ptmv_study_offsets(test.sizes, object$group.baseline[test.study.names])
  }
  
  overall.lambda.value <- ptmv_resolve_s(object$fitoverall, overall.lambda)
  overall.link.stacked <- as.numeric(predict(
    object$fitoverall,
    newx = x.list.test,
    s = overall.lambda.value,
    type = "link",
    newoffset = overall.baseline.offset
  ))
  overall.resp.stacked <- if (type == "class") {
    as.numeric(predict(
      object$fitoverall,
      newx = x.list.test,
      s = overall.lambda.value,
      type = "response",
      newoffset = overall.baseline.offset
    ))
  } else {
    as.numeric(predict(
      object$fitoverall,
      newx = x.list.test,
      s = overall.lambda.value,
      type = type,
      newoffset = overall.baseline.offset
    ))
  }
  
  overall.link <- ptmv_split_vector_by_study(overall.link.stacked, test.sizes, test.study.names)
  overall.resp <- ptmv_split_vector_by_study(overall.resp.stacked, test.sizes, test.study.names)
  overall.pred <- if (type == "class") {
    lapply(overall.resp, function(x) ifelse(x >= 0.5, 1, 0))
  } else {
    overall.resp
  }
  
  stage1.link <- ptmv_split_vector_by_study(
    as.numeric(predict(
      object$fitoverall,
      newx = x.list.test,
      s = object$fitoverall.lambda,
      type = "link",
      newoffset = overall.baseline.offset
    )),
    test.sizes,
    test.study.names
  )
  
  pre.offsets <- lapply(stage1.link, function(x) (1 - object$alpha.ptlasso) * x)
  linkpre <- ptmv_predict_by_study(object$fitpre, xtest.split, s = pre.lambda, type = "link", offsets = pre.offsets)
  yhatpre.resp <- ptmv_predict_by_study(object$fitpre, xtest.split, s = pre.lambda, type = "response", offsets = pre.offsets)
  yhatpre <- if (type == "class") {
    yhatpre.resp
  } else {
    ptmv_predict_by_study(object$fitpre, xtest.split, s = pre.lambda, type = type, offsets = pre.offsets)
  }
  if (type == "class") {
    yhatpre <- lapply(yhatpre, function(x) ifelse(x >= 0.5, 1, 0))
  }
  
  linkind <- ptmv_predict_by_study(object$fitind, xtest.split, s = ind.lambda, type = "link")
  yhatind.resp <- ptmv_predict_by_study(object$fitind, xtest.split, s = ind.lambda, type = "response")
  yhatind <- if (type == "class") {
    yhatind.resp
  } else {
    ptmv_predict_by_study(object$fitind, xtest.split, s = ind.lambda, type = type)
  }
  if (type == "class") {
    yhatind <- lapply(yhatind, function(x) ifelse(x >= 0.5, 1, 0))
  }
  
  supoverall <- ptmv_get_support(object$fitoverall, overall.lambda)
  supind <- ptmv_get_union_support(object$fitind, ind.lambda)
  suppre.common <- ptmv_get_support(object$fitoverall, object$fitoverall.lambda)
  suppre.individual <- setdiff(ptmv_get_union_support(object$fitpre, pre.lambda), suppre.common)
  
  metrics <- erroverall <- errind <- errpre <- NULL
  if (!is.null(ytest)) {
    erroverall <- ptmv_summarize_metric(overall.resp, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    errind <- ptmv_summarize_metric(yhatind.resp, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    errpre <- ptmv_summarize_metric(yhatpre.resp, ytest, object$family, object$type.measure, add_r2 = object$family == "gaussian")
    metrics <- ptmv_build_metric_report(
      overall_preds = overall.resp,
      ind_preds = yhatind.resp,
      pre_preds = yhatpre.resp,
      y = ytest,
      family = object$family,
      type.measure = object$type.measure
    )
  }
  
  ptmv_build_prediction_object(
    call = this.call,
    alpha.ptlasso = object$alpha.ptlasso,
    type.measure = object$type.measure,
    yhatoverall = overall.pred,
    yhatind = yhatind,
    yhatpre = yhatpre,
    supoverall = supoverall,
    supind = supind,
    suppre.common = suppre.common,
    suppre.individual = suppre.individual,
    linkoverall = if (return.link) overall.link else NULL,
    linkind = if (return.link) linkind else NULL,
    linkpre = if (return.link) linkpre else NULL,
    metrics = metrics,
    erroverall = erroverall,
    errind = errind,
    errpre = errpre,
    metric_predictions = if (!is.null(ytest)) {
      list(overall = overall.resp, ind = yhatind.resp, pre = yhatpre.resp)
    } else {
      NULL
    },
    class_name = "predict.gptLasso"
  )
}

#' Predict from a cross-validated multistudy multiview gptLasso fit
#'
#' Resolve one fixed or varying transfer-learning choice from a `cv.gptLasso`
#' object and generate held-out predictions from the corresponding overall,
#' individual, and pretrained models.
#'
#' @param object A fitted `cv.gptLasso` object.
#' @param xtest A named list with `feature_table`, `sample_metadata`, and `feature_metadata`.
#' @param ytest Optional list of response vectors used to compute performance summaries.
#' @param alpha.ptlasso Optional user-specified transfer-learning choice. May be one value or one per study.
#' @param alpha.ptlasso.type Either `"fixed"` or `"varying"` when `alpha.ptlasso` is not supplied.
#' @param type Prediction scale: `"link"`, `"response"`, or `"class"` for binomial fits.
#' @param overall.lambda Lambda rule used for overall-model prediction.
#' @param ind.lambda Lambda rule used for individual-model prediction.
#' @param pre.lambda Lambda rule used for pretrained-model prediction.
#' @param return.link Should link-scale predictions also be returned?
#' @param ... Reserved for future extensions.
#'
#' @return A `predict.cv.gptLasso` object with the same prediction components as
#'   `predict.gptLasso()`, plus:
#' \itemize{
#'   \item `alpha.ptlasso`: the chosen transfer-learning value, either fixed or study-specific.
#'   \item `fit`: the originating `cv.gptLasso` object.
#' }
#' When `ytest` is supplied, the returned object also includes `metrics`,
#' `erroverall`, `errind`, and `errpre`. Each metric table has rows `overall`,
#' `ind`, and `pre`, and columns `mean` plus one column per study.
#'
#' @examples
#' # Gaussian cross-validated prediction example
#' set.seed(1234)
#' sim_dat <- sim.gaussian.data()
#' x_train <- sim_dat$x_train
#' x_test <- sim_dat$x_test
#' ytest <- split(x_test$sample_metadata$Y, x_test$sample_metadata$study)
#'
#' cv_fit_gaussian <- cv.gptLasso(
#'   x = x_train,
#'   family = "gaussian",
#'   type.measure = "mse",
#'   alpha.ptlasso.list = c(0, 0.5, 1),
#'   nfolds = 3
#' )
#'
#' pred_cv_gaussian <- predict(
#'   cv_fit_gaussian,
#'   xtest = x_test,
#'   ytest = ytest,
#'   alpha.ptlasso.type = "fixed",
#'   type = "response"
#' )
#'
#' pred_cv_gaussian$alpha.ptlasso
#' pred_cv_gaussian$errpre
#'
#' # Binomial cross-validated prediction example
#' set.seed(5678)
#' sim_dat_bin <- sim.binary.data()
#' x_train_bin <- sim_dat_bin$x_train
#' x_test_bin <- sim_dat_bin$x_test
#' ytest_bin <- split(x_test_bin$sample_metadata$Y, x_test_bin$sample_metadata$study)
#'
#' cv_fit_binomial <- cv.gptLasso(
#'   x = x_train_bin,
#'   family = "binomial",
#'   type.measure = "auc",
#'   alpha.ptlasso.list = c(0, 0.5, 1),
#'   nfolds = 3
#' )
#'
#' pred_cv_binomial <- predict(
#'   cv_fit_binomial,
#'   xtest = x_test_bin,
#'   ytest = ytest_bin,
#'   alpha.ptlasso.type = "fixed",
#'   type = "response"
#' )
#'
#' pred_cv_binomial$fit$alpha.ptlasso.hat
#' @export
predict.cv.gptLasso <- function(object, xtest, ytest = NULL,
                                alpha.ptlasso = NULL,
                                alpha.ptlasso.type = c("fixed", "varying"),
                                type = c("link", "response", "class"),
                                overall.lambda = c("lambda.min", "lambda.1se"),
                                ind.lambda = c("lambda.min", "lambda.1se"),
                                pre.lambda = c("lambda.min", "lambda.1se"),
                                return.link = FALSE, ...) {
  if (missing(xtest)) {
    stop("Please supply xtest.")
  }
  if (!inherits(object, "cv.gptLasso")) {
    stop("object must be a cv.gptLasso fit.")
  }
  
  dot.args <- list(...)
  legacy.names <- c("alpha_ptlasso", "alpha_ptlasso_type", "s")
  bad.args <- intersect(names(dot.args), legacy.names)
  if (length(bad.args) > 0L) {
    stop(
      sprintf(
        "Legacy underscore argument(s) not supported in predict.cv.gptLasso(): %s. Use dot-style names instead.",
        paste(bad.args, collapse = ", ")
      )
    )
  }
  
  this.call <- match.call()
  alpha.ptlasso.type <- match.arg(alpha.ptlasso.type)
  type <- match.arg(type)
  overall.lambda <- match.arg(overall.lambda, c("lambda.min", "lambda.1se"))
  ind.lambda <- match.arg(ind.lambda, c("lambda.min", "lambda.1se"))
  pre.lambda <- match.arg(pre.lambda, c("lambda.min", "lambda.1se"))
  
  close.enough <- 1e-6
  if (is.null(alpha.ptlasso)) {
    alpha.ptlasso <- if (alpha.ptlasso.type == "fixed") {
      object$alpha.ptlasso.hat
    } else {
      object$varying.alpha.ptlasso.hat
    }
  }
  
  if (length(alpha.ptlasso) == 1L) {
    model.idx <- which(abs(object$alpha.ptlasso.list - alpha.ptlasso) < close.enough)
    if (length(model.idx) == 0L) {
      stop("Not a valid choice of alpha.ptlasso. Please choose from object$alpha.ptlasso.list.")
    }
    fit <- object$fit[[model.idx[1]]]
    out <- predict.gptLasso(
      fit,
      xtest = xtest,
      ytest = ytest,
      type = type,
      overall.lambda = overall.lambda,
      ind.lambda = ind.lambda,
      pre.lambda = pre.lambda,
      return.link = return.link,
      ...
    )
    out$call <- this.call
    out$fit <- object
    class(out) <- c("predict.cv.gptLasso", "predict.gptLasso")
    return(out)
  }
  
  if (is.null(names(alpha.ptlasso))) {
    if (length(alpha.ptlasso) != length(object$fit[[1]]$study.names)) {
      stop("Must have one alpha.ptlasso for each study.")
    }
    names(alpha.ptlasso) <- object$fit[[1]]$study.names
  }
  if (!setequal(names(alpha.ptlasso), object$fit[[1]]$study.names)) {
    stop("alpha.ptlasso vector names must match the training study names exactly.")
  }
  alpha.ptlasso <- alpha.ptlasso[object$fit[[1]]$study.names]
  
  if (!all(vapply(alpha.ptlasso, function(a) any(abs(object$alpha.ptlasso.list - a) < close.enough), logical(1)))) {
    stop("Includes at least one invalid alpha.ptlasso choice. Please choose from object$alpha.ptlasso.list.")
  }
  
  pred.by.alpha <- lapply(object$fit, function(fit) {
    predict.gptLasso(
      fit,
      xtest = xtest,
      ytest = ytest,
      type = type,
      overall.lambda = overall.lambda,
      ind.lambda = ind.lambda,
      pre.lambda = pre.lambda,
      return.link = return.link,
      ...
    )
  })
  names(pred.by.alpha) <- as.character(object$alpha.ptlasso.list)
  
  yhatpre <- yhatind <- yhatoverall <- vector("list", length(alpha.ptlasso))
  names(yhatpre) <- names(yhatind) <- names(yhatoverall) <- names(alpha.ptlasso)
  if (return.link) {
    linkpre <- linkind <- linkoverall <- vector("list", length(alpha.ptlasso))
    names(linkpre) <- names(linkind) <- names(linkoverall) <- names(alpha.ptlasso)
  } else {
    linkpre <- linkind <- linkoverall <- NULL
  }
  
  for (study.name in names(alpha.ptlasso)) {
    key <- as.character(alpha.ptlasso[[study.name]])
    pred <- pred.by.alpha[[key]]
    yhatpre[[study.name]] <- pred$yhatpre[[study.name]]
    yhatind[[study.name]] <- pred$yhatind[[study.name]]
    yhatoverall[[study.name]] <- pred$yhatoverall[[study.name]]
    if (return.link) {
      linkpre[[study.name]] <- pred$linkpre[[study.name]]
      linkind[[study.name]] <- pred$linkind[[study.name]]
      linkoverall[[study.name]] <- pred$linkoverall[[study.name]]
    }
  }
  
  supoverall <- pred.by.alpha[[as.character(object$alpha.ptlasso.hat)]]$supoverall
  supind <- sort(unique(unlist(lapply(names(alpha.ptlasso), function(study.name) {
    pred.by.alpha[[as.character(alpha.ptlasso[[study.name]])]]$supind
  }))))
  suppre.common <- pred.by.alpha[[as.character(object$alpha.ptlasso.hat)]]$suppre.common
  suppre.individual <- sort(unique(unlist(lapply(names(alpha.ptlasso), function(study.name) {
    pred.by.alpha[[as.character(alpha.ptlasso[[study.name]])]]$suppre.individual
  }))))
  
  metrics <- erroverall <- errind <- errpre <- NULL
  if (!is.null(ytest)) {
    xtest.norm <- ptmv_normalize_newdata(xtest, object$fit[[1]]$training.layout)
    ytest.norm <- ptmv_prepare_ytest(ytest, names(alpha.ptlasso), study_sizes = xtest.norm$n_by_study)
    family <- object$family
    type.measure <- object$type.measure
    metric.preds <- lapply(pred.by.alpha, function(pred) pred$.metric_predictions)
    overall.metric.preds <- lapply(names(alpha.ptlasso), function(study.name) {
      key <- as.character(alpha.ptlasso[[study.name]])
      metric.preds[[key]]$overall[[study.name]]
    })
    ind.metric.preds <- lapply(names(alpha.ptlasso), function(study.name) {
      key <- as.character(alpha.ptlasso[[study.name]])
      metric.preds[[key]]$ind[[study.name]]
    })
    pre.metric.preds <- lapply(names(alpha.ptlasso), function(study.name) {
      key <- as.character(alpha.ptlasso[[study.name]])
      metric.preds[[key]]$pre[[study.name]]
    })
    names(overall.metric.preds) <- names(ind.metric.preds) <- names(pre.metric.preds) <- names(alpha.ptlasso)
    erroverall <- ptmv_summarize_metric(overall.metric.preds, ytest.norm, family, type.measure, add_r2 = family == "gaussian")
    errind <- ptmv_summarize_metric(ind.metric.preds, ytest.norm, family, type.measure, add_r2 = family == "gaussian")
    errpre <- ptmv_summarize_metric(pre.metric.preds, ytest.norm, family, type.measure, add_r2 = family == "gaussian")
    metrics <- ptmv_build_metric_report(
      overall_preds = overall.metric.preds,
      ind_preds = ind.metric.preds,
      pre_preds = pre.metric.preds,
      y = ytest.norm,
      family = family,
      type.measure = type.measure
    )
  }
  
  ptmv_build_prediction_object(
    call = this.call,
    alpha.ptlasso = alpha.ptlasso,
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
    metrics = metrics,
    erroverall = erroverall,
    errind = errind,
    errpre = errpre,
    fit = object,
    metric_predictions = if (!is.null(ytest)) {
      list(overall = overall.metric.preds, ind = ind.metric.preds, pre = pre.metric.preds)
    } else {
      NULL
    },
    class_name = "predict.cv.gptLasso"
  )
}
