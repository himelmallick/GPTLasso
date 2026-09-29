cvar.multiview <- function(x.list, y,
                           family = gaussian(),
                           alpha = c(0, seq(0.05, 0.95, 0.05), 1),
                           rho = c(0, 0.1, 0.25, 0.5, 1, 5, 10),
                           lambda = NULL,
                           s = c("lambda.1se", "lambda.min"),
                           nfolds = 10,
                           foldid = NULL,
                           weights = NULL,
                           offset = NULL,
                           cox.groups = NULL,
                           penalty.factor = NULL,
                           type.measure = NULL,
                           alignment = c("lambda", "fraction"),
                           keep = FALSE,
                           trace.it = 0,
                           seed = NULL,
                           verbose = FALSE,
                           ...) {
  mv_fit <- function(x.list, y, family, alpha, rho, lambda = NULL,
                     weights = NULL, offset = NULL,
                     penalty.factor = NULL, trace.it = 0, ...) {
    args <- list(
      x_list = x.list, y = y, family = family, alpha = alpha, rho = rho,
      weights = weights, offset = offset, lambda = lambda,
      trace.it = trace.it, ...
    )
    if (!is.null(penalty.factor)) args$penalty.factor <- penalty.factor
    fit_once <- function(extra = list()) do.call(multiview::multiview, c(args, extra))
    out <- tryCatch(fit_once(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    if (grepl("cannot correct step size", conditionMessage(out), fixed = TRUE)) {
      out2 <- tryCatch(
        fit_once(list(lambda.min.ratio = 0.05, nlambda = 50, maxit = 1e5)),
        error = function(e) e
      )
      if (!inherits(out2, "error")) return(out2)
    }
    stop(conditionMessage(out), call. = FALSE)
  }
  
  make_stratified_foldid <- function(y01, K) {
    y01 <- as.integer(y01)
    if (!all(y01 %in% c(0L, 1L))) stop("Stratification indicator must be 0/1.")
    fid <- integer(length(y01))
    for (value in c(0L, 1L)) {
      idx <- which(y01 == value)
      fid[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
    }
    fid
  }
  
  validate_binomial_folds <- function(y01, fid) {
    for (f in sort(unique(fid))) {
      if (length(unique(y01[fid != f])) < 2L) {
        stop(sprintf(
          "Invalid CV split for binomial: training set in fold %d has only one class. Use stratified foldid or reduce nfolds.",
          f
        ))
      }
    }
    invisible(TRUE)
  }
  
  build_predmat_mv <- function(outlist, lambda, x.list, foldid,
                               fold.offsets = NULL,
                               alignment = c("lambda", "fraction"),
                               family, type = "response", ...) {
    alignment <- match.arg(alignment)
    N <- nrow(x.list[[1L]])
    predmat <- matrix(NA_real_, N, length(lambda))
    nlambda <- length(lambda)
    for (f in sort(unique(foldid))) {
      test <- foldid == f
      fitobj <- outlist[[f]]
      x_sub <- lapply(x.list, function(x) x[test, , drop = FALSE])
      offset_sub <- if (is.null(fold.offsets)) NULL else fold.offsets[[f]][test]
      preds <- switch(
        alignment,
        fraction = predict(fitobj, newx = x_sub, newoffset = offset_sub,
                           type = type, ...),
        lambda = predict(fitobj, newx = x_sub, s = lambda,
                         newoffset = offset_sub, type = type, ...)
      )
      preds <- as.matrix(preds)
      nlami <- min(ncol(preds), nlambda)
      predmat[test, seq_len(nlami)] <- preds[, seq_len(nlami), drop = FALSE]
      if (nlami < nlambda) predmat[test, (nlami + 1L):nlambda] <- preds[, nlami]
    }
    dimnames(predmat) <- list(rownames(x.list[[1L]]), paste0("s", seq_len(nlambda)))
    attr(predmat, "family") <- family
    predmat
  }
  
  score_fold_metric <- function(y_test, preds_fold, family_name, type.measure,
                                weights_test = NULL) {
    if (family_name == "cox") {
      if (type.measure != "C") stop("Internal Cox scoring expects type.measure = 'C'.")
      return(apply(preds_fold, 2, function(p) glmnet::Cindex(p, y_test, weights = weights_test)))
    }
    if (type.measure == "mse") return(colMeans((as.numeric(y_test) - preds_fold)^2))
    if (type.measure == "deviance") {
      if (family_name == "gaussian") return(colMeans((as.numeric(y_test) - preds_fold)^2))
      if (family_name == "binomial") {
        eps <- 1e-8
        yv <- as.numeric(y_test)
        return(apply(preds_fold, 2, function(p) {
          pp <- pmin(pmax(as.numeric(p), eps), 1 - eps)
          -2 * mean(yv * log(pp) + (1 - yv) * log(1 - pp))
        }))
      }
    }
    if (type.measure == "auc") {
      return(apply(preds_fold, 2, function(p) {
        if (length(unique(stats::na.omit(y_test))) < 2L) return(NA_real_)
        as.numeric(pROC::auc(pROC::roc(
          response = y_test, predictor = p, levels = c(0, 1),
          direction = "<", quiet = TRUE
        )))
      }))
    }
    stop("Unsupported type.measure in cvar.multiview().")
  }
  
  choose_best_index <- function(vals, metric) {
    if (metric %in% c("auc", "C")) which.max(vals) else which.min(vals)
  }
  choose_1se <- function(lambda, cvmean, cvse, idx_best, metric) {
    eligible <- if (metric %in% c("auc", "C")) {
      cvmean >= cvmean[idx_best] - cvse[idx_best]
    } else {
      cvmean <= cvmean[idx_best] + cvse[idx_best]
    }
    max(lambda[eligible], na.rm = TRUE)
  }
  weighted_summary <- function(cvm, fold_weights) {
    means <- apply(cvm, 1, function(z) {
      ok <- is.finite(z)
      if (!any(ok)) return(NA_real_)
      stats::weighted.mean(z[ok], fold_weights[ok])
    })
    ses <- vapply(seq_len(nrow(cvm)), function(j) {
      z <- cvm[j, ]
      ok <- is.finite(z)
      if (sum(ok) <= 1L) return(NA_real_)
      sqrt(stats::weighted.mean((z[ok] - means[j])^2, fold_weights[ok]) /
             (sum(ok) - 1L))
    }, numeric(1))
    list(mean = means, se = ses)
  }
  
  alignment <- match.arg(alignment)
  if (!is.null(seed)) set.seed(seed)
  family_name <- ptmv_family_name(family)
  family_arg <- if (family_name == "cox") "cox" else family
  N <- nrow(x.list[[1L]])
  if (any(vapply(x.list, nrow, integer(1)) != N)) {
    stop("All views in x.list must have the same number of rows.")
  }
  if (ptmv_response_nobs(y) != N) {
    stop("Number of observations in y must equal nrow(x.list[[1]]).")
  }
  if (family_name == "cox") ptmv_validate_surv(y, context = "y") else y <- drop(y)
  if (is.null(weights)) weights <- rep(1, N)
  if (length(weights) != N || any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must contain one finite non-negative value per observation.")
  }
  if (!is.null(offset) && length(offset) != N) stop("Length of offset must equal N when provided.")
  base_offset <- if (is.null(offset)) rep(0, N) else as.numeric(offset)
  if (!is.null(cox.groups)) {
    if (family_name != "cox") stop("cox.groups is only available for family = 'cox'.")
    if (length(cox.groups) != N) stop("cox.groups must have length N.")
    cox.groups <- factor(cox.groups)
  }
  
  if (is.null(type.measure)) {
    type.measure <- if (family_name == "gaussian") "mse" else if (family_name == "binomial") "auc" else "deviance"
  }
  allowed <- switch(family_name,
                    gaussian = c("mse", "deviance"), binomial = c("auc", "deviance"),
                    cox = c("C", "deviance"), character()
  )
  if (!type.measure %in% allowed) {
    stop(sprintf("For %s family, type.measure must be one of: %s.", family_name, paste(allowed, collapse = ", ")))
  }
  if (length(s) != 1L) s <- s[1L]
  if (is.character(s) && identical(s, "lambda.lse")) s <- "lambda.1se"
  s <- match.arg(s, c("lambda.1se", "lambda.min"))
  alpha <- sort(unique(alpha))
  rho <- sort(unique(rho))
  
  if (is.null(foldid)) {
    if (family_name == "binomial") {
      foldid <- make_stratified_foldid(as.integer(as.numeric(y) > 0), nfolds)
    } else if (family_name == "cox") {
      strata <- if (is.null(cox.groups)) factor(rep("all", N)) else cox.groups
      event_counts <- tapply(as.integer(y[, 2]), strata, sum)
      if (any(event_counts < 3L)) {
        stop("Cox cross-validation requires at least 3 observed events in every study.")
      }
      nfolds <- min(nfolds, min(event_counts), min(table(strata)))
      foldid <- integer(N)
      for (group_name in levels(strata)) {
        idx <- which(strata == group_name)
        foldid[idx] <- make_stratified_foldid(as.integer(y[idx, 2]), nfolds)
      }
    } else {
      foldid <- sample(rep(seq_len(nfolds), length.out = N))
    }
  } else {
    if (length(foldid) != N) stop("foldid must have length N (nrow of data).")
    foldid <- ptmv_renumber_foldid(foldid)
  }
  K <- length(unique(foldid))
  if (K < 2L) stop("foldid must contain at least 2 unique folds.")
  if (family_name == "binomial") validate_binomial_folds(as.integer(as.numeric(y) > 0), foldid)
  if (family_name == "cox") {
    validation_groups <- if (is.null(cox.groups)) factor(rep("all", N)) else cox.groups
    ptmv_validate_cox_folds(y, validation_groups, foldid)
  }
  fold_levels <- sort(unique(foldid))
  if (any(tapply(weights, foldid, sum) <= 0)) {
    stop("Every validation fold must have positive total weight.")
  }
  
  full_shift <- list(offset = rep(0, N), shift = NULL, fit = NULL)
  if (family_name == "cox" && !is.null(cox.groups)) {
    full_shift <- ptmv_fit_cox_study_shift(y, cox.groups, weights = weights)
  }
  full_offset <- if (family_name == "cox") base_offset + full_shift$offset else offset
  fold_offsets <- NULL
  if (family_name == "cox") {
    fold_offsets <- vector("list", K)
    for (f in fold_levels) {
      shift_f <- if (is.null(cox.groups)) list(offset = rep(0, N), shift = NULL) else {
        ptmv_fit_cox_study_shift(y, cox.groups, train = foldid != f, weights = weights)
      }
      fold_offsets[[f]] <- base_offset + shift_f$offset
      attr(fold_offsets[[f]], "study.shift") <- shift_f$shift
    }
  } else if (!is.null(offset)) {
    fold_offsets <- rep(list(as.numeric(offset)), K)
  }
  
  grid <- expand.grid(alpha = alpha, rho = rho)
  combo_results <- vector("list", nrow(grid))
  score_choice <- lambda_choice_vec <- lambda_min_vec <- lambda_1se_vec <- rep(NA_real_, nrow(grid))
  x_all <- do.call(cbind, x.list)
  
  for (g in seq_len(nrow(grid))) {
    a <- grid$alpha[g]
    r <- grid$rho[g]
    if (trace.it || verbose) cat(sprintf("Tuning alpha = %.3f, rho = %.3f (%d/%d)\n", a, r, g, nrow(grid)))
    combo_attempt <- tryCatch({
      if (is.null(lambda)) {
        master <- mv_fit(x.list, y, family_arg, a, r, weights = weights,
                         offset = full_offset, penalty.factor = penalty.factor, ...)
        lambda_g <- master$lambda
      } else lambda_g <- lambda
      
      outlist <- vector("list", K)
      for (f in fold_levels) {
        test <- foldid == f
        outlist[[f]] <- mv_fit(
          lapply(x.list, function(m) m[!test, , drop = FALSE]),
          ptmv_subset_response(y, !test), family_arg, a, r,
          lambda = lambda_g, weights = weights[!test],
          offset = if (is.null(fold_offsets)) NULL else fold_offsets[[f]][!test],
          penalty.factor = penalty.factor, ...
        )
      }
      pred_type <- if (family_name == "cox") "link" else "response"
      predmat <- build_predmat_mv(outlist, lambda_g, x.list, foldid,
                                  fold.offsets = fold_offsets, alignment = alignment,
                                  family = family_arg, type = pred_type)
      cvm <- matrix(NA_real_, nrow = length(lambda_g), ncol = K)
      fold_weights <- numeric(K)
      for (i in seq_along(fold_levels)) {
        f <- fold_levels[i]
        test <- foldid == f
        fold_weights[i] <- sum(weights[test])
        if (family_name == "cox" && type.measure == "deviance") {
          coefmat <- as.matrix(if (alignment == "fraction") {
            stats::coef(outlist[[f]])
          } else {
            stats::coef(outlist[[f]], s = lambda_g)
          })
          if (nrow(coefmat) != ncol(x_all)) stop("Cox coefficient path does not match the concatenated feature count.")
          nlami <- min(ncol(coefmat), length(lambda_g))
          dev_full <- glmnet::coxnet.deviance(
            x = x_all, y = y, offset = fold_offsets[[f]], weights = weights,
            beta = coefmat[, seq_len(nlami), drop = FALSE],
            std.weights = FALSE, cox.ties = "breslow"
          )
          dev_train <- glmnet::coxnet.deviance(
            x = x_all[!test, , drop = FALSE], y = ptmv_subset_response(y, !test),
            offset = fold_offsets[[f]][!test], weights = weights[!test],
            beta = coefmat[, seq_len(nlami), drop = FALSE],
            std.weights = FALSE, cox.ties = "breslow"
          )
          values <- (dev_full - dev_train) / fold_weights[i]
          cvm[seq_len(nlami), i] <- values
          if (nlami < length(lambda_g)) cvm[(nlami + 1L):length(lambda_g), i] <- values[nlami]
        } else {
          cvm[, i] <- score_fold_metric(
            ptmv_subset_response(y, test), predmat[test, , drop = FALSE],
            family_name, type.measure, weights[test]
          )
        }
      }
      summary <- weighted_summary(cvm, fold_weights)
      idx_best <- choose_best_index(summary$mean, type.measure)
      lambda_min <- lambda_g[idx_best]
      lambda_1se <- choose_1se(lambda_g, summary$mean, summary$se, idx_best, type.measure)
      lambda_choice <- if (s == "lambda.min") lambda_min else lambda_1se
      choice_idx <- which.min(abs(lambda_g - lambda_choice))
      list(
        score = summary$mean[choice_idx], lambda.choice = lambda_choice,
        lambda.min = lambda_min, lambda.1se = lambda_1se,
        combo.result = list(
          alpha = a, rho = r, lambda = lambda_g, cvm = cvm,
          cvmean = summary$mean, cvse = summary$se,
          lambda.min = lambda_min, lambda.1se = lambda_1se,
          lambda.choice = lambda_choice,
          fit.preval = if (keep) predmat else NULL,
          fit.preval.offset = if (keep && family_name == "cox") fold_offsets else NULL,
          fit.preval.study.shift = if (keep && family_name == "cox") {
            lapply(fold_offsets, attr, which = "study.shift")
          } else NULL,
          foldid = foldid, type.measure = type.measure, s = s, alignment = alignment
        )
      )
    }, error = function(e) e)
    if (inherits(combo_attempt, "error")) {
      if (trace.it || verbose) message(sprintf(
        "Skipping alpha = %.3f, rho = %.3f because fitting failed: %s",
        a, r, conditionMessage(combo_attempt)
      ))
      next
    }
    score_choice[g] <- combo_attempt$score
    lambda_choice_vec[g] <- combo_attempt$lambda.choice
    lambda_min_vec[g] <- combo_attempt$lambda.min
    lambda_1se_vec[g] <- combo_attempt$lambda.1se
    combo_results[[g]] <- combo_attempt$combo.result
  }
  
  valid_idx <- which(is.finite(score_choice) & !vapply(combo_results, is.null, logical(1)))
  if (length(valid_idx) == 0L) stop("All alpha/rho combinations failed in cvar.multiview().")
  best_idx <- valid_idx[choose_best_index(score_choice[valid_idx], type.measure)]
  alpha_choice <- grid$alpha[best_idx]
  rho_choice <- grid$rho[best_idx]
  lambda_choice <- lambda_choice_vec[best_idx]
  final_fit <- mv_fit(x.list, y, family_arg, alpha_choice, rho_choice,
                      lambda = lambda_choice, weights = weights,
                      offset = full_offset, penalty.factor = penalty.factor, ...)
  selected <- combo_results[[best_idx]]
  idx_min <- which.min(abs(selected$lambda - selected$lambda.min))
  idx_1se <- which.min(abs(selected$lambda - selected$lambda.1se))
  idx_alpha_choice <- valid_idx[grid$alpha[valid_idx] == alpha_choice]
  idx_alpha_choice <- idx_alpha_choice[order(grid$rho[idx_alpha_choice])]
  cv_by_rho <- lapply(idx_alpha_choice, function(i) {
    ri <- combo_results[[i]]
    list(rho = ri$rho, lambda = ri$lambda, cvm = ri$cvm,
         cvmean = ri$cvmean, cvse = ri$cvse, lambda.min = ri$lambda.min,
         lambda.1se = ri$lambda.1se, lambda.choice = ri$lambda.choice,
         type.measure = ri$type.measure, s = ri$s, alignment = ri$alignment,
         keep = keep, fit.preval = ri$fit.preval,
         fit.preval.offset = ri$fit.preval.offset,
         fit.preval.study.shift = ri$fit.preval.study.shift,
         foldid = ri$foldid)
  })
  rho_seq <- vapply(cv_by_rho, `[[`, numeric(1), "rho")
  rho_lambda_min <- vapply(cv_by_rho, `[[`, numeric(1), "lambda.min")
  rho_lambda_1se <- vapply(cv_by_rho, `[[`, numeric(1), "lambda.1se")
  rho_lambda_choice <- vapply(cv_by_rho, `[[`, numeric(1), "lambda.choice")
  rho_cvmean <- vapply(cv_by_rho, function(z) z$cvmean[which.min(abs(z$lambda - z$lambda.choice))], numeric(1))
  rho_cvse <- vapply(cv_by_rho, function(z) z$cvse[which.min(abs(z$lambda - z$lambda.choice))], numeric(1))
  
  out <- list(
    call = match.call(), multiview.fit = final_fit, lambda = selected$lambda,
    cvm = selected$cvmean, cvsd = selected$cvse,
    cvup = selected$cvmean + selected$cvse, cvlo = selected$cvmean - selected$cvse,
    lambda.min = selected$lambda.min, lambda.1se = selected$lambda.1se,
    index = matrix(c(idx_min, idx_1se), ncol = 1,
                   dimnames = list(c("min", "1se"), "Lambda")),
    fit.preval = if (keep) selected$fit.preval else NULL,
    fit.preval.offset = if (keep) selected$fit.preval.offset else NULL,
    fit.preval.study.shift = if (keep) selected$fit.preval.study.shift else NULL,
    foldid = if (keep) foldid else NULL, type.measure = type.measure, s = s,
    alpha = alpha, rho = rho_seq, cv_by_rho = cv_by_rho,
    rho.cvmean = rho_cvmean, rho.cvse = rho_cvse,
    rho.lambda.min = rho_lambda_min, rho.lambda.1se = rho_lambda_1se,
    rho.lambda.choice = rho_lambda_choice, alpha.choice = alpha_choice,
    rho.choice = rho_choice, lambda.choice = lambda_choice,
    alpha.rho.grid = grid, combo.results = combo_results,
    score.choice = score_choice, alpha.rho.choice.score = score_choice[best_idx],
    alpha.rho.lambda.min = lambda_min_vec, alpha.rho.lambda.1se = lambda_1se_vec,
    alpha.rho.lambda.choice = lambda_choice_vec,
    keep = keep, alignment = alignment, nfolds = K, family = family_arg,
    cox.study.effect = if (family_name == "cox" && !is.null(cox.groups)) "fixed_shift" else NULL,
    cox.group.shift = full_shift$shift, offset.full = full_offset
  )
  class(out) <- c("cvar.multiview", "cv.multiview", "cv.glmnet")
  out
}