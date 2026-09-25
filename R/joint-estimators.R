.update_linear_joint <- function(data, z, y, x, x_star, complete, strata,
                                 x_models, h, initial, intercept, n_quad, nboot,
                                 se_method = "ordinary") {
  call <- match.call()
  .joint_validate_nboot(nboot)
  z <- .joint_z_names(data, z)
  x_models <- .joint_x_models(x_models)
  quadrature <- if (any(x_models %in% c("linear", "gam"))) {
    .joint_quadrature(n_quad)
  } else NULL

  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .regression_joint_once(
    data, family = "linear", z, y, x, x_star, complete, strata,
    x_models, quadrature, tune = NULL, h = h, initial = initial,
    intercept = intercept,
    compute_update = !fully_bootstrap
  )
  fit_once <- function(boot_data) {
    .regression_joint_once(
      boot_data, family = "linear", z, y, x, x_star, complete, strata,
      x_models, quadrature, tune = NULL, h = h, initial = initial,
      intercept = intercept,
      compute_update = !fully_bootstrap
    )[[if (fully_bootstrap) "components" else "estimates"]]
  }
  if (fully_bootstrap) {
    bootstrap_fit <- .fully_bootstrap_fit(
      data, nboot, fit_once, fit$components
    )
    fit$estimates <- bootstrap_fit$estimates
    se <- bootstrap_fit$se
    fit$details$projection <- stats::setNames(
      rep("full_bootstrap", ncol(fit$estimates) - 1L),
      colnames(fit$estimates)[-1L]
    )
    fit$details$full_bootstrap_covariance <- bootstrap_fit$covariance
    fit$details$bootstrap_attempts <- bootstrap_fit$attempts
    fit$details$bootstrap_failures <- bootstrap_fit$failures
  } else {
    se <- .joint_bootstrap_se(data, nboot, fit_once, fit$estimates)
  }
  .new_joint_fit(call, "linear", fit$estimates, se, fit$details)
}

.update_logistic_joint <- function(data, z, y, x, x_star, complete, strata,
                                   x_models, tune, h, initial, intercept,
                                   n_quad, nboot,
                                   se_method = "ordinary") {
  call <- match.call()
  .joint_validate_nboot(nboot)
  z <- .joint_z_names(data, z)
  x_models <- .joint_x_models(x_models)
  quadrature <- if (any(x_models %in% c("linear", "gam"))) {
    .joint_quadrature(n_quad)
  } else NULL
  if (length(x_models) > 0L &&
      (!is.numeric(tune) || length(tune) < 1L || anyNA(tune) ||
       any(!is.finite(tune)) || any(tune < 0))) {
    stop("`tune` must contain non-negative finite numbers.", call. = FALSE)
  }

  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .regression_joint_once(
    data, family = "logistic", z, y, x, x_star, complete, strata,
    x_models, quadrature, tune, h, initial, intercept,
    compute_update = !fully_bootstrap
  )
  fit_once <- function(boot_data) {
    .regression_joint_once(
      boot_data, family = "logistic", z, y, x, x_star, complete, strata,
      x_models, quadrature, tune, h, initial, intercept,
      compute_update = !fully_bootstrap
    )[[if (fully_bootstrap) "components" else "estimates"]]
  }
  if (fully_bootstrap) {
    bootstrap_fit <- .fully_bootstrap_fit(
      data, nboot, fit_once, fit$components
    )
    fit$estimates <- bootstrap_fit$estimates
    se <- bootstrap_fit$se
    fit$details$projection <- stats::setNames(
      rep("full_bootstrap", ncol(fit$estimates) - 1L),
      colnames(fit$estimates)[-1L]
    )
    fit$details$full_bootstrap_covariance <- bootstrap_fit$covariance
    fit$details$bootstrap_attempts <- bootstrap_fit$attempts
    fit$details$bootstrap_failures <- bootstrap_fit$failures
  } else {
    se <- .joint_bootstrap_se(data, nboot, fit_once, fit$estimates)
  }
  .new_joint_fit(call, "logistic", fit$estimates, se, fit$details)
}

.joint_x_models <- function(x_models) {
  if (is.null(x_models) || length(x_models) == 0L) return(character())
  x_models <- unique(match.arg(
    x_models, c("kernel", "linear", "gam"), several.ok = TRUE
  ))
  x_models
}

.regression_joint_once <- function(data, family, z, y, x, x_star, complete,
                                   strata, x_models, quadrature, tune, h,
                                   initial, intercept,
                                   compute_update = TRUE) {
  auxiliary <- .joint_auxiliary_names(x_star, z)
  .check_required_columns(data, c(y, x, auxiliary, complete, strata, z))
  complete_index <- data[[complete]] == 1
  data_compl <- data[complete_index, , drop = FALSE]
  terms <- c(if (intercept) "(Intercept)", x, z)
  p <- length(terms)
  q <- length(auxiliary) + as.integer(intercept)
  if (nrow(data_compl) <= max(p, q)) {
    stop("The Phase II sample is too small for the joint regression model.", call. = FALSE)
  }
  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[complete_index]
  centered <- data[[complete]] * ipw_all - 1
  original_formula <- .joint_formula(y, c(x, z), intercept = intercept)
  working_formula <- .joint_formula(y, auxiliary, intercept = intercept)

  if (family == "linear") {
    original_fit <- stats::lm(original_formula, data = data_compl, weights = ipw_compl)
    working_fit <- stats::lm(working_formula, data = data)
    working_fit_compl <- stats::lm(working_formula, data = data_compl,
                                   weights = ipw_compl)
  } else {
    original_fit <- suppressWarnings(stats::glm(
      original_formula, data = data_compl, family = stats::binomial(),
      weights = ipw_compl
    ))
    working_fit <- stats::glm(working_formula, data = data,
                              family = stats::binomial())
    working_fit_compl <- suppressWarnings(stats::glm(
      working_formula, data = data_compl, family = stats::binomial(),
      weights = ipw_compl
    ))
  }
  est_compl <- as.numeric(stats::coef(original_fit))
  est_work <- as.numeric(stats::coef(working_fit))
  est_work_compl <- as.numeric(stats::coef(working_fit_compl))
  if (any(!is.finite(c(est_compl, est_work, est_work_compl)))) {
    stop("A joint outcome model is rank deficient or did not converge.", call. = FALSE)
  }
  solver_initial <- if (is.null(initial)) est_compl else initial

  kernel_scales <- NULL
  if ("kernel" %in% x_models) {
    standardized_kernel <- length(z) > 0L
    h <- .kernel_bandwidth(
      data_compl, auxiliary, h, standardized = standardized_kernel
    )
    if (standardized_kernel) {
      kernel_scales <- .kernel_scales(data_compl, auxiliary)
    }
  }

  sdy_compl <- if (family == "linear") {
    .weighted_resid_sd(stats::residuals(original_fit), ipw_compl)
  } else {
    NULL
  }
  differences <- list(default = est_work_compl - est_work)
  methods <- list(complete_case = est_compl)
  projection_methods <- character()

  if (compute_update) {
    c_compl <- .joint_design(data_compl, c(x, z), intercept)
    w_all <- .joint_design(data, auxiliary, intercept)
    if (family == "linear") {
      residual_compl <- data_compl[[y]] - as.vector(c_compl %*% est_compl)
      residual_work <- data[[y]] - as.vector(w_all %*% est_work)
      info_compl <- crossprod(c_compl, ipw_compl * c_compl)
      info_work <- crossprod(w_all)
      raw_compl_part <- c_compl * (ipw_compl * residual_compl)
      raw_work <- w_all * (centered * residual_work)
    } else {
      prob_compl <- stats::plogis(as.vector(c_compl %*% est_compl))
      prob_work <- stats::plogis(as.vector(w_all %*% est_work))
      raw_compl_part <- c_compl * (ipw_compl * (data_compl[[y]] - prob_compl))
      raw_work <- w_all * (centered * (data[[y]] - prob_work))
      info_compl <- crossprod(
        c_compl,
        (ipw_compl * prob_compl * (1 - prob_compl)) * c_compl
      )
      info_work <- crossprod(
        w_all, (prob_work * (1 - prob_work)) * w_all
      )
    }
    raw_compl <- matrix(0, nrow(data), p)
    raw_compl[complete_index, ] <- raw_compl_part
    influence_compl <- .joint_influence(raw_compl, info_compl)
    influence_work <- .joint_influence(raw_work, info_work)
    outcome_projection <- .joint_projection(influence_work, influence_compl)
    methods$default <- as.numeric(
      est_compl - outcome_projection$coef %*% differences$default
    )
    projection_methods["default"] <- outcome_projection$method
  }

  add_conditional_update <- function(label, score_full, score_compl, score_vec) {
    est_full <- .joint_solve_score(solver_initial, score_full)
    est_phase2 <- .joint_solve_score(solver_initial, score_compl)
    difference <- est_phase2 - est_full
    joint_label <- paste0("joint_", label)
    differences[[label]] <<- difference
    differences[[joint_label]] <<- c(differences$default, difference)
    if (!compute_update) return(invisible(NULL))

    info <- -numDeriv::jacobian(score_vec$sum, est_full)
    vector <- score_vec$vector(est_full)
    influence <- .joint_influence(vector * centered, info)
    single_projection <- .joint_projection(influence, influence_compl)
    joint_influence <- cbind(influence_work, influence)
    joint_projection <- .joint_projection(joint_influence, influence_compl)
    methods[[label]] <<- as.numeric(
      est_compl - single_projection$coef %*% difference
    )
    methods[[joint_label]] <<- as.numeric(
      est_compl - joint_projection$coef %*%
        differences[[joint_label]]
    )
    projection_methods[label] <<- single_projection$method
    projection_methods[joint_label] <<- joint_projection$method
  }

  for (model_type in x_models) {
    model_tune <- if (length(tune) <= 1L) tune else unname(tune[model_type])
    if (model_type == "kernel") {
      score_full <- function(beta) .regression_kernel_score_impl(
        beta, data, data_compl, auxiliary, h, ipw_compl,
        family = family, sdy = sdy_compl, tune = model_tune,
        y = y, x = x, z = z, intercept = intercept,
        kernel_scales = kernel_scales
      )
      score_compl <- function(beta) .regression_kernel_score_impl(
        beta, data_compl, data_compl, auxiliary, h, ipw_compl,
        ipw = ipw_compl, family = family, sdy = sdy_compl,
        tune = model_tune, y = y, x = x, z = z,
        intercept = intercept, kernel_scales = kernel_scales
      )
      score_vec <- list(
        sum = function(beta) .regression_kernel_score_impl(
          beta, data, data_compl, auxiliary, h, ipw_compl,
          family = family, sdy = sdy_compl,
          y = y, x = x, z = z, intercept = intercept,
          kernel_scales = kernel_scales
        ),
        vector = function(beta) .regression_kernel_score_impl(
          beta, data, data_compl, auxiliary, h, ipw_compl,
          family = family, sdy = sdy_compl,
          y = y, x = x, z = z, intercept = intercept, vector = TRUE,
          kernel_scales = kernel_scales
        )
      )
    } else if (family == "linear") {
      x_fit <- .joint_fit_x_model(
        data_compl, x, auxiliary, ipw_compl, model_type
      )
      score_full <- function(beta) linear_score_modelx(
        beta, data, x_fit$fit, x_fit$se, sdy_compl,
        quadrature$nodes, quadrature$weights, z = z,
        intercept = intercept, y = y
      )
      score_compl <- function(beta) linear_score_modelx(
        beta, data_compl, x_fit$fit, x_fit$se, sdy_compl,
        quadrature$nodes, quadrature$weights, ipw_compl, z,
        intercept, y
      )
      score_vec <- list(
        sum = score_full,
        vector = function(beta) .linear_score_modelx_impl(
          beta, data, x_fit$fit, x_fit$se, sdy_compl,
          quadrature$nodes, quadrature$weights, NULL, vector = TRUE, z = z,
          intercept = intercept, y = y
        )
      )
    } else {
      x_fit <- .joint_fit_x_model(
        data_compl, x, auxiliary, ipw_compl, model_type
      )
      score_full <- function(beta) logistic_score_modelx(
        beta, data, x_fit$fit, x_fit$se,
        quadrature$nodes, quadrature$weights, tune = model_tune, z = z,
        intercept = intercept, y = y
      )
      score_compl <- function(beta) logistic_score_modelx(
        beta, data_compl, x_fit$fit, x_fit$se,
        quadrature$nodes, quadrature$weights, ipw_compl, model_tune, z,
        intercept, y
      )
      score_vec <- list(
        sum = function(beta) logistic_score_modelx(
          beta, data, x_fit$fit, x_fit$se,
          quadrature$nodes, quadrature$weights, z = z,
          intercept = intercept, y = y
        ),
        vector = function(beta) .logistic_score_modelx_impl(
          beta, data, x_fit$fit, x_fit$se,
          quadrature$nodes, quadrature$weights, NULL, z, intercept, y
        )
      )
    }
    add_conditional_update(model_type, score_full, score_compl, score_vec)
  }

  components <- .bootstrap_components(
    stats::setNames(est_compl, terms), differences
  )
  estimates <- do.call(cbind, methods)
  rownames(estimates) <- terms
  list(
    estimates = estimates,
    components = components,
    details = list(
      z = z, ipw = ipw_all, projection = projection_methods,
      x_models = x_models, tune = tune, h = h, initial = initial,
      intercept = intercept,
      sdy = sdy_compl,
      quadrature_points = if (is.null(quadrature)) 0L else length(quadrature$nodes)
    )
  )
}
