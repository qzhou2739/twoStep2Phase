.update_cox_joint <- function(data, z, time, event, x, x_star, complete,
                              strata, x_models, h, initial, n_quad, nboot,
                              se_method = "ordinary") {
  call <- match.call()
  .joint_validate_nboot(nboot)
  z <- .joint_z_names(data, z)
  x_models <- .joint_x_models(x_models)
  quadrature <- if (any(x_models %in% c("linear", "gam"))) {
    .joint_quadrature(n_quad)
  } else NULL

  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .cox_joint_once(
    data, z, time, event, x, x_star, complete, strata, x_models,
    quadrature, h, initial, compute_update = !fully_bootstrap
  )
  fit_once <- function(boot_data) {
    .cox_joint_once(
      boot_data, z, time, event, x, x_star, complete, strata, x_models,
      quadrature, h, initial, compute_update = !fully_bootstrap
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
  .new_joint_fit(call, "cox", fit$estimates, se, fit$details)
}

.cox_joint_once <- function(data, z, time, event, x, x_star, complete,
                            strata, x_models, quadrature, h, initial,
                            compute_update = TRUE) {
  auxiliary_names <- .joint_auxiliary_names(x_star, z)
  .check_required_columns(data, c(
    time, event, x, auxiliary_names, complete, strata, z
  ))
  complete_index <- data[[complete]] == 1
  data_compl <- data[complete_index, , drop = FALSE]
  p <- length(z) + 1L
  if (nrow(data_compl) <= max(p, length(auxiliary_names)) ||
      sum(data_compl[[event]] == 1) == 0) {
    stop("The Phase II sample is too small for the joint Cox model.", call. = FALSE)
  }
  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[complete_index]
  centered <- data[[complete]] * ipw_all - 1
  terms <- c(x, z)
  response <- paste0("survival::Surv(", time, ", ", event, ")")
  original_formula <- .joint_formula(response, terms)
  working_formula <- .joint_formula(response, auxiliary_names)

  original_fit <- survival::coxph(original_formula, data = data_compl,
                                  weights = ipw_compl)
  working_fit <- survival::coxph(working_formula, data = data)
  working_fit_compl <- survival::coxph(working_formula, data = data_compl,
                                       weights = ipw_compl)
  est_compl <- as.numeric(stats::coef(original_fit))
  est_work <- as.numeric(stats::coef(working_fit))
  est_work_compl <- as.numeric(stats::coef(working_fit_compl))
  if (any(!is.finite(c(est_compl, est_work, est_work_compl)))) {
    stop("A joint Cox model is rank deficient or did not converge.", call. = FALSE)
  }
  solver_initial <- if (is.null(initial)) est_compl else initial

  kernel_scales <- NULL
  if ("kernel" %in% x_models) {
    standardized_kernel <- length(z) > 0L
    h <- .kernel_bandwidth(
      data_compl, auxiliary_names, h, standardized = standardized_kernel
    )
    if (standardized_kernel) {
      kernel_scales <- .kernel_scales(data_compl, auxiliary_names)
    }
  }

  c_compl <- .joint_design(data_compl, c(x, z))
  a_all <- .joint_design(data, auxiliary_names)
  a_compl <- .joint_design(data_compl, auxiliary_names)
  kernel_a_all <- a_all
  kernel_a_compl <- a_compl
  if (!is.null(kernel_scales)) {
    kernel_a_all <- sweep(kernel_a_all, 2L, kernel_scales, FUN = "/")
    kernel_a_compl <- sweep(kernel_a_compl, 2L, kernel_scales, FUN = "/")
  }
  differences <- list(default = est_work_compl - est_work)
  methods <- list(complete_case = est_compl)
  projection_methods <- character()

  if (compute_update) {
    info_compl <- .cox_general_info(
      est_compl, data, data_compl, c_compl, ipw_compl, time, event
    )
    info_work <- .cox_general_info(
      est_work, data, data, a_all, rep(1, nrow(data)), time, event
    )
    score_compl <- .cox_general_score_vec(
      est_compl, data_compl, data, c_compl, data_compl, c_compl,
      ipw_compl, rep(1, nrow(data_compl)), time, event
    )
    raw_compl <- matrix(0, nrow(data), p)
    raw_compl[complete_index, ] <- score_compl * ipw_compl
    score_work <- .cox_general_score_vec(
      est_work, data, data, a_all, data, a_all,
      rep(1, nrow(data)), rep(1, nrow(data)), time, event
    )
    influence_compl <- .joint_influence(raw_compl, info_compl)
    influence_work <- .joint_influence(score_work * centered, info_work)
    outcome_projection <- .joint_projection(influence_work, influence_compl)
    methods$default <- as.numeric(
      est_compl - outcome_projection$coef %*% differences$default
    )
    projection_methods["default"] <- outcome_projection$method
  }

  add_conditional_update <- function(label, vector_full, vector_compl) {
    score_full <- function(beta) colSums(vector_full(beta))
    score_phase2 <- function(beta) colSums(vector_compl(beta))
    est_full <- .joint_solve_score(solver_initial, score_full)
    est_phase2 <- .joint_solve_score(solver_initial, score_phase2)
    difference <- est_phase2 - est_full
    joint_label <- paste0("joint_", label)
    differences[[label]] <<- difference
    differences[[joint_label]] <<- c(differences$default, difference)
    if (!compute_update) return(invisible(NULL))

    information <- -numDeriv::jacobian(score_full, est_full)
    influence <- .joint_influence(vector_full(est_full) * centered, information)
    single_projection <- .joint_projection(influence, influence_compl)
    joint_projection <- .joint_projection(
      cbind(influence_work, influence), influence_compl
    )
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
    if (model_type == "kernel") {
      vector_full <- function(beta) .cox_general_kernel_score_vec(
        beta, data, data, data_compl, z, c_compl, kernel_a_all, kernel_a_compl,
        h, est_compl, ipw_compl, rep(1, nrow(data)), time, event
      )
      vector_compl <- function(beta) .cox_general_kernel_score_vec(
        beta, data_compl, data, data_compl, z, c_compl,
        kernel_a_compl, kernel_a_compl, h, est_compl, ipw_compl, ipw_compl,
        time, event
      )
    } else {
      x_fit <- .joint_fit_x_model(
        data_compl, x, auxiliary_names, ipw_compl, model_type
      )
      pred_all <- as.numeric(stats::predict(x_fit$fit, data))
      pred_compl <- as.numeric(stats::predict(x_fit$fit, data_compl))
      if (!is.null(x_star) && length(z) > 0L) {
        vector_full <- function(beta) .cox_joint_conditional_score_vec(
          beta, data, data, data_compl, c_compl, ipw_compl,
          rep(1, nrow(data)), quadrature$nodes, quadrature$weights,
          time, event, z, pred_all, x_fit$se, est_compl
        )
        vector_compl <- function(beta) .cox_joint_conditional_score_vec(
          beta, data_compl, data, data_compl, c_compl, ipw_compl,
          ipw_compl, quadrature$nodes, quadrature$weights,
          time, event, z, pred_compl, x_fit$se, est_compl
        )
      } else {
        vector_full <- function(beta) .cox_general_model_score_vec(
          beta, data, data, data_compl, z, c_compl, pred_all, x_fit$se,
          quadrature$nodes, quadrature$weights, est_compl, ipw_compl,
          rep(1, nrow(data)), time, event
        )
        vector_compl <- function(beta) .cox_general_model_score_vec(
          beta, data_compl, data, data_compl, z, c_compl,
          pred_compl, x_fit$se, quadrature$nodes, quadrature$weights,
          est_compl, ipw_compl, ipw_compl, time, event
        )
      }
    }
    add_conditional_update(model_type, vector_full, vector_compl)
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
      x_models = x_models, h = h, initial = initial,
      quadrature_points = if (is.null(quadrature)) 0L else length(quadrature$nodes)
    )
  )
}
