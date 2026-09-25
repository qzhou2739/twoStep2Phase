.update_linear_optimal <- function(data, y = "Y", x = "X", x_star = "X_star",
                                  complete = "R", strata = "ind_stra",
                                  beta_interval = NULL, h = NULL, sdy = NULL,
                                  estimate_kernel = TRUE, nboot = 0,
                                  se_method = "ordinary") {
  call <- match.call()
  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .linear_update_once(
    data, y, x, x_star, complete, strata, beta_interval, h, sdy,
    estimate_kernel, compute_update = !fully_bootstrap
  )

  se <- NULL
  if (fully_bootstrap) {
    bootstrap_fit <- .fully_bootstrap_fit(
      data, nboot,
      function(boot_data) {
        .linear_update_once(
          boot_data, y, x, x_star, complete, strata, beta_interval, h, sdy,
          estimate_kernel, compute_update = FALSE
        )$components
      },
      fit$components
    )
    fit$estimates <- stats::setNames(
      as.numeric(bootstrap_fit$estimates[1, ]),
      colnames(bootstrap_fit$estimates)
    )
    se <- stats::setNames(
      as.numeric(bootstrap_fit$se[1, ]), colnames(bootstrap_fit$se)
    )
    fit$details$full_bootstrap_covariance <- bootstrap_fit$covariance
    fit$details$bootstrap_attempts <- bootstrap_fit$attempts
    fit$details$bootstrap_failures <- bootstrap_fit$failures
  } else if (nboot > 0) {
    boot <- replicate(nboot, {
      idx <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
      out <- try(.linear_update_once(data[idx, , drop = FALSE], y, x, x_star, complete, strata,
                                     beta_interval, h, sdy, estimate_kernel), silent = TRUE)
      if (inherits(out, "try-error")) rep(NA_real_, length(fit$estimates)) else out$estimates
    })
    se <- sqrt(apply(boot, 1, .safe_var))
    names(se) <- names(fit$estimates)
  }

  .new_twostep_fit(call, "linear", fit$estimates, se, fit$details)
}

.linear_update_once <- function(data, y, x, x_star, complete, strata,
                                beta_interval, h, sdy, estimate_kernel,
                                compute_update = TRUE) {
  .check_required_columns(data, c(y, x, x_star, complete, strata))
  data_compl <- .complete_data(data, complete)
  if (nrow(data_compl) < 2) {
    stop("At least two complete observations are required.", call. = FALSE)
  }

  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[data[[complete]] == 1]
  f_x <- stats::as.formula(paste(y, "~ 0 +", x))
  f_xstar <- stats::as.formula(paste(y, "~ 0 +", x_star))

  res_compl <- stats::lm(f_x, data = data_compl, weights = ipw_compl)
  est_compl <- unname(stats::coef(res_compl)[1])
  if (is.null(sdy)) {
    sdy <- .weighted_resid_sd(stats::residuals(res_compl), ipw_compl)
  }

  res_star <- stats::lm(f_xstar, data = data)
  res_star_compl <- stats::lm(f_xstar, data = data_compl, weights = ipw_compl)
  est_star <- unname(stats::coef(res_star)[1])
  est_star_compl <- unname(stats::coef(res_star_compl)[1])

  details <- list(ipw = ipw_all, auxiliary_estimates = c(full = est_star, complete = est_star_compl),
                  h = h, sdy = sdy)
  differences <- list(default = est_star_compl - est_star)

  if (estimate_kernel) {
    if (is.null(h)) {
      h <- 1.06 * stats::sd(data_compl[[x_star]]) * (nrow(data_compl)^(-1 / 5))
    }
    if (is.null(beta_interval)) {
      beta_interval <- est_compl + c(-0.3, 0.3)
    }
    est_star_kernel <- stats::uniroot(
      function(beta) .linear_score_kernel(
        beta, data, data_compl, h, sdy, ipw_compl, NULL, y, x, x_star
      ),
      beta_interval, extendInt = "yes"
    )$root
    est_star_compl_kernel <- stats::uniroot(
      function(beta) .linear_score_kernel(
        beta, data_compl, data_compl, h, sdy, ipw_compl, ipw_compl,
        y, x, x_star
      ),
      beta_interval, extendInt = "yes"
    )$root
    details$kernel_auxiliary_estimates <- c(full = est_star_kernel, complete = est_star_compl_kernel)
    details$h <- h
    differences$kernel <- est_star_compl_kernel - est_star_kernel
  }

  components <- .bootstrap_components(est_compl, differences)
  if (!compute_update) {
    return(list(estimates = NULL, details = details, components = components))
  }

  info_compl <- sum(ipw_compl * data_compl[[x]]^2)
  info_linear <- sum(data[[x_star]]^2)
  sigma12_mid <- sum(ipw_compl * (ipw_compl - 1) * data_compl[[x]] *
                       (data_compl[[y]] - est_compl * data_compl[[x]]) *
                       data_compl[[x_star]] *
                       (data_compl[[y]] - est_star * data_compl[[x_star]]))
  sigma22_mid <- sum(((data[[complete]] * ipw_all - 1) * data[[x_star]] *
                        (data[[y]] - est_star * data[[x_star]]))^2)
  sigma12 <- sigma12_mid / (info_compl * info_linear)
  sigma22 <- sigma22_mid / (info_linear^2)
  estimates <- c(
    complete_case = est_compl,
    default = est_compl - sigma12 / sigma22 * differences$default
  )

  if (estimate_kernel) {
    info_kernel <- -.finite_grad(
      function(beta) .linear_score_kernel(
        beta, data, data_compl, h, sdy, ipw_compl, NULL, y, x, x_star
      ),
      est_star_kernel
    )
    sigma12_kernel <- .linear_score2_org_kernel(est_compl, est_star_kernel, data_compl, data_compl,
                                                h, sdy, ipw_compl, ipw_compl,
                                                y, x, x_star, complete) /
      (info_compl * info_kernel)
    sigma22_kernel <- .linear_score2_kernel(est_star_kernel, data, data_compl, h, sdy,
                                            ipw_compl, ipw_all,
                                            y, x, x_star, complete) /
      (info_kernel^2)
    estimates <- c(
      estimates,
      kernel = est_compl - sigma12_kernel / sigma22_kernel * differences$kernel
    )
  }

  list(estimates = estimates, details = details, components = components)
}

.update_logistic_optimal <- function(data, y = "Y", x = "X", x_star = "X_star",
                                    complete = "R", strata = "ind_stra",
                                    beta_interval = NULL, h = NULL, tune = 0.01,
                                    estimate_kernel = TRUE, nboot = 0,
                                    se_method = "ordinary") {
  call <- match.call()
  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .logistic_update_once(
    data, y, x, x_star, complete, strata, beta_interval, h, tune,
    estimate_kernel, compute_update = !fully_bootstrap
  )

  se <- NULL
  if (fully_bootstrap) {
    bootstrap_fit <- .fully_bootstrap_fit(
      data, nboot,
      function(boot_data) {
        .logistic_update_once(
          boot_data, y, x, x_star, complete, strata, beta_interval, h, tune,
          estimate_kernel, compute_update = FALSE
        )$components
      },
      fit$components
    )
    fit$estimates <- stats::setNames(
      as.numeric(bootstrap_fit$estimates[1, ]),
      colnames(bootstrap_fit$estimates)
    )
    se <- stats::setNames(
      as.numeric(bootstrap_fit$se[1, ]), colnames(bootstrap_fit$se)
    )
    fit$details$full_bootstrap_covariance <- bootstrap_fit$covariance
    fit$details$bootstrap_attempts <- bootstrap_fit$attempts
    fit$details$bootstrap_failures <- bootstrap_fit$failures
  } else if (nboot > 0) {
    boot <- replicate(nboot, {
      idx <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
      out <- try(.logistic_update_once(data[idx, , drop = FALSE], y, x, x_star, complete, strata,
                                       beta_interval, h, tune, estimate_kernel), silent = TRUE)
      if (inherits(out, "try-error")) rep(NA_real_, length(fit$estimates)) else out$estimates
    })
    se <- sqrt(apply(boot, 1, .safe_var))
    names(se) <- names(fit$estimates)
  }

  .new_twostep_fit(call, "logistic", fit$estimates, se, fit$details)
}

.logistic_update_once <- function(data, y, x, x_star, complete, strata,
                                  beta_interval, h, tune, estimate_kernel,
                                  compute_update = TRUE) {
  .check_required_columns(data, c(y, x, x_star, complete, strata))
  data_compl <- .complete_data(data, complete)
  if (nrow(data_compl) < 2) {
    stop("At least two complete observations are required.", call. = FALSE)
  }

  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[data[[complete]] == 1]
  f_x <- stats::as.formula(paste(y, "~ 0 +", x))
  f_xstar <- stats::as.formula(paste(y, "~ 0 +", x_star))

  res_compl <- stats::glm(f_x, data = data_compl, family = "binomial", weights = ipw_compl)
  est_compl <- unname(stats::coef(res_compl)[1])
  res_star <- stats::glm(f_xstar, data = data, family = "binomial")
  res_star_compl <- stats::glm(f_xstar, data = data_compl, family = "binomial", weights = ipw_compl)
  est_star <- unname(stats::coef(res_star)[1])
  est_star_compl <- unname(stats::coef(res_star_compl)[1])

  details <- list(ipw = ipw_all, auxiliary_estimates = c(full = est_star, complete = est_star_compl),
                  h = h, tune = tune)
  differences <- list(default = est_star_compl - est_star)

  if (estimate_kernel) {
    if (is.null(h)) {
      h <- 1.06 * stats::sd(data_compl[[x_star]]) * (nrow(data_compl)^(-1 / 5))
    }
    if (is.null(beta_interval)) {
      beta_interval <- est_compl + c(-0.3, 0.3)
    }
    est_star_kernel <- stats::uniroot(
      function(beta) .logistic_score_kernel(
        beta, data, data_compl, h, ipw_compl, NULL, tune, y, x, x_star
      ),
      beta_interval, extendInt = "yes"
    )$root
    est_star_compl_kernel <- stats::uniroot(
      function(beta) .logistic_score_kernel(
        beta, data_compl, data_compl, h, ipw_compl, ipw_compl,
        tune, y, x, x_star
      ),
      beta_interval, extendInt = "yes"
    )$root
    details$kernel_auxiliary_estimates <- c(full = est_star_kernel, complete = est_star_compl_kernel)
    details$h <- h
    differences$kernel <- est_star_compl_kernel - est_star_kernel
  }

  components <- .bootstrap_components(est_compl, differences)
  if (!compute_update) {
    return(list(estimates = NULL, details = details, components = components))
  }

  p_compl <- stats::plogis(est_compl * data_compl[[x]])
  p_star <- stats::plogis(est_star * data[[x_star]])
  p_star_compl <- stats::plogis(est_star * data_compl[[x_star]])
  info_compl <- sum(ipw_compl * data_compl[[x]]^2 * p_compl * (1 - p_compl))
  info_logis <- sum(data[[x_star]]^2 * p_star * (1 - p_star))
  sigma12_mid <- sum(ipw_compl * (ipw_compl - 1) * data_compl[[x]] *
                       (data_compl[[y]] - p_compl) *
                       data_compl[[x_star]] *
                       (data_compl[[y]] - p_star_compl))
  sigma22_mid <- sum(((data[[complete]] * ipw_all - 1) * data[[x_star]] *
                        (data[[y]] - p_star))^2)
  sigma12 <- sigma12_mid / (info_compl * info_logis)
  sigma22 <- sigma22_mid / (info_logis^2)
  estimates <- c(
    complete_case = est_compl,
    default = est_compl - sigma12 / sigma22 * differences$default
  )

  if (estimate_kernel) {
    info_kernel <- -.finite_grad(
      function(beta) .logistic_score_kernel(
        beta, data, data_compl, h, ipw_compl, NULL, NULL, y, x, x_star
      ),
      est_star_kernel
    )
    sigma12_kernel <- .logistic_score2_org_kernel(est_compl, est_star_kernel, data_compl, data_compl,
                                                  h, ipw_compl, ipw_compl,
                                                  y, x, x_star, complete) /
      (info_compl * info_kernel)
    sigma22_kernel <- .logistic_score2_kernel(est_star_kernel, data, data_compl, h,
                                              ipw_compl, ipw_all,
                                              y, x, x_star, complete) /
      (info_kernel^2)
    estimates <- c(
      estimates,
      kernel = est_compl - sigma12_kernel / sigma22_kernel * differences$kernel
    )
  }

  list(estimates = estimates, details = details, components = components)
}
