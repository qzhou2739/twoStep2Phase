.update_cox_optimal <- function(data, time = "Y", event = "D", x = "X",
                               x_star = "X_star", complete = "R",
                               strata = "ind_stra", beta_interval = NULL,
                               h = NULL, estimate_kernel = TRUE, nboot = 0,
                               se_method = "ordinary") {
  call <- match.call()
  fully_bootstrap <- identical(se_method, "fully_bootstrap")
  fit <- .cox_update_once(data, time, event, x, x_star, complete, strata,
                          beta_interval, h, estimate_kernel,
                          compute_update = !fully_bootstrap)

  se <- NULL
  if (fully_bootstrap) {
    bootstrap_fit <- .fully_bootstrap_fit(
      data, nboot,
      function(boot_data) {
        .cox_update_once(
          boot_data, time, event, x, x_star, complete, strata,
          beta_interval, h, estimate_kernel, compute_update = FALSE
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
      out <- try(.cox_update_once(data[idx, , drop = FALSE], time, event, x, x_star,
                                  complete, strata, beta_interval, h, estimate_kernel),
                 silent = TRUE)
      if (inherits(out, "try-error")) rep(NA_real_, length(fit$estimates)) else out$estimates
    })
    se <- sqrt(apply(boot, 1, .safe_var))
    names(se) <- names(fit$estimates)
  }

  .new_twostep_fit(call, "cox", fit$estimates, se, fit$details)
}

.cox_update_once <- function(data, time, event, x, x_star, complete, strata,
                             beta_interval, h, estimate_kernel,
                             compute_update = TRUE) {
  .check_required_columns(data, c(time, event, x, x_star, complete, strata))
  data_compl <- .complete_data(data, complete)
  if (nrow(data_compl) < 2 || sum(data_compl[[event]] == 1) == 0) {
    stop("At least two complete observations and one complete event are required.", call. = FALSE)
  }

  ipw_all <- .ipw_by_stratum(data, complete, strata)
  ipw_compl <- ipw_all[data[[complete]] == 1]
  f_x <- stats::as.formula(paste0("Surv(", time, ", ", event, ") ~ ", x))
  f_xstar <- stats::as.formula(paste0("Surv(", time, ", ", event, ") ~ ", x_star))

  res_compl <- survival::coxph(f_x, data = data_compl, weights = ipw_compl)
  est_compl <- unname(stats::coef(res_compl)[1])
  res_star <- survival::coxph(f_xstar, data = data)
  res_star_compl <- survival::coxph(f_xstar, data = data_compl, weights = ipw_compl)
  est_star <- unname(stats::coef(res_star)[1])
  est_star_compl <- unname(stats::coef(res_star_compl)[1])

  details <- list(ipw = ipw_all, auxiliary_estimates = c(full = est_star, complete = est_star_compl),
                  h = h)
  differences <- list(default = est_star_compl - est_star)

  if (estimate_kernel) {
    if (is.null(h)) {
      h <- 1.06 * stats::sd(data_compl[[x_star]]) * (nrow(data_compl)^(-1 / 5))
    }
    if (is.null(beta_interval)) {
      beta_interval <- est_compl + c(-0.3, 0.3)
    }
    est_star_kernel <- stats::uniroot(
      function(beta) cox_score_kernel(beta, data, data, data_compl, h, est_compl,
                                      ipw_compl, rep(1, nrow(data)), time, event, x, x_star),
      beta_interval, extendInt = "yes"
    )$root
    est_star_compl_kernel <- stats::uniroot(
      function(beta) cox_score_kernel(beta, data_compl, data, data_compl, h, est_compl,
                                      ipw_compl, ipw_compl, time, event, x, x_star),
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

  info_compl <- .cox_info_org(est_compl, data, data_compl, ipw_compl, time, event, x)
  info_coxph <- .cox_info_work(est_star, data, time, event, x_star)
  complete_index <- data[[complete]] == 1
  score_vec_org <- numeric(nrow(data))
  score_vec_org[complete_index] <- .cox_score_org_vec(
    est_compl, data_compl, data, data_compl, ipw_compl,
    rep(1, nrow(data_compl)), time, event, x
  )
  score_vec_coxph <- .cox_score_coxph_vec(est_star, data, data, rep(1, nrow(data)),
                                          time, event, x_star)
  centered <- data[[complete]] * ipw_all - 1
  sigma12_mid <- sum(data[[complete]] * ipw_all * score_vec_org * centered * score_vec_coxph)
  sigma22_mid <- sum((centered * score_vec_coxph)^2)
  sigma12 <- sigma12_mid / (info_compl * info_coxph)
  sigma22 <- sigma22_mid / (info_coxph^2)
  estimates <- c(
    complete_case = est_compl,
    default = est_compl - sigma12 / sigma22 * differences$default
  )

  if (estimate_kernel) {
    info_kernel <- -.finite_grad(
      function(beta) cox_score_kernel(beta, data, data, data_compl, h, est_compl,
                                      ipw_compl, rep(1, nrow(data)), time, event, x, x_star),
      est_star_kernel
    )
    score_vec_kernel <- .cox_score_kernel_vec(est_star_kernel, data, data, data_compl,
                                              h, est_compl, ipw_compl, rep(1, nrow(data)),
                                              time, event, x, x_star)
    sigma12_mid_kernel <- sum(data[[complete]] * ipw_all * score_vec_org * centered * score_vec_kernel)
    sigma22_mid_kernel <- sum((centered * score_vec_kernel)^2)
    sigma12_kernel <- sigma12_mid_kernel / (info_compl * info_kernel)
    sigma22_kernel <- sigma22_mid_kernel / (info_kernel^2)
    estimates <- c(
      estimates,
      kernel = est_compl - sigma12_kernel / sigma22_kernel * differences$kernel
    )
  }

  list(estimates = estimates, details = details, components = components)
}
