.validate_generator_inputs <- function(N, n, beta, mu, Sigma, sde) {
  if (!is.numeric(N) || length(N) != 1L || is.na(N) ||
      N != as.integer(N) || N < 2L) {
    stop("`N` must be a whole number of at least 2.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != as.integer(n) || n < 1L || n > N) {
    stop("`n` must be a whole number between 1 and `N`.", call. = FALSE)
  }
  if (!is.numeric(beta) || length(beta) < 1L || anyNA(beta) ||
      any(!is.finite(beta))) {
    stop("`beta` must be a finite numeric vector.", call. = FALSE)
  }
  p <- length(beta)
  Sigma <- as.matrix(Sigma)
  if (!is.numeric(mu) || length(mu) != p || anyNA(mu) ||
      any(!is.finite(mu))) {
    stop(sprintf("`mu` must be a finite numeric vector of length %d.", p),
         call. = FALSE)
  }
  if (!is.numeric(Sigma) || !identical(dim(Sigma), c(p, p)) ||
      anyNA(Sigma) || any(!is.finite(Sigma))) {
    stop(sprintf("`Sigma` must be a finite %d x %d covariance matrix.", p, p),
         call. = FALSE)
  }
  if (!isTRUE(all.equal(Sigma, t(Sigma)))) {
    stop("`Sigma` must be symmetric.", call. = FALSE)
  }
  tryCatch(chol(Sigma), error = function(e) {
    stop("`Sigma` must be positive definite.", call. = FALSE)
  })
  if (!is.numeric(sde) || length(sde) != 1L ||
      is.na(sde) || !is.finite(sde) || sde < 0) {
    stop("`sde` must be one non-negative finite standard deviation.", call. = FALSE)
  }
  invisible(TRUE)
}

.generate_covariates <- function(N, beta, mu, Sigma, sde) {
  covariates <- .rmvnorm_chol(N, mu, Sigma)
  X <- covariates[, 1L]
  Z <- if (length(beta) > 1L) {
    covariates[, -1L, drop = FALSE]
  } else {
    matrix(numeric(), nrow = N, ncol = 0L)
  }
  X_star <- X + stats::rnorm(N, mean = 0, sd = sde)
  list(X = X, Z = Z, X_star = X_star)
}

.add_z_columns <- function(data, Z) {
  if (ncol(Z) == 0L) return(data)
  z_data <- as.data.frame(Z)
  names(z_data) <- paste0("Z", seq_len(ncol(Z)))
  cbind(data, z_data)
}

.linear_strata <- function(Y, quantile) {
  if (is.null(quantile)) return(rep(0L, length(Y)))
  if (!is.numeric(quantile) || length(quantile) < 1L || anyNA(quantile) ||
      any(!is.finite(quantile)) || any(quantile <= 0 | quantile >= 1) ||
      is.unsorted(quantile, strictly = TRUE)) {
    stop(
      "`quantile` must be NULL or strictly increasing probabilities between 0 and 1.",
      call. = FALSE
    )
  }
  cut_points <- stats::quantile(Y, probs = quantile, names = FALSE)
  findInterval(Y, cut_points)
}

.validate_generator_ns <- function(ns, n) {
  if (!is.numeric(ns) || length(ns) < 1L || anyNA(ns) ||
      any(!is.finite(ns)) || any(ns < 0) || any(ns != as.integer(ns))) {
    stop("`ns` must contain non-negative whole-number sample sizes.", call. = FALSE)
  }
  if (sum(ns) != n) {
    stop("The entries of `ns` must sum to `n`.", call. = FALSE)
  }
  as.integer(ns)
}

.binary_sampling_strata <- function(indicator, ns) {
  if (length(ns) == 1L) return(rep(0L, length(indicator)))
  if (length(ns) == 2L) return(as.integer(indicator))
  stop(
    "For logistic and Cox data, `ns` must have length 1 for MCAR sampling or length 2 for outcome/event-stratified sampling.",
    call. = FALSE
  )
}

generate_linear_data <- function(N, n, beta, mu, Sigma, sde,
                                 sdy = 0.5, ns = n, quantile = NULL) {
  .validate_generator_inputs(N, n, beta, mu, Sigma, sde)
  if (!is.numeric(sdy) || length(sdy) != 1L || is.na(sdy) ||
      !is.finite(sdy) || sdy <= 0) {
    stop("`sdy` must be one positive finite standard deviation.", call. = FALSE)
  }
  generated <- .generate_covariates(N, beta, mu, Sigma, sde)
  design <- cbind(generated$X, generated$Z)
  Y <- as.vector(design %*% beta) + stats::rnorm(N, mean = 0, sd = sdy)
  stratum <- .linear_strata(Y, quantile)
  ns <- .validate_generator_ns(ns, n)
  R <- .sample_two_phase(stratum, ns, decreasing = FALSE)
  generated$X[R == 0] <- NA_real_
  out <- data.frame(
    R = R, Y = Y, X = generated$X,
    X_star = generated$X_star
  )
  out <- .add_z_columns(out, generated$Z)
  out$stratum <- stratum
  out
}

generate_logistic_data <- function(N, n, beta, mu, Sigma, sde, ns = n) {
  .validate_generator_inputs(N, n, beta, mu, Sigma, sde)
  generated <- .generate_covariates(N, beta, mu, Sigma, sde)
  design <- cbind(generated$X, generated$Z)
  Y <- stats::rbinom(N, size = 1, prob = stats::plogis(as.vector(design %*% beta)))
  ns <- .validate_generator_ns(ns, n)
  stratum <- .binary_sampling_strata(Y, ns)
  R <- .sample_two_phase(stratum, ns)
  generated$X[R == 0] <- NA_real_
  out <- data.frame(
    R = R, Y = Y, X = generated$X,
    X_star = generated$X_star
  )
  out <- .add_z_columns(out, generated$Z)
  out$stratum <- stratum
  out
}

generate_cox_data <- function(N, n, beta, mu, Sigma, sde,
                              tau = 2.5, ns = n) {
  .validate_generator_inputs(N, n, beta, mu, Sigma, sde)
  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) ||
      !is.finite(tau) || tau <= 0) {
    stop("`tau` must be one positive finite number.", call. = FALSE)
  }
  generated <- .generate_covariates(N, beta, mu, Sigma, sde)
  design <- cbind(generated$X, generated$Z)
  lp <- as.vector(design %*% beta)
  T_true <- sqrt(-log1p(-stats::runif(N)) * exp(-lp) * 4)
  C <- pmin(stats::runif(N, min = 0, max = 5 * tau / 3), tau)
  Y <- pmin(T_true, C)
  D <- as.numeric(T_true <= C)
  ns <- .validate_generator_ns(ns, n)
  stratum <- .binary_sampling_strata(D, ns)
  R <- .sample_two_phase(stratum, ns)
  generated$X[R == 0] <- NA_real_
  out <- data.frame(
    R = R, Y = Y, D = D, X = generated$X,
    X_star = generated$X_star
  )
  out <- .add_z_columns(out, generated$Z)
  out$stratum <- stratum
  out
}
