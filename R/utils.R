.check_required_columns <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.ipw_by_stratum <- function(data, complete = "R", strata = "ind_stra") {
  .check_required_columns(data, c(complete, strata))
  r <- data[[complete]]
  s <- data[[strata]]
  if (any(is.na(r)) || any(is.na(s))) {
    stop("Complete-case and stratum indicators cannot contain missing values.", call. = FALSE)
  }

  ipw <- rep(NA_real_, length(r))
  for (level in unique(s)) {
    idx <- s == level
    phase2 <- sum(r[idx] == 1)
    if (phase2 == 0) {
      stop("At least one stratum has no complete Phase II observations.", call. = FALSE)
    }
    ipw[idx] <- sum(idx) / phase2
  }
  ipw
}

.complete_data <- function(data, complete = "R") {
  data[data[[complete]] == 1, , drop = FALSE]
}

.safe_var <- function(x) {
  if (length(x) < 2 || all(is.na(x))) {
    return(NA_real_)
  }
  stats::var(x, na.rm = TRUE)
}

.finite_grad <- function(f, x, eps = 1e-5) {
  (f(x + eps) - f(x - eps)) / (2 * eps)
}

.log_kernel <- function(t, h) {
  -((t / h)^2) / 2 - log(2 * pi) / 2 - log(h)
}

.log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

.weighted_resid_sd <- function(residuals, weights) {
  sqrt(sum(weights * residuals^2) / sum(weights))
}

.validate_ns <- function(ns, strata_size) {
  if (length(ns) != length(strata_size)) {
    stop("`ns` must contain one sample size per stratum.", call. = FALSE)
  }
  if (any(ns < 0) || any(ns > strata_size)) {
    stop("Each `ns` entry must be between 0 and the corresponding stratum size.", call. = FALSE)
  }
  invisible(TRUE)
}

.sample_two_phase <- function(ind_stra, ns, decreasing = TRUE) {
  strata <- sort(unique(ind_stra), decreasing = decreasing)
  strata_size <- vapply(strata, function(level) sum(ind_stra == level), numeric(1))
  .validate_ns(ns, strata_size)

  ind_compl <- rep(0, length(ind_stra))
  for (j in seq_along(strata)) {
    idx <- which(ind_stra == strata[j])
    if (ns[j] > 0) {
      ind_compl[sample(idx, ns[j])] <- 1
    }
  }
  ind_compl
}

.rmvnorm_chol <- function(n, mu, Sigma) {
  z <- matrix(stats::rnorm(n * length(mu)), nrow = n)
  sweep(z %*% chol(Sigma), 2, mu, "+")
}

.new_twostep_fit <- function(call, family, estimates, se = NULL, details = list()) {
  results <- data.frame(
    estimate = unname(estimates),
    std.error = if (is.null(se)) rep(NA_real_, length(estimates)) else unname(se),
    row.names = names(estimates)
  )
  structure(
    list(call = call, family = family, estimate = estimates, se = se,
         results = results, details = details),
    class = "twostep2phase"
  )
}

print.twostep2phase <- function(x, ...) {
  cat("Two-step two-phase estimator\n")
  cat("Family:", x$family, "\n\n")
  print(x$results)
  invisible(x)
}
