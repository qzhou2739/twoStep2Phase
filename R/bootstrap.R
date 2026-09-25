.validate_var_method <- function(var_method) {
  match.arg(var_method, c("standard", "full_bootstrap"))
}

.bootstrap_components <- function(theta, differences) {
  theta_names <- names(theta)
  theta <- as.numeric(theta)
  names(theta) <- theta_names
  differences <- lapply(differences, as.numeric)
  list(theta = theta, differences = differences)
}

.valid_bootstrap_components <- function(draw, expected) {
  if (!is.list(draw) || !is.numeric(draw$theta) ||
      length(draw$theta) != length(expected$theta) ||
      !identical(names(draw$differences), names(expected$differences)) ||
      any(!is.finite(draw$theta))) {
    return(FALSE)
  }
  expected_lengths <- lengths(expected$differences)
  draw_lengths <- lengths(draw$differences)
  identical(draw_lengths, expected_lengths) &&
    all(vapply(draw$differences, function(x) {
      is.numeric(x) && all(is.finite(x))
    }, logical(1)))
}

.fully_bootstrap_fit <- function(data, nboot, fit_once, observed) {
  largest_difference <- max(lengths(observed$differences))
  if (nboot <= largest_difference) {
    stop(
      sprintf(
        paste0(
          "For `var_method = \"full_bootstrap\"`, `nboot` must be greater ",
          "than the largest auxiliary-vector dimension (%d)."
        ),
        largest_difference
      ),
      call. = FALSE
    )
  }
  draws <- vector("list", nboot)
  completed <- 0L
  attempts <- 0L
  max_attempts <- nboot + max(10L, ceiling(nboot / 4))
  while (completed < nboot && attempts < max_attempts) {
    attempts <- attempts + 1L
    index <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
    draw <- try(fit_once(data[index, , drop = FALSE]), silent = TRUE)
    if (!inherits(draw, "try-error") &&
        .valid_bootstrap_components(draw, observed)) {
      completed <- completed + 1L
      draws[[completed]] <- draw
    }
  }
  if (completed <= largest_difference) {
    stop(
      sprintf(
        paste0(
          "Only %d fully bootstrap fits succeeded; more than %d are required ",
          "to estimate the update covariance."
        ),
        completed, largest_difference
      ),
      call. = FALSE
    )
  }
  if (completed < nboot) {
    warning(
      sprintf(
        "Only %d of %d requested fully bootstrap fits succeeded after %d attempts.",
        completed, nboot, attempts
      ),
      call. = FALSE
    )
  }
  draws <- draws[seq_len(completed)]

  p <- length(observed$theta)
  methods <- names(observed$differences)
  theta_draws <- do.call(rbind, lapply(draws, `[[`, "theta"))
  theta_draws <- matrix(theta_draws, ncol = p)
  sigma11 <- matrix(stats::cov(theta_draws), nrow = p, ncol = p)

  estimates <- matrix(
    observed$theta, nrow = p, ncol = length(methods) + 1L,
    dimnames = list(names(observed$theta), c("complete_case", methods))
  )
  standard_errors <- matrix(
    NA_real_, nrow = p, ncol = length(methods) + 1L,
    dimnames = dimnames(estimates)
  )
  standard_errors[, "complete_case"] <- sqrt(pmax(diag(sigma11), 0))
  covariance <- list(complete_case = sigma11)

  for (method in methods) {
    q <- length(observed$differences[[method]])
    difference_draws <- do.call(
      rbind,
      lapply(draws, function(draw) draw$differences[[method]])
    )
    difference_draws <- matrix(difference_draws, ncol = q)
    sigma12 <- matrix(
      stats::cov(theta_draws, difference_draws), nrow = p, ncol = q
    )
    sigma22 <- matrix(stats::cov(difference_draws), nrow = q, ncol = q)
    sigma22_inv <- try(solve(sigma22), silent = TRUE)
    if (inherits(sigma22_inv, "try-error")) {
      stop(
        sprintf(
          paste0(
            "The fully bootstrap covariance matrix is singular for method ",
            "`%s`. Increase `nboot` or remove a highly correlated working model."
          ),
          method
        ),
        call. = FALSE
      )
    }

    projection <- sigma12 %*% sigma22_inv
    estimates[, method] <- observed$theta - as.numeric(
      projection %*% observed$differences[[method]]
    )
    variance <- sigma11 - projection %*% t(sigma12)
    variance <- (variance + t(variance)) / 2
    variance_diagonal <- diag(variance)
    tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(diag(sigma11))))
    if (any(variance_diagonal < -tolerance)) {
      stop(
        sprintf(
          "The fully bootstrap variance is negative for method `%s`.", method
        ),
        call. = FALSE
      )
    }
    standard_errors[, method] <- sqrt(pmax(variance_diagonal, 0))
    covariance[[method]] <- rbind(
      cbind(sigma11, sigma12),
      cbind(t(sigma12), sigma22)
    )
  }

  list(
    estimates = estimates,
    se = standard_errors,
    covariance = covariance,
    attempts = attempts,
    failures = attempts - completed
  )
}
