.d_log_normal_pdf <- function(t, m, s) -(t - m) / (s^2)
.dd_log_normal_pdf <- function(s) -1 / (s^2)
.dm_log_normal_pdf <- function(t, m, s) (t - m) / (s^2)
.ddm_log_normal_pdf <- function(s) -1 / (s^2)

linear_score_modelx <- function(beta, data, resX, seX, sdy, grid, weight,
                                ipw = NULL, z = NULL, intercept = FALSE,
                                y = "Y") {
  .linear_score_modelx_impl(
    beta, data, resX, seX, sdy, grid, weight, ipw,
    vector = FALSE, z = z, intercept = intercept, y = y
  )
}

logistic_score_modelx <- function(beta, data, resX, seX, grid, weight, ipw = NULL,
                                  tune = NULL, z = NULL, intercept = FALSE,
                                  y = "Y") {
  res <- .logistic_score_modelx_impl(
    beta, data, resX, seX, grid, weight, ipw, z, intercept, y
  )
  if (is.null(tune)) {
    colSums(res)
  } else {
    colMeans(res) - 2 * beta * (nrow(data)^(-1 / 3)) * tune
  }
}

.z_matrix <- function(data, z = NULL) {
  if (is.null(z)) z <- grep("^Z[0-9]+$", names(data), value = TRUE)
  out <- as.matrix(data[, z, drop = FALSE])
  storage.mode(out) <- "double"
  out
}

.conditional_outcome_matrix <- function(x_values, z_values, intercept = FALSE) {
  x_values <- as.numeric(x_values)
  out <- cbind(x_values)
  if (length(z_values)) {
    out <- cbind(
      out,
      matrix(
        rep(z_values, each = length(x_values)),
        nrow = length(x_values), ncol = length(z_values)
      )
    )
  }
  if (intercept) out <- cbind(1, out)
  storage.mode(out) <- "double"
  out
}

.linear_score_modelx_impl <- function(beta, data, resX, seX, sdy, grid, weight,
                                      ipw, vector, z = NULL,
                                      intercept = FALSE, y = "Y") {
  res <- matrix(0, nrow = nrow(data), ncol = length(beta))
  beta_x <- beta[1L + as.integer(intercept)]
  for (i in seq_len(nrow(data))) {
    Zi <- as.numeric(.z_matrix(data[i, , drop = FALSE], z))
    predXi <- as.numeric(stats::predict(resX, data[i, , drop = FALSE]))
    b <- 0
    for (j in seq_len(100)) {
      old_b <- b
      cova_b <- c(if (intercept) 1, b, Zi)
      score <- .dm_log_normal_pdf(data[[y]][i], sum(cova_b * beta), sdy) * beta_x +
        .d_log_normal_pdf(b, predXi, seX)
      hes <- .ddm_log_normal_pdf(sdy) * (beta_x^2) + .dd_log_normal_pdf(seX)
      b <- b - score / hes
      if (abs(b - old_b) < 1e-3) break
    }
    gridb <- grid / sqrt(-hes) + b
    cova <- .conditional_outcome_matrix(gridb, Zi, intercept)
    score_org <- cova * (data[[y]][i] - as.vector(cova %*% beta))
    log_fY_XZ <- stats::dnorm(data[[y]][i], as.vector(cova %*% beta), sdy, log = TRUE)
    log_fX_XstarZ <- stats::dnorm(gridb, predXi, seX, log = TRUE)
    weight1 <- exp(log_fY_XZ - max(log_fY_XZ) + log_fX_XstarZ - max(log_fX_XstarZ) +
                     grid^2 + log(weight))
    res[i, ] <- colSums(score_org * weight1) / sum(weight1)
  }
  if (!is.null(ipw)) {
    res <- ipw * res
  }
  if (vector) res else colSums(res)
}

.logistic_score_modelx_impl <- function(beta, data, resX, seX, grid, weight,
                                        ipw, z = NULL, intercept = FALSE,
                                        y = "Y") {
  res <- matrix(0, nrow = nrow(data), ncol = length(beta))
  for (i in seq_len(nrow(data))) {
    Zi <- as.numeric(.z_matrix(data[i, , drop = FALSE], z))
    predXi <- as.numeric(stats::predict(resX, data[i, , drop = FALSE]))
    cova <- .conditional_outcome_matrix(grid, Zi, intercept)
    eta <- as.vector(cova %*% beta)
    score_org <- cova * (data[[y]][i] - stats::plogis(eta))
    log_fY_XZ <- eta * data[[y]][i] - .log1pexp(eta)
    log_fX_XstarZ <- stats::dnorm(grid, predXi, seX, log = TRUE)
    weight1 <- exp(log_fY_XZ - max(log_fY_XZ) + log_fX_XstarZ - max(log_fX_XstarZ) +
                     grid^2 + log(weight))
    res[i, ] <- colSums(score_org * weight1) / sum(weight1)
  }
  if (!is.null(ipw)) {
    res <- ipw * res
  }
  res
}

.regression_kernel_score_impl <- function(
    beta, data, data_compl, auxiliary, h, ipw_compl, ipw = NULL,
    family = c("linear", "logistic"), sdy = NULL, tune = NULL,
    y = "Y", x = "X", z = character(), intercept = FALSE,
    vector = FALSE, kernel_scales = NULL) {
  family <- match.arg(family)
  auxiliary_all <- as.matrix(data[, auxiliary, drop = FALSE])
  auxiliary_compl <- as.matrix(data_compl[, auxiliary, drop = FALSE])
  storage.mode(auxiliary_all) <- "double"
  storage.mode(auxiliary_compl) <- "double"
  if (!is.null(kernel_scales)) {
    if (length(kernel_scales) != ncol(auxiliary_all) ||
        any(!is.finite(kernel_scales)) || any(kernel_scales <= 0)) {
      stop("Invalid standardized-kernel scales.", call. = FALSE)
    }
    kernel_scales <- as.numeric(kernel_scales)
  }
  x_compl <- data_compl[[x]]
  res <- matrix(0, nrow = nrow(data), ncol = length(beta))

  for (i in seq_len(nrow(data))) {
    Zi <- as.numeric(.z_matrix(data[i, , drop = FALSE], z))
    cova <- .conditional_outcome_matrix(x_compl, Zi, intercept)
    eta <- as.vector(cova %*% beta)
    if (family == "linear") {
      residual <- data[[y]][i] - eta
      score_org <- cova * residual
      log_outcome <- stats::dnorm(data[[y]][i], eta, sdy, log = TRUE)
    } else {
      score_org <- cova * (data[[y]][i] - stats::plogis(eta))
      log_outcome <- eta * data[[y]][i] - .log1pexp(eta)
    }
    differences <- sweep(
      auxiliary_compl, 2L, auxiliary_all[i, ], FUN = "-"
    )
    if (is.null(kernel_scales)) {
      log_kernel <- -rowSums(differences^2) / (2 * h^2)
    } else {
      standardized_differences <- sweep(
        differences, 2L, kernel_scales, FUN = "/"
      )
      log_kernel <- -rowSums(standardized_differences^2) / (2 * h^2)
    }
    log_weights <- log(ipw_compl) + log_kernel + log_outcome
    weights <- exp(log_weights - max(log_weights))
    res[i, ] <- colSums(score_org * weights) / sum(weights)
  }
  if (!is.null(ipw)) res <- ipw * res
  if (vector) return(res)
  if (family == "logistic" && !is.null(tune)) {
    colMeans(res) - 2 * beta * nrow(data)^(-1 / 3) * tune
  } else {
    colSums(res)
  }
}
