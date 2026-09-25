.linear_score_kernel <- function(beta, data, data_compl, h, sdy, ipw_compl,
                                 ipw = NULL,
                                 y = "Y", x = "X", x_star = "X_star") {
  res <- numeric(nrow(data))
  xc <- data_compl[[x]]
  xsc <- data_compl[[x_star]]
  for (i in seq_len(nrow(data))) {
    score_org <- xc * (data[[y]][i] - beta * xc)
    log_fY_X <- stats::dnorm(data[[y]][i], beta * xc, sdy, log = TRUE)
    log_weight <- .log_kernel(data[[x_star]][i] - xsc, h)
    fY_X <- exp(log_fY_X - max(log_fY_X))
    weight <- ipw_compl * exp(log_weight - max(log_weight))
    res[i] <- sum(score_org * fY_X * weight) / sum(fY_X * weight)
  }
  if (!is.null(ipw)) {
    res <- ipw * res
  }
  sum(res)
}

.linear_score2_kernel <- function(beta, data, data_compl, h, sdy, ipw_compl,
                                  ipw,
                                  y = "Y", x = "X", x_star = "X_star",
                                  complete = "R") {
  res <- numeric(nrow(data))
  xc <- data_compl[[x]]
  xsc <- data_compl[[x_star]]
  for (i in seq_len(nrow(data))) {
    score_org <- xc * (data[[y]][i] - beta * xc)
    log_fY_X <- stats::dnorm(data[[y]][i], beta * xc, sdy, log = TRUE)
    log_weight <- .log_kernel(data[[x_star]][i] - xsc, h)
    fY_X <- exp(log_fY_X - max(log_fY_X))
    weight <- ipw_compl * exp(log_weight - max(log_weight))
    score_worki <- sum(score_org * fY_X * weight) / sum(fY_X * weight)
    res[i] <- ((data[[complete]][i] * ipw[i] - 1) * score_worki)^2
  }
  sum(res)
}

.linear_score2_org_kernel <- function(beta_true, beta_work, data, data_compl,
                                      h, sdy, ipw_compl, ipw,
                                      y = "Y", x = "X", x_star = "X_star",
                                      complete = "R") {
  res <- numeric(nrow(data))
  xc <- data_compl[[x]]
  xsc <- data_compl[[x_star]]
  for (i in seq_len(nrow(data))) {
    score_org <- xc * (data[[y]][i] - beta_work * xc)
    log_fY_X <- stats::dnorm(data[[y]][i], beta_work * xc, sdy, log = TRUE)
    log_weight <- .log_kernel(data[[x_star]][i] - xsc, h)
    fY_X <- exp(log_fY_X - max(log_fY_X))
    weight <- ipw_compl * exp(log_weight - max(log_weight))
    score_worki <- sum(score_org * fY_X * weight) / sum(fY_X * weight)
    score_orgi <- data[[x]][i] * (data[[y]][i] - beta_true * data[[x]][i])
    res[i] <- (data[[complete]][i] * ipw[i] * score_orgi) *
      ((data[[complete]][i] * ipw[i] - 1) * score_worki)
  }
  sum(res)
}

.logistic_score_kernel <- function(beta, data, data_compl, h, ipw_compl,
                                   ipw = NULL,
                                   tune = NULL, y = "Y", x = "X", x_star = "X_star") {
  res <- numeric(nrow(data))
  xc <- data_compl[[x]]
  xsc <- data_compl[[x_star]]
  for (i in seq_len(nrow(data))) {
    eta <- beta * xc
    score_org <- xc * (data[[y]][i] - stats::plogis(eta))
    log_fY_X <- eta * data[[y]][i] - .log1pexp(eta)
    log_weight <- .log_kernel(data[[x_star]][i] - xsc, h)
    fY_X <- exp(log_fY_X - max(log_fY_X))
    weight <- ipw_compl * exp(log_weight - max(log_weight))
    res[i] <- sum(score_org * fY_X * weight) / sum(fY_X * weight)
  }
  if (!is.null(ipw)) {
    res <- ipw * res
  }
  if (is.null(tune)) {
    sum(res)
  } else {
    mean(res) - 2 * beta * (nrow(data)^(-1 / 3)) * tune
  }
}

.logistic_score2_kernel <- function(beta, data, data_compl, h, ipw_compl,
                                    ipw,
                                    y = "Y", x = "X", x_star = "X_star",
                                    complete = "R") {
  res <- numeric(nrow(data))
  xc <- data_compl[[x]]
  xsc <- data_compl[[x_star]]
  for (i in seq_len(nrow(data))) {
    eta <- beta * xc
    score_org <- xc * (data[[y]][i] - stats::plogis(eta))
    log_fY_X <- eta * data[[y]][i] - .log1pexp(eta)
    log_weight <- .log_kernel(data[[x_star]][i] - xsc, h)
    fY_X <- exp(log_fY_X - max(log_fY_X))
    weight <- ipw_compl * exp(log_weight - max(log_weight))
    score_worki <- sum(score_org * fY_X * weight) / sum(fY_X * weight)
    res[i] <- ((data[[complete]][i] * ipw[i] - 1) * score_worki)^2
  }
  sum(res)
}

.logistic_score2_org_kernel <- function(beta_true, beta_work, data, data_compl,
                                        h, ipw_compl, ipw,
                                        y = "Y", x = "X", x_star = "X_star",
                                        complete = "R") {
  res <- numeric(nrow(data))
  xc <- data_compl[[x]]
  xsc <- data_compl[[x_star]]
  for (i in seq_len(nrow(data))) {
    eta <- beta_work * xc
    score_org <- xc * (data[[y]][i] - stats::plogis(eta))
    log_fY_X <- eta * data[[y]][i] - .log1pexp(eta)
    log_weight <- .log_kernel(data[[x_star]][i] - xsc, h)
    fY_X <- exp(log_fY_X - max(log_fY_X))
    weight <- ipw_compl * exp(log_weight - max(log_weight))
    score_worki <- sum(score_org * fY_X * weight) / sum(fY_X * weight)
    score_orgi <- data[[x]][i] * (data[[y]][i] - stats::plogis(beta_true * data[[x]][i]))
    res[i] <- (data[[complete]][i] * ipw[i] * score_orgi) *
      ((data[[complete]][i] * ipw[i] - 1) * score_worki)
  }
  sum(res)
}
