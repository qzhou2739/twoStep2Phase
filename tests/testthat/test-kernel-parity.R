test_that("linear kernel update uses IPW in the reference sample", {
  set.seed(1)
  N <- 1000
  X <- stats::rnorm(N)
  Y <- 0.5 * X + stats::rnorm(N, sd = 0.5)
  X_star <- X + stats::rnorm(N, sd = 3.1)
  ind_stra <- as.numeric(Y > stats::quantile(Y, 0.7))
  selected_high <- numeric(sum(ind_stra == 1))
  selected_low <- numeric(sum(ind_stra == 0))
  selected_high[sample.int(length(selected_high), 280)] <- 1
  selected_low[sample.int(length(selected_low), 120)] <- 1
  ind_compl <- numeric(N)
  ind_compl[ind_stra == 1] <- selected_high
  ind_compl[ind_stra == 0] <- selected_low
  data <- data.frame(ind_compl, Y, X, X_star, ind_stra)

  fit <- update_linear(
    data, z = NULL, complete = "ind_compl", stratum = "ind_stra",
    x_models = "kernel", nboot = 0
  )

  expect_gt(fit$details$kernel_auxiliary_estimates[["full"]], 0)
  expect_equal(unname(fit$estimate[["kernel"]]), 0.5444516, tolerance = 5e-6)
})

test_that("logistic kernel update uses IPW in the reference sample", {
  set.seed(2)
  N <- 1000
  X <- stats::rnorm(N)
  Y <- stats::rbinom(N, 1, stats::plogis(0.5 * X))
  X_star <- X + stats::rnorm(N)
  ind_stra <- as.numeric(X_star > stats::quantile(X_star, 0.7))
  selected_high <- numeric(sum(ind_stra == 1))
  selected_low <- numeric(sum(ind_stra == 0))
  selected_high[sample.int(length(selected_high), 280)] <- 1
  selected_low[sample.int(length(selected_low), 120)] <- 1
  ind_compl <- numeric(N)
  ind_compl[ind_stra == 1] <- selected_high
  ind_compl[ind_stra == 0] <- selected_low
  data <- data.frame(ind_compl, Y, X, X_star, ind_stra)

  fit <- suppressWarnings(update_logistic(
    data, z = NULL, complete = "ind_compl", stratum = "ind_stra",
    x_models = "kernel", tune_para = 0.01, nboot = 0
  ))

  expect_gt(fit$details$kernel_auxiliary_estimates[["full"]], 0)
  expect_equal(unname(fit$estimate[["kernel"]]), 0.5078009, tolerance = 5e-6)
})
