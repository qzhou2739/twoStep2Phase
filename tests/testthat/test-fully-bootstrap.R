test_that("fully bootstrap matches the supplied covariance formula", {
  set.seed(201)
  data <- generate_linear_data(
    120, 70, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  observed <- twoStep2Phase:::.linear_update_once(
    data, "Y", "X", "X_star", "R", "stratum",
    NULL, NULL, NULL, FALSE, compute_update = FALSE
  )$components

  set.seed(19)
  draws <- lapply(seq_len(4), function(i) {
    index <- sample(seq_len(nrow(data)), nrow(data), replace = TRUE)
    twoStep2Phase:::.linear_update_once(
      data[index, , drop = FALSE], "Y", "X", "X_star", "R", "stratum",
      NULL, NULL, NULL, FALSE, compute_update = FALSE
    )$components
  })
  bootstrap_values <- do.call(
    rbind,
    lapply(draws, function(draw) {
      c(draw$theta, draw$differences$default)
    })
  )
  covariance <- stats::cov(bootstrap_values)
  expected_estimate <- observed$theta - covariance[1, 2] /
    covariance[2, 2] * observed$differences$default
  expected_se <- c(
    sqrt(covariance[1, 1]),
    sqrt(covariance[1, 1] - covariance[1, 2]^2 / covariance[2, 2])
  )

  fit <- update_linear(
    data, z = NULL, stratum = "stratum", x_models = NULL,
    nboot = 4, seed = 19, var_method = "full_bootstrap"
  )

  expect_equal(unname(fit$estimate["default"]), unname(expected_estimate))
  expect_equal(unname(fit$se), expected_se)
  expect_identical(fit$details$var_method, "full_bootstrap")
  expect_equal(fit$details$bootstrap_failures, 0L)
})

test_that("fully bootstrap supports all three outcome families without Z", {
  set.seed(202)
  linear_data <- generate_linear_data(
    120, 70, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  logistic_data <- generate_logistic_data(
    180, 110, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  cox_data <- generate_cox_data(
    120, 80, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )

  linear_fit <- update_linear(
    linear_data, z = NULL, nboot = 3, seed = 20,
    var_method = "full_bootstrap"
  )
  logistic_fit <- suppressWarnings(update_logistic(
    logistic_data, z = NULL, x_models = NULL, nboot = 3, seed = 21,
    var_method = "full_bootstrap"
  ))
  cox_fit <- update_cox(
    cox_data, z = NULL, x_models = NULL, nboot = 3, seed = 22,
    var_method = "full_bootstrap"
  )

  expect_identical(names(linear_fit$estimate),
                   c("complete_case", "default", "kernel"))
  expect_true(all(is.finite(linear_fit$se)))
  expect_true(all(is.finite(logistic_fit$se)))
  expect_true(all(is.finite(cox_fit$se)))
})

test_that("fully bootstrap supports joint working models", {
  set.seed(203)
  data <- generate_linear_data(
    90, 60, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )
  fit <- update_linear(
    data, x_models = "linear", n_quad = 5, nboot = 6, seed = 23,
    var_method = "full_bootstrap"
  )

  expect_identical(
    colnames(fit$estimate),
    c("complete_case", "default", "linear", "joint_linear")
  )
  expect_true(all(is.finite(fit$se)))
  expect_true(all(fit$details$projection == "full_bootstrap"))
})

test_that("fully bootstrap supports an intercept with a kernel update", {
  set.seed(205)
  data <- generate_linear_data(
    110, 75, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 0.5
  )
  fit <- update_linear(
    data, z = NULL, intercept = TRUE, x_models = "kernel", h = 0.5,
    nboot = 7, seed = 25, var_method = "full_bootstrap"
  )

  expect_identical(rownames(fit$estimate), c("(Intercept)", "X"))
  expect_identical(
    colnames(fit$estimate),
    c("complete_case", "default", "kernel", "joint_kernel")
  )
  expect_true(all(is.finite(fit$se)))
})

test_that("fully bootstrap requires enough replicates for the covariance inverse", {
  set.seed(204)
  data <- generate_linear_data(
    100, 70, beta = c(0.5, 0.25), mu = c(0, 0), Sigma = diag(2),
    sde = 0.5
  )

  expect_error(
    update_linear(
      data, x_models = NULL, nboot = 2, seed = 24,
      var_method = "full_bootstrap"
    ),
    "largest auxiliary-vector dimension (2)", fixed = TRUE
  )
  expect_true("var_method" %in% names(formals(update_linear)))
  expect_true("var_method" %in% names(formals(update_logistic)))
  expect_true("var_method" %in% names(formals(update_cox)))
})
