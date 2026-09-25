test_that("Cox estimators use registered native routines", {
  dll <- getLoadedDLLs()[["twoStep2Phase"]]
  expect_false(dll[["dynamicLookup"]])
  expect_s3_class(
    getNativeSymbolInfo(
      "_twoStep2Phase_cpp_cox_opt_score_work_kernel_vec",
      PACKAGE = "twoStep2Phase"
    ),
    "NativeSymbolInfo"
  )
  expect_s3_class(
    getNativeSymbolInfo(
      "_twoStep2Phase_cpp_cox_joint_score_work_modelX_vec",
      PACKAGE = "twoStep2Phase"
    ),
    "NativeSymbolInfo"
  )
  expect_s3_class(
    getNativeSymbolInfo(
      "_twoStep2Phase_cpp_cox_general_kernel_score_vec",
      PACKAGE = "twoStep2Phase"
    ),
    "NativeSymbolInfo"
  )
})

test_that("compiled Cox paths return finite estimates", {
  set.seed(25)
  no_z <- generate_cox_data(
    80, 40, beta = 0.5, mu = 0, Sigma = matrix(1), sde = 1
  )
  fit_no_z <- update_cox(no_z, z = NULL, x_models = "kernel", nboot = 0)
  expect_identical(
    names(fit_no_z$estimate),
    c("complete_case", "default", "kernel")
  )
  expect_true(all(is.finite(fit_no_z$estimate)))

  set.seed(26)
  with_z <- generate_cox_data(
    100, 50, beta = c(0.5, 0.25), mu = c(0, 0),
    Sigma = diag(2), sde = 1
  )
  set.seed(27)
  fit_with_z <- update_cox(
    with_z, z = "Z1", x_models = "linear", n_quad = 5, nboot = 0
  )
  expect_identical(
    colnames(fit_with_z$estimate),
    c("complete_case", "default", "linear", "joint_linear")
  )
  expect_true(all(is.finite(fit_with_z$estimate)))
})
