.cox_optimal_matrix <- function(data, time, event, x = NULL, x_star = NULL) {
  n <- nrow(data)
  x_value <- if (is.null(x)) numeric(n) else data[[x]]
  x_star_value <- if (is.null(x_star)) numeric(n) else data[[x_star]]
  out <- cbind(
    complete = numeric(n), time = data[[time]], event = data[[event]],
    x = x_value, x_star = x_star_value
  )
  storage.mode(out) <- "double"
  out
}

.cox_info_org <- function(beta, data_all, data_compl, ipw_compl,
                          time, event, x) {
  cpp_cox_opt_info_org(
    beta,
    .cox_optimal_matrix(data_all, time, event),
    .cox_optimal_matrix(data_compl, time, event, x = x),
    ipw_compl
  )
}

.cox_info_work <- function(beta, data_all, time, event, x_star) {
  cpp_cox_opt_info_work_coxph(
    beta,
    .cox_optimal_matrix(data_all, time, event, x_star = x_star)
  )
}

.cox_score_org_vec <- function(beta, data, data_all, data_compl, ipw_compl, ipw,
                               time, event, x) {
  as.numeric(cpp_cox_opt_score_org_vec(
    beta,
    .cox_optimal_matrix(data, time, event, x = x),
    .cox_optimal_matrix(data_all, time, event),
    .cox_optimal_matrix(data_compl, time, event, x = x),
    ipw_compl, ipw
  ))
}

.cox_score_coxph_vec <- function(beta, data, data_all, ipw,
                                 time, event, x_star) {
  as.numeric(cpp_cox_opt_score_work_coxph_vec(
    beta,
    .cox_optimal_matrix(data, time, event, x_star = x_star),
    .cox_optimal_matrix(data_all, time, event, x_star = x_star),
    ipw
  ))
}

cox_score_coxph <- function(beta, data, data_all = data, ipw = NULL,
                            time = "Y", event = "D", x_star = "X_star") {
  .check_required_columns(data, c(time, event, x_star))
  .check_required_columns(data_all, c(time, event, x_star))
  if (is.null(ipw)) ipw <- rep(1, nrow(data))
  sum(.cox_score_coxph_vec(beta, data, data_all, ipw, time, event, x_star))
}

.cox_score_kernel_vec <- function(beta, data, data_all, data_compl, h, est_compl,
                                  ipw_compl, ipw, time, event, x, x_star) {
  as.numeric(cpp_cox_opt_score_work_kernel_vec(
    beta,
    .cox_optimal_matrix(data, time, event, x_star = x_star),
    .cox_optimal_matrix(data_all, time, event, x_star = x_star),
    .cox_optimal_matrix(data_compl, time, event, x = x, x_star = x_star),
    h, pi, est_compl, ipw_compl, ipw
  ))
}

cox_score_kernel <- function(beta, data, data_all, data_compl, h, est_compl,
                             ipw_compl, ipw = NULL, time = "Y", event = "D",
                             x = "X", x_star = "X_star") {
  .check_required_columns(data, c(time, event, x_star))
  .check_required_columns(data_all, c(time, event, x_star))
  .check_required_columns(data_compl, c(time, event, x, x_star))
  if (is.null(ipw)) ipw <- rep(1, nrow(data))
  sum(.cox_score_kernel_vec(
    beta, data, data_all, data_compl, h, est_compl,
    ipw_compl, ipw, time, event, x, x_star
  ))
}
