#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

namespace {

struct RiskMoments {
  arma::vec s0;
  arma::mat s1;
};

void check_reference(const arma::vec& beta,
                     const arma::vec& reference_time,
                     const arma::mat& reference_covariates,
                     const arma::vec& reference_weights) {
  if (reference_time.n_elem != reference_covariates.n_rows ||
      reference_time.n_elem != reference_weights.n_elem ||
      reference_covariates.n_cols != beta.n_elem) {
    Rcpp::stop("Incompatible Cox reference-data dimensions.");
  }
  if (arma::any(reference_weights <= 0.0)) {
    Rcpp::stop("Cox reference weights must be positive.");
  }
}

RiskMoments risk_moments(const arma::vec& beta,
                         const arma::vec& evaluation_time,
                         const arma::vec& reference_time,
                         const arma::mat& reference_covariates,
                         const arma::vec& reference_weights) {
  const arma::uword n = evaluation_time.n_elem;
  const arma::uword m = reference_time.n_elem;
  const arma::uword p = beta.n_elem;
  const double weight_sum = arma::sum(reference_weights);
  const arma::vec relative_risk = arma::exp(reference_covariates * beta);
  RiskMoments result{arma::zeros<arma::vec>(n),
                     arma::zeros<arma::mat>(n, p)};

  for (arma::uword i = 0; i < n; ++i) {
    double s0 = 0.0;
    arma::rowvec s1 = arma::zeros<arma::rowvec>(p);
    for (arma::uword j = 0; j < m; ++j) {
      if (reference_time[j] >= evaluation_time[i]) {
        const double contribution = reference_weights[j] * relative_risk[j];
        s0 += contribution;
        s1 += contribution * reference_covariates.row(j);
      }
    }
    result.s0[i] = s0 / weight_sum;
    result.s1.row(i) = s1 / weight_sum;
  }
  return result;
}

double baseline_cumulative_hazard(
    double time, const arma::vec& reference_time,
    const arma::vec& reference_event,
    const arma::mat& reference_covariates,
    const arma::vec& estimate,
    const arma::vec& reference_weights) {
  const arma::uword m = reference_time.n_elem;
  const arma::vec relative_risk = arma::exp(reference_covariates * estimate);
  double result = 0.0;
  for (arma::uword k = 0; k < m; ++k) {
    if (reference_event[k] == 1.0 && reference_time[k] <= time) {
      double denominator = 0.0;
      for (arma::uword j = 0; j < m; ++j) {
        if (reference_time[j] >= reference_time[k]) {
          denominator += reference_weights[j] * relative_risk[j];
        }
      }
      result += reference_weights[k] / denominator;
    }
  }
  return result;
}

arma::rowvec event_correction(
    double subject_time, const arma::rowvec& subject_covariates,
    double subject_relative_risk, const arma::vec& all_time,
    const arma::vec& all_event, const RiskMoments& all_moments) {
  const arma::uword p = subject_covariates.n_elem;
  arma::rowvec correction = arma::zeros<arma::rowvec>(p);
  for (arma::uword l = 0; l < all_time.n_elem; ++l) {
    if (all_event[l] == 1.0 && all_time[l] <= subject_time) {
      const arma::rowvec mean = all_moments.s1.row(l) / all_moments.s0[l];
      correction += subject_relative_risk *
        (subject_covariates - mean) / all_moments.s0[l];
    }
  }
  return correction / static_cast<double>(all_time.n_elem);
}

arma::mat conditional_kernel_score_vec(
    const arma::vec& beta, const arma::vec& time,
    const arma::vec& event, const arma::mat& z,
    const arma::vec& all_time, const arma::vec& all_event,
    const arma::vec& complete_time, const arma::vec& complete_event,
    const arma::mat& complete_covariates,
    const arma::mat& auxiliary, const arma::mat& complete_auxiliary,
    double h, const arma::vec& complete_estimate,
    const arma::vec& complete_weights, const arma::vec& ipw) {
  const arma::uword n = time.n_elem;
  const arma::uword m = complete_time.n_elem;
  const arma::uword p = beta.n_elem;
  if (event.n_elem != n || z.n_rows != n || z.n_cols + 1 != p ||
      auxiliary.n_rows != n || complete_covariates.n_rows != m ||
      complete_covariates.n_cols != p || complete_auxiliary.n_rows != m ||
      auxiliary.n_cols != complete_auxiliary.n_cols || ipw.n_elem != n ||
      complete_event.n_elem != m || complete_weights.n_elem != m || h <= 0.0) {
    Rcpp::stop("Incompatible Cox kernel-data dimensions.");
  }

  const RiskMoments moments = risk_moments(
    beta, time, complete_time, complete_covariates, complete_weights
  );
  const RiskMoments all_moments = risk_moments(
    beta, all_time, complete_time, complete_covariates, complete_weights
  );
  arma::mat result = arma::zeros<arma::mat>(n, p);

  for (arma::uword i = 0; i < n; ++i) {
    const double lambda = baseline_cumulative_hazard(
      time[i], complete_time, complete_event, complete_covariates,
      complete_estimate, complete_weights
    );
    const arma::rowvec risk_mean = moments.s1.row(i) / moments.s0[i];
    arma::vec log_weights(m);
    arma::mat scores(m, p);

    for (arma::uword j = 0; j < m; ++j) {
      arma::rowvec candidate(p);
      candidate[0] = complete_covariates(j, 0);
      if (p > 1) candidate.cols(1, p - 1) = z.row(i);
      const double eta = arma::dot(candidate, beta);
      const arma::rowvec correction = event_correction(
        time[i], candidate, std::exp(eta), all_time, all_event, all_moments
      );
      scores.row(j) = event[i] * (candidate - risk_mean) - correction;
      const double distance_squared = arma::accu(
        arma::square(auxiliary.row(i) - complete_auxiliary.row(j))
      );
      log_weights[j] = event[i] * eta - lambda * std::exp(eta) -
        distance_squared / (2.0 * h * h) + std::log(complete_weights[j]);
    }
    const double maximum = log_weights.max();
    const arma::vec weights = arma::exp(log_weights - maximum);
    result.row(i) = ipw[i] * (weights.t() * scores) / arma::sum(weights);
  }
  return result;
}

arma::mat conditional_model_score_vec(
    const arma::vec& beta, const arma::vec& time,
    const arma::vec& event, const arma::mat& z,
    const arma::vec& all_time, const arma::vec& all_event,
    const arma::vec& complete_time, const arma::vec& complete_event,
    const arma::mat& complete_covariates,
    const arma::vec& predicted_x, double se_x,
    const arma::vec& grid, const arma::vec& quadrature_weights,
    const arma::vec& complete_estimate,
    const arma::vec& complete_weights, const arma::vec& ipw) {
  const arma::uword n = time.n_elem;
  const arma::uword p = beta.n_elem;
  const arma::uword q = grid.n_elem;
  if (event.n_elem != n || z.n_rows != n || z.n_cols + 1 != p ||
      predicted_x.n_elem != n || ipw.n_elem != n || se_x <= 0.0 ||
      quadrature_weights.n_elem != q || complete_covariates.n_cols != p) {
    Rcpp::stop("Incompatible Cox working-model dimensions.");
  }

  const RiskMoments moments = risk_moments(
    beta, time, complete_time, complete_covariates, complete_weights
  );
  const RiskMoments all_moments = risk_moments(
    beta, all_time, complete_time, complete_covariates, complete_weights
  );
  arma::mat result = arma::zeros<arma::mat>(n, p);
  const double beta_x = beta[0];

  for (arma::uword i = 0; i < n; ++i) {
    const double lambda = baseline_cumulative_hazard(
      time[i], complete_time, complete_event, complete_covariates,
      complete_estimate, complete_weights
    );
    double fixed_eta = 0.0;
    if (p > 1) fixed_eta = arma::dot(z.row(i), beta.subvec(1, p - 1));
    double mode = 0.0;
    double hessian = -1.0 / (se_x * se_x);
    for (int iteration = 0; iteration < 100; ++iteration) {
      const double old_mode = mode;
      const double eta = beta_x * mode + fixed_eta;
      const double score = event[i] * beta_x -
        lambda * std::exp(eta) * beta_x -
        (mode - predicted_x[i]) / (se_x * se_x);
      hessian = -lambda * std::exp(eta) * beta_x * beta_x -
        1.0 / (se_x * se_x);
      mode -= score / hessian;
      if (std::abs(mode - old_mode) < 1e-3) break;
    }

    const arma::vec candidate_x = grid / std::sqrt(-hessian) + mode;
    const arma::rowvec risk_mean = moments.s1.row(i) / moments.s0[i];
    arma::vec log_weights(q);
    arma::mat scores(q, p);
    for (arma::uword j = 0; j < q; ++j) {
      arma::rowvec candidate(p);
      candidate[0] = candidate_x[j];
      if (p > 1) candidate.cols(1, p - 1) = z.row(i);
      const double eta = arma::dot(candidate, beta);
      const arma::rowvec correction = event_correction(
        time[i], candidate, std::exp(eta), all_time, all_event, all_moments
      );
      scores.row(j) = event[i] * (candidate - risk_mean) - correction;
      const double normal_part = -0.5 * std::pow(
        (candidate_x[j] - predicted_x[i]) / se_x, 2.0
      );
      log_weights[j] = event[i] * eta - lambda * std::exp(eta) +
        normal_part + grid[j] * grid[j] + std::log(quadrature_weights[j]);
    }
    const double maximum = log_weights.max();
    const arma::vec weights = arma::exp(log_weights - maximum);
    result.row(i) = ipw[i] * (weights.t() * scores) / arma::sum(weights);
  }
  return result;
}

} // namespace

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_general_info(
    arma::vec beta, arma::vec evaluation_time, arma::vec evaluation_event,
    arma::vec reference_time, arma::mat reference_covariates,
    arma::vec reference_weights) {
  check_reference(beta, reference_time, reference_covariates, reference_weights);
  if (evaluation_time.n_elem != evaluation_event.n_elem) {
    Rcpp::stop("Incompatible Cox evaluation-data dimensions.");
  }
  const arma::uword p = beta.n_elem;
  arma::mat result = arma::zeros<arma::mat>(p, p);
  const arma::vec relative_risk = arma::exp(reference_covariates * beta);

  for (arma::uword i = 0; i < evaluation_time.n_elem; ++i) {
    if (evaluation_event[i] != 1.0) continue;
    double s0 = 0.0;
    arma::rowvec s1 = arma::zeros<arma::rowvec>(p);
    arma::mat s2 = arma::zeros<arma::mat>(p, p);
    for (arma::uword j = 0; j < reference_time.n_elem; ++j) {
      if (reference_time[j] >= evaluation_time[i]) {
        const double contribution = reference_weights[j] * relative_risk[j];
        s0 += contribution;
        s1 += contribution * reference_covariates.row(j);
        s2 += contribution * reference_covariates.row(j).t() *
          reference_covariates.row(j);
      }
    }
    result += s2 / s0 - s1.t() * s1 / (s0 * s0);
  }
  return result;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_general_score_vec(
    arma::vec beta, arma::vec time, arma::vec event,
    arma::mat covariates, arma::vec all_time, arma::vec all_event,
    arma::vec reference_time, arma::mat reference_covariates,
    arma::vec reference_weights, arma::vec ipw) {
  check_reference(beta, reference_time, reference_covariates, reference_weights);
  if (time.n_elem != event.n_elem || time.n_elem != covariates.n_rows ||
      covariates.n_cols != beta.n_elem || ipw.n_elem != time.n_elem ||
      all_time.n_elem != all_event.n_elem) {
    Rcpp::stop("Incompatible Cox score-data dimensions.");
  }
  const RiskMoments moments = risk_moments(
    beta, time, reference_time, reference_covariates, reference_weights
  );
  const RiskMoments all_moments = risk_moments(
    beta, all_time, reference_time, reference_covariates, reference_weights
  );
  arma::mat result = arma::zeros<arma::mat>(time.n_elem, beta.n_elem);
  for (arma::uword i = 0; i < time.n_elem; ++i) {
    const arma::rowvec covariate = covariates.row(i);
    const arma::rowvec mean = moments.s1.row(i) / moments.s0[i];
    const double relative_risk = std::exp(arma::dot(covariate, beta));
    result.row(i) = ipw[i] * (
      event[i] * (covariate - mean) -
      event_correction(
        time[i], covariate, relative_risk,
        all_time, all_event, all_moments
      )
    );
  }
  return result;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_general_kernel_score_vec(
    arma::vec beta, arma::vec time, arma::vec event, arma::mat z,
    arma::vec all_time, arma::vec all_event,
    arma::vec complete_time, arma::vec complete_event,
    arma::mat complete_covariates, arma::mat auxiliary,
    arma::mat complete_auxiliary, double h, arma::vec complete_estimate,
    arma::vec complete_weights, arma::vec ipw) {
  return conditional_kernel_score_vec(
    beta, time, event, z, all_time, all_event,
    complete_time, complete_event, complete_covariates,
    auxiliary, complete_auxiliary, h, complete_estimate,
    complete_weights, ipw
  );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cpp_cox_general_model_score_vec(
    arma::vec beta, arma::vec time, arma::vec event, arma::mat z,
    arma::vec all_time, arma::vec all_event,
    arma::vec complete_time, arma::vec complete_event,
    arma::mat complete_covariates, arma::vec predicted_x, double se_x,
    arma::vec grid, arma::vec quadrature_weights,
    arma::vec complete_estimate, arma::vec complete_weights,
    arma::vec ipw) {
  return conditional_model_score_vec(
    beta, time, event, z, all_time, all_event,
    complete_time, complete_event, complete_covariates,
    predicted_x, se_x, grid, quadrature_weights,
    complete_estimate, complete_weights, ipw
  );
}
