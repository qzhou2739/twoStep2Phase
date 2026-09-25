\name{generate_linear_data}
\alias{generate_linear_data}
\alias{generate_logistic_data}
\alias{generate_cox_data}
\title{Generate Two-Phase Data for Linear, Logistic, and Cox Models}
\usage{
generate_linear_data(N, n, beta, mu, Sigma, sde,
  sdy = 0.5, ns = n, quantile = NULL)

generate_logistic_data(N, n, beta, mu, Sigma, sde, ns = n)

generate_cox_data(N, n, beta, mu, Sigma, sde,
  tau = 2.5, ns = n)
}
\arguments{
\item{N}{Size of the full sample.}
\item{n}{Size of the Phase II subsample.}
\item{beta}{Regression coefficient vector for \code{(X, Z)} in the
outcome model. Its length determines the total number of covariates: the
first element corresponds to \code{X}, and the remaining elements correspond
to the components of \code{Z}.}
\item{mu}{Mean vector of \code{(X, Z)}. Its length must equal
\code{length(beta)}.}
\item{Sigma}{Variance-covariance matrix of \code{(X, Z)}. It must be a
positive-definite matrix with \code{length(beta)} rows and columns.}
\item{sde}{Standard deviation of the measurement error in \code{X_star}.}
\item{sdy}{Standard deviation of the error term in the linear outcome model.}
\item{ns}{Number of subjects selected into Phase II from each sampling
stratum. Its entries must sum to \code{n}.

For logistic and Cox outcomes, \code{ns} must be either the single value
\code{n} or a vector of length two. Setting \code{ns = n} gives simple random
sampling of \code{n} subjects from the full sample. For a vector of length
two, the first entry is the number sampled from subjects with \code{Y = 1}
for logistic data or \code{D = 1} for Cox data; the second is the number
sampled from subjects with \code{Y = 0} or \code{D = 0}, respectively.

For linear outcomes, setting \code{ns = n} gives simple random sampling.
With outcome-dependent sampling, \code{length(ns)} must equal
\code{length(quantile) + 1}, and the entries specify the numbers sampled from
the outcome strata in increasing order of \code{Y}. The default is \code{n}.}
\item{quantile}{For linear outcomes, a vector of strictly increasing
probabilities between 0 and 1 used to define the strata for outcome-dependent
sampling. If \code{quantile = c(q1, q2, ...)}, the strata are formed by
partitioning \code{Y} at its empirical \code{q1}, \code{q2}, ... quantiles,
in increasing order. If \code{quantile = NULL}, all subjects belong to a
single stratum and Phase II subjects are selected by simple random sampling.
The default is \code{NULL}.}
\item{tau}{Administrative study end time for the Cox model.}
}
\description{
Generate two-phase data under a linear, logistic, or Cox proportional hazards
outcome model. The expensive covariate \code{X} and the cheap covariates
\code{Z} are generated jointly from a multivariate normal distribution. A
surrogate measurement \code{X_star} is observed for every subject, whereas
\code{X} is observed only for subjects selected into the Phase II subsample.
Phase II subjects may be selected by simple random sampling or by
outcome-dependent stratified sampling.
}
\details{
Let
\deqn{W = (X, Z^{\mathsf T})^{\mathsf T}, \qquad W \sim N_p(\mu, \Sigma),}{
W = (X, Z')', W is distributed as N_p(mu, Sigma),}
where \eqn{p = \operatorname{length}(\beta)}{p = length(beta)}. If
\eqn{p = 1}, then
\eqn{W = X} and \code{Z} is absent. The surrogate is generated as
\deqn{X^* = X + \epsilon, \qquad \epsilon \sim N(0, \mathit{sde}^2),}{
X* = X + epsilon, epsilon is distributed as N(0, sde^2),}
independently of the other variables.

For linear data,
\deqn{Y = W^{\mathsf T}\beta + \epsilon_Y, \qquad
\epsilon_Y \sim N(0, \mathit{sdy}^2).}{
Y = W' beta + epsilon_Y, epsilon_Y is distributed as N(0, sdy^2).}

For logistic data,
\deqn{\Pr(Y = 1 \mid W) = \operatorname{expit}(W^{\mathsf T}\beta).}{
Pr(Y = 1 | W) = expit(W' beta).}

For Cox data, the hazard function is
\deqn{\lambda(t \mid W) = \frac{t}{2}\exp(W^{\mathsf T}\beta).}{
lambda(t | W) = (t / 2) exp(W' beta).}
Let \eqn{T} denote the failure time. The censoring time is generated as
\eqn{C = \min(C^*, \tau)}{C = min(C*, tau)}, where
\eqn{C^* \sim \operatorname{Uniform}(0, 5\tau/3)}{
C* is distributed as Uniform(0, 5 tau / 3)}. The observed follow-up time and
event indicator are
\deqn{Y = \min(T, C), \qquad D = I(T \le C).}{
Y = min(T, C), D = I(T <= C).}
}
\value{
A data frame with the following columns:
\itemize{
\item \code{R}: Phase II selection indicator, equal to 1 if \code{X} is
observed and 0 otherwise.
\item \code{Y}: Continuous or binary outcome for linear and logistic data,
or observed follow-up time for Cox data.
\item \code{D}: Event indicator for Cox data only.
\item \code{X}: Expensive covariate, set to \code{NA} when \code{R = 0}.
\item \code{X_star}: Surrogate measurement of \code{X}, observed for all
subjects.
\item \code{Z1}, ..., \code{Z(p-1)}: Cheap covariates, observed for all
subjects; absent when \code{length(beta) = 1}.
\item \code{stratum}: Sampling-stratum indicator. For linear data, strata
are numbered from 0 in increasing order of \code{Y}. For stratified logistic
and Cox data, it equals \code{Y} and \code{D}, respectively. Under simple
random sampling, it is 0 for every subject.
}
}
