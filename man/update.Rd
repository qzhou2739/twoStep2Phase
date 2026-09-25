\name{update_linear}
\alias{update_linear}
\alias{update_logistic}
\alias{update_cox}
\title{Update Estimation for Two-Phase Data}
\usage{
update_linear(data,
  z = grep("^Z[0-9]+$", names(data), value = TRUE), y = "Y",
  x = "X", x_star = "X_star", complete = "R", stratum = NULL,
  intercept = FALSE, x_models = NULL, h = NULL, initial = NULL,
  n_quad = 30L, nboot = 0, seed = NULL,
  var_method = "standard")

update_logistic(data,
  z = grep("^Z[0-9]+$", names(data), value = TRUE), y = "Y",
  x = "X", x_star = "X_star", complete = "R", stratum = NULL,
  intercept = FALSE, x_models = NULL, h = NULL, initial = NULL,
  tune_para = 0.01, nboot_tune = 50L,
  n_quad = 30L, nboot = 0, seed = NULL,
  var_method = "standard")

update_cox(data,
  z = grep("^Z[0-9]+$", names(data), value = TRUE), time = "Y",
  event = "D", x = "X", x_star = "X_star", complete = "R",
  stratum = NULL, x_models = NULL,
  h = NULL, initial = NULL, n_quad = 30L, nboot = 0, seed = NULL,
  var_method = "standard")
}
\arguments{
\item{data}{Data frame containing the full Phase I sample, the Phase II
selection indicator, and the observed Phase II values of the expensive
covariate.}
\item{z}{Character vector containing the names of the cheap-covariate
columns. By default, columns whose names consist of \code{"Z"} followed by
one or more digits, such as \code{"Z1"} and \code{"Z2"}, are selected
automatically. Use \code{NULL} or \code{character(0)} for an outcome model
without cheap covariates.}
\item{y}{Name of the outcome column for linear and logistic models. The
default is \code{"Y"}.}
\item{time}{Name of the column containing the observed event or censoring
time for the Cox model. The default is \code{"Y"}.}
\item{event}{Name of the event-indicator column for the Cox model, coded 1
for an event and 0 for censoring. The default is \code{"D"}.}
\item{x}{Name of the expensive-covariate column. The default is \code{"X"}.}
\item{x_star}{Name of the Phase I surrogate-measurement column for \code{x}.
The default is \code{"X_star"}. Use \code{NULL} when no surrogate measurement
is available. At least one of \code{x_star} and \code{z} must be supplied.}
\item{complete}{Name of the Phase II selection-indicator column, coded 1
when \code{x} is observed and 0 otherwise. The default is \code{"R"}.}
\item{stratum}{Name of the sampling-stratum column for outcome-dependent
Phase II sampling. The inverse-probability weight within each stratum is
estimated as the full-sample stratum size divided by its Phase II sample
size. If \code{NULL}, all subjects are treated as belonging to one stratum.
For outcome-dependent data produced by the package generators, specify
\code{stratum = "stratum"}.}
\item{intercept}{Logical value indicating whether an intercept is included in
the linear or logistic outcome model. The default is \code{FALSE}. When
\code{TRUE}, an intercept is included in the outcome model of interest
and in all working outcome models used to construct the update estimators.
This argument is not available for the Cox model.}
\item{x_models}{Working models for the conditional distribution of
\code{x} given the available Phase I auxiliary variables: \code{x_star}, when
supplied, together with the variables in \code{z}.

The available values are \code{"kernel"}, \code{"linear"}, and \code{"gam"},
with any combination allowed. The kernel option estimates the conditional
distribution nonparametrically. The linear option uses a linear Gaussian
working model, and the GAM option uses a Gaussian generalized additive model
with a separate additive smooth effect for each available auxiliary variable.

Use \code{NULL} to omit these optional conditional-\code{x} updates. The
complete-case analysis and default update are still computed. The default
is \code{"kernel"} without \code{z} and \code{"linear"} with \code{z}.}
\item{h}{Kernel bandwidth used when \code{x_models = "kernel"}. When one or
more \code{z} variables are present, the Gaussian product kernel uses a
separate bandwidth for each available auxiliary variable. Specifically, if
\eqn{s_j}{s_j} is the Phase II standard deviation of auxiliary variable
\eqn{j}{j}, its bandwidth is \eqn{h_j=s_j h}{h_j = s_j h}. If \code{h = NULL},
the common rule-of-thumb factor is
\eqn{h = \{4/(d+2)\}^{1/(d+4)} n^{-1/(d+4)}}{h = (4/(d+2))^(1/(d+4)) n^(-1/(d+4))},
where \eqn{d}{d} is the number of available auxiliary variables and \eqn{n}{n}
is the Phase II sample size. A supplied \code{h} is interpreted on this
standardized scale. Without \code{z}, the original one-dimensional bandwidth
calculation based on \code{x_star} is retained.}
\item{initial}{Numerical values used to initialize the solution of the working
estimating equations for the optional update methods.

For the original scalar kernel case (no intercept, no \code{z}, and no other
requested conditional-\code{x} model), supply two finite numbers in increasing
order. These give the initial search interval for the scalar coefficient. If
\code{initial = NULL}, the interval is formed by adding \code{c(-0.3, 0.3)} to
the complete-case estimate.

For every vector estimating equation, supply one finite starting value per
outcome coefficient. Order the values as the intercept, when requested,
followed by \code{x} and then the variables in \code{z}. Thus, the required
length is
\eqn{1 + \operatorname{length}(z) + I_{\mathrm{intercept}}}{
1 + length(z) + I(intercept)}. The same starting values are used for the Phase II
and full-sample equations under every requested kernel, linear, or GAM working
model. If \code{initial = NULL}, the complete-case coefficient estimates are
used.

\code{initial} is not used when the optional conditional-\code{x} updates are
omitted.}
\item{n_quad}{Number of Gauss-Hermite quadrature nodes used to approximate
integrals when \code{x_models} includes \code{"linear"} or \code{"gam"}.}
\item{tune_para}{Vector of non-negative candidate values for the
\eqn{L_2}{L2} penalty used to regularize the logistic working estimating
equations. If one value is supplied, it is used directly. If multiple values
are supplied, the bootstrap selection procedure described in the article is
used to select among them.}
\item{nboot_tune}{Number of bootstrap samples used to select among multiple
values of \code{tune_para}. It is ignored when only one value is supplied.}
\item{nboot}{Number of bootstrap samples used for standard-error or
covariance estimation. With \code{var_method = "standard"}, use 0 to omit
standard-error estimation or a whole number of at least 2. With
\code{var_method = "full_bootstrap"}, \code{nboot} must be greater than the
number of components in the largest working-estimator contrast
\eqn{\widetilde{\vartheta}_S-\widetilde{\vartheta}_F}{
vartheta_S - vartheta_F} used in any requested update.}
\item{seed}{Optional random seed set at the beginning of the analysis. It
affects logistic tuning-parameter selection, cross-validation for the ridge
update projection, and bootstrap calculations.}
\item{var_method}{Method used to estimate the covariance quantities required
for the update and its standard error. With \code{"standard"}, the update
projection is estimated from subject-level influence-function quantities,
using ridge regression where applicable, and the standard error is estimated
by repeating the complete update procedure on bootstrap samples.

With \code{"full_bootstrap"}, the covariance matrices
\eqn{\Sigma_{11}}, \eqn{\Sigma_{12}}, and \eqn{\Sigma_{22}} are estimated
jointly from bootstrap replicates of
\eqn{(\widehat{\theta}_S,
\widetilde{\vartheta}_S-\widetilde{\vartheta}_F)}{
(theta_S, vartheta_S - vartheta_F)}. These estimates are used for both the
update projection and the closed-form standard error.}
}
\description{
Fit complete-case and two-step update estimators for linear, logistic, or Cox
proportional hazards models using two-phase data. The expensive covariate
\code{x} is observed only for subjects selected into Phase II, whereas the
outcome and the available Phase I auxiliary variables are observed in the full
sample. These auxiliary variables consist of \code{x_star}, when supplied, and
the cheap covariates \code{z}.

Every analysis computes the complete-case estimator and the default update
estimator. Optional kernel, linear, or generalized additive working models
for the conditional distribution of \code{x} can be used to construct
additional update estimators.
}
\details{
The complete-case estimator, labelled \code{complete_case}, is computed from
the Phase II observations. When \code{stratum} is supplied, it is an
inverse-probability weighted estimator using the empirically estimated
stratum-specific sampling probabilities.

An update estimator has the form
\deqn{\widehat{\theta} =
\widehat{\theta}_S-\widehat{\Sigma}_{12}\widehat{\Sigma}_{22}^{-1}
(\widetilde{\vartheta}_S-\widetilde{\vartheta}_F),}{
theta_hat = theta_S - Sigma_12 Sigma_22^(-1)
(vartheta_S - vartheta_F),}
where \eqn{\widehat{\theta}_S}{theta_S} is the complete-case estimator and
\eqn{\widetilde{\vartheta}_S-\widetilde{\vartheta}_F}{
vartheta_S - vartheta_F} is a working-estimator contrast obtained from the
Phase II and full samples.

The \code{default} update uses the same model family as the outcome model and
the available Phase I auxiliary variables. When \code{x_star} is supplied, it
replaces \code{x}; when it is absent, the working outcome model uses \code{z}
alone. The
optional method labels have the following meanings:
\itemize{
\item \code{kernel}: kernel-based update using a nonparametric model for
\code{x} given the available auxiliary variables. With one or more \code{z}
variables, it uses the Euclidean distance between the auxiliary vectors after
each component has been divided by its Phase II standard deviation. The
auxiliary vector is \eqn{A=(X^*,Z)}{A = (X_star, Z)}, with unavailable
components omitted;
\item \code{joint_kernel}: joint update combining the default and kernel
working-estimator contrasts;
\item \code{linear}: update using a linear Gaussian working model for
\code{x} given the available auxiliary variables;
\item \code{joint_linear}: joint update combining the default and linear
working-estimator contrasts;
\item \code{gam}: update using a Gaussian generalized additive working model
for \code{x} given the available auxiliary variables; and
\item \code{joint_gam}: joint update combining the default and GAM
working-estimator contrasts.
}
For every requested conditional-\code{x} model, both its individual update and
the corresponding joint update are returned when the vector estimator is
used. The original scalar kernel path returns \code{kernel} without a separate
\code{joint_kernel} column.

}
\value{
An object containing:
\itemize{
\item \code{call}: the matched function call;
\item \code{family}: the outcome-model family;
\item \code{estimate}: point estimates for every coefficient and analysis
method;
\item \code{se}: corresponding standard-error estimates when requested;
\item \code{results}: a printable data frame containing the method,
coefficient, estimate, and standard error;
\item \code{details}: calculation details, including the empirical
inverse-probability weights, requested working models, sampling-stratum
specification, quadrature setting, and bootstrap information; and
\item \code{tune_parameters}: calculated or selected tuning quantities,
including the kernel bandwidth and logistic penalty parameter when applicable.
}

For the original no-intercept scalar analysis with no \code{z} and only the
kernel optional method, \code{estimate} is a named vector. Otherwise it is a
matrix whose rows correspond to the intercept, when requested, followed by
\code{x} and the variables in \code{z}; its columns correspond to the analysis
methods. When standard errors are not requested, they are displayed as
\code{NA} in \code{results}.
}
\examples{
library(twoStep2Phase)

##### Examples of optimal update estimation

## Linear model

N = 1000
n = 400
beta = 0.5
p = length(beta)
mu = rep(0,p)
sigma = diag(p)
sde = 0.5
ns = c(120,280)
sdy = 0.5
qq = 0.7

data = generate_linear_data(N, n, beta, mu, sigma, sde, sdy, ns, qq)

res = update_linear(data, z = NULL, y = "Y",
                    x = "X", x_star = "X_star", complete = "R", stratum = "stratum",
                    x_models = "kernel", h = NULL, initial = NULL,
                    n_quad = 30L, nboot = 0, seed = NULL,
                    var_method = "standard")

res$results



## Logistic model

N = 1000
n = 400
beta = 0.5
p = length(beta)
mu = rep(0,p)
sigma = diag(p)
sde = 0.5
ns = c(120,280)

data = generate_logistic_data(N, n, beta, mu, sigma, sde, ns)

res = update_logistic(data, z = NULL, y = "Y",
                      x = "X", x_star = "X_star", complete = "R", stratum = "stratum",
                      x_models = "kernel", h = NULL, initial = NULL,
                      tune_para = 0.01, nboot_tune = 50L,
                      n_quad = 30L, nboot = 0, seed = NULL,
                      var_method = "standard")

res$results



## Cox model

N = 1000
n = 400
beta = 0.5
p = length(beta)
mu = rep(0,p)
sigma = diag(p)
sde = 0.5
ns = c(120,280)
tau = 2.5

data = generate_cox_data(N, n, beta, mu, sigma, sde, tau, ns)

res = update_cox(data, z = NULL, time = "Y",
                 event = "D", x = "X", x_star = "X_star", complete = "R",
                 stratum = "stratum", x_models = "kernel", h = NULL, initial = NULL,
                 n_quad = 30L, nboot = 0, seed = NULL,
                 var_method = "standard")

res$results




##### Examples of joint update estimation

## Linear model

N = 1000
n = 400
beta = c(1,0.5)
p = length(beta)
mu = rep(0,p)
sigma = diag(p)
sde = 0.5
ns = c(120,280)
sdy = 0.5
qq = 0.7

data = generate_linear_data(N, n, beta, mu, sigma, sde, sdy, ns, qq)

res = update_linear(data, z = grep("^Z[0-9]+$", names(data), value = TRUE), y = "Y",
                    x = "X", x_star = "X_star", complete = "R", stratum = "stratum",
                    x_models = c("linear","gam"), h = NULL, initial = NULL,
                    n_quad = 30L, nboot = 0, seed = NULL,
                    var_method = "standard")

res$results



## Logistic model

N = 1000
n = 400
beta = c(1,0.5)
p = length(beta)
mu = rep(0,p)
sigma = diag(p)
sde = 0.5
ns = c(120,280)

data = generate_logistic_data(N, n, beta, mu, sigma, sde, ns)

res = update_logistic(data, z = grep("^Z[0-9]+$", names(data), value = TRUE), y = "Y",
                      x = "X", x_star = "X_star", complete = "R", stratum = "stratum",
                      x_models = c("linear","gam"), h = NULL, initial = NULL,
                      tune_para = 0.01, nboot_tune = 50L,
                      n_quad = 30L, nboot = 0, seed = NULL,
                      var_method = "standard")

res$results



## Cox model

N = 1000
n = 400
beta = c(1,0.5)
p = length(beta)
mu = rep(0,p)
sigma = diag(p)
sde = 0.5
ns = c(120,280)
tau = 2.5

data = generate_cox_data(N, n, beta, mu, sigma, sde, tau, ns)

res = update_cox(data, z = grep("^Z[0-9]+$", names(data), value = TRUE), time = "Y",
                 event = "D", x = "X", x_star = "X_star", complete = "R",
                 stratum = "stratum", x_models = c("linear","gam"), h = NULL, initial = NULL,
                 n_quad = 30L, nboot = 0, seed = NULL,
                 var_method = "standard")

res$results



##### Additional auxiliary-variable configurations

set.seed(7)
data_aux = generate_linear_data(
  200, 120, beta = c(0.5, 0.25), mu = c(0, 0),
  Sigma = diag(2), sde = 0.5
)

## Include an intercept and use a kernel based on (X_star, Z1)
update_linear(data_aux, intercept = TRUE, x_models = "kernel", nboot = 0)

## Fit a linear working model using X_star but no Z
update_linear(data_aux, z = NULL, x_models = "linear", n_quad = 5, nboot = 0)

## Analyze data without X_star by using Z1 as the Phase I auxiliary variable
data_aux$X_star = NULL
update_linear(
  data_aux, z = "Z1", x_star = NULL, x_models = "linear",
  n_quad = 5, nboot = 0
)
}
