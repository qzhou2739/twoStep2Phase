# twoStep2Phase

`twoStep2Phase` implements the two-step update estimation approach described in
*An optimal two-step estimation approach for two-phase studies* by Zhou and
Wong. It supports linear, logistic, and Cox proportional hazards outcome
models.

In a two-phase study, the outcome and inexpensive variables are observed for
the full Phase I sample, while an expensive covariate, `X`, is measured only in
a Phase II subsample. The package starts with an inverse-probability-weighted
complete-case estimate and updates it using information from variables observed
in both phases. This can improve efficiency without requiring the working model
used in the update to be correctly specified.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("alexwky/twoStep2Phase")
```

## Main functions

| Outcome model | Analysis function | Data generator |
|---|---|---|
| Linear | `update_linear()` | `generate_linear_data()` |
| Logistic | `update_logistic()` | `generate_logistic_data()` |
| Cox proportional hazards | `update_cox()` | `generate_cox_data()` |

## Input data

The analysis functions accept a data frame containing:

- the outcome (`Y` by default), or follow-up time `Y` and event indicator `D`
  for a Cox model;
- the expensive covariate `X`, observed when `R = 1` and set to `NA` otherwise;
- an optional Phase I surrogate measurement `X_star`, observed for everyone;
- optional inexpensive covariates such as `Z1`, `Z2`, and so on; and
- the Phase II indicator `R`, plus an optional sampling-stratum variable.

Column names can be changed through the function arguments. Columns named
`Z1`, `Z2`, ... are detected automatically. For outcome-dependent Phase II
sampling, supply the stratum column using `stratum = "stratum"`. At least one
of `X_star` and the `Z` variables must be available; use `x_star = NULL` when
there is no surrogate measurement.

Linear and logistic outcome models can include an intercept by setting
`intercept = TRUE`. Each of `"kernel"`, `"linear"`, and `"gam"` can be used
with or without `Z`, and multiple working models can be requested together.
The kernel method uses all available auxiliary variables (`X_star` and `Z`).
When `Z` is present, it uses a multivariate Gaussian product kernel with a
rule-of-thumb bandwidth for each dimension.

## Example

The package includes generators for creating example two-phase data. The code
below simulates a linear outcome with one expensive covariate and one cheap
covariate, then fits the update estimators.

```r
library(twoStep2Phase)

set.seed(123)
dat <- generate_linear_data(
  N = 1000,
  n = 400,
  beta = c(1, 0.5),
  mu = c(0, 0),
  Sigma = diag(2),
  sde = 0.5,
  sdy = 0.5,
  ns = c(120, 280),
  quantile = 0.7
)

fit <- update_linear(
  data = dat,
  stratum = "stratum",
  x_models = c("linear", "gam"),
  nboot = 0,
  seed = 123
)

fit
fit$results
```

The same workflow applies to binary and time-to-event outcomes:

```r
## Logistic outcome
set.seed(123)
dat_logistic <- generate_logistic_data(
  N = 1000, n = 400,
  beta = c(1, 0.5), mu = c(0, 0), Sigma = diag(2), sde = 0.5,
  ns = c(120, 280)
)
fit_logistic <- update_logistic(
  dat_logistic,
  stratum = "stratum",
  x_models = "linear",
  tune_para = 0.01,
  nboot = 0
)

## Cox proportional hazards outcome
set.seed(123)
dat_cox <- generate_cox_data(
  N = 1000, n = 400,
  beta = c(1, 0.5), mu = c(0, 0), Sigma = diag(2), sde = 0.5,
  ns = c(120, 280)
)
fit_cox <- update_cox(
  dat_cox,
  stratum = "stratum",
  x_models = "linear",
  nboot = 0
)

fit_logistic$results
fit_cox$results
```

Additional configurations include an outcome-model intercept, a kernel with
both `X_star` and `Z`, linear or GAM working models without `Z`, and an analysis
without `X_star`:

```r
## Intercept and a multivariate kernel based on (X_star, Z1)
fit_intercept <- update_linear(
  dat, intercept = TRUE, x_models = "kernel", nboot = 0
)

## Linear working model based only on X_star
fit_no_z <- update_linear(
  dat, z = NULL, x_models = "linear", nboot = 0
)

## No X_star: use Z1 as the Phase I auxiliary variable
dat_no_xstar <- dat
dat_no_xstar$X_star <- NULL
fit_no_xstar <- update_linear(
  dat_no_xstar, x_star = NULL, z = "Z1",
  x_models = c("kernel", "linear"), nboot = 0
)
```

## Output

Each analysis returns a `twostep2phase` object. Its main components are:

- `estimate`: coefficient estimates for the complete-case and update methods;
- `se`: bootstrap standard errors when `nboot > 0`;
- `results`: a convenient table of estimates and standard errors;
- `details`: sampling weights and calculation details; and
- `tune_parameters`: selected bandwidth or penalty values, when applicable.

Every fit includes the `complete_case` and `default` estimates. Depending on
`x_models`, the output may also include `kernel`, `linear`, `gam`, and their
corresponding joint updates. For a multivariate kernel, `tune_parameters$h`
contains the common rule-of-thumb factor applied after each auxiliary variable
is scaled by its Phase II standard deviation. Set `nboot` to the desired number
of bootstrap replicates when standard errors are required.

For full argument descriptions, see `?update_linear`, `?update_logistic`,
`?update_cox`, and the corresponding data-generator help pages.
