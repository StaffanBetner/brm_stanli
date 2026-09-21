# `brm_stanli`

> **Status: experimental.** This is an unofficial compatibility layer between
> [`brms`](https://paulbuerkner.com/brms/) and
> [`stanli`](https://github.com/seantalts/stanli).
>
> It relies partly on non-exported `brms` internals and on Stanli's
> RStan-compatible adapter. Pin compatible package versions and run the test
> file before relying on it for an important analysis.

`brm_stanli()` uses `brms` to generate Stan code and Stan data, samples through
`stanli`, and converts the result into a `brmsfit`-compatible object.

Returned model objects inherit from both:

```r
class(fit)
# [1] "brm_stanli_fit" "brmsfit"
```

The `"brm_stanli_fit"` class provides Stanli-aware methods for LOO corrections
and model updates, while `"brmsfit"` retains ordinary `brms` post-processing.

## Installation

Install the required R packages:

```r
install.packages(c(
  "brms",
  "loo",
  "rstan",
  "testthat",
  "lme4"
))

install.packages(
  "stanli",
  repos = "https://seantalts.r-universe.dev"
)
```

Install the Stanli runtime once per machine/runtime version:

```r
stanli::stanli_install()
```

Then source the implementation:

```r
source("brm_stanli.R")
```

## Basic use

```r
library(brms)
library(stanli)

m0 <- brm_stanli(
  Reaction ~ 1 + Days + (1 + Days | Subject),
  data = lme4::sleepstudy,
  family = student(),
  chains = 4,
  cores = 4,
  seed = 1234
)

summary(m0)
pp_check(m0)
```

Most ordinary `brms` post-processing that depends on stored posterior draws and
model metadata should work, including:

```r
summary(m0)
plot(m0)
fixef(m0)
ranef(m0)
posterior_predict(m0)
posterior_epred(m0)
predict(m0)
fitted(m0)
log_lik(m0)
loo(m0)
waic(m0)
```

## PSIS-LOO cross-validation

### Ordinary PSIS-LOO

```r
m0_loo <- loo(m0)
```

### PSIS-LOO with moment matching

```r
m0_loo_mm <- loo(
  m0,
  moment_match = TRUE
)
```

Moment matching is applied to observations whose Pareto-$k$ values exceed
`k_threshold`, which defaults to `0.7`.

### Exact refits for problematic observations

```r
m0_loo_reloo <- loo(
  m0,
  reloo = TRUE
)
```

This computes ordinary PSIS-LOO first, identifies observations with problematic
Pareto-$k$ diagnostics, and exactly refits the model while omitting each such
observation.

### Combined correction workflow

Use both options to apply moment matching first and exact LOO refits only for
observations that remain problematic:

```r
m0_loo_corrected <- loo(
  m0,
  moment_match = TRUE,
  reloo = TRUE
)
```

The workflow is:

1. Compute ordinary PSIS-LOO.
2. Apply moment matching to problematic observations.
3. Identify any observations still above `k_threshold`.
4. Refit only those remaining observations exactly.

For Stanli fits, moment matching and exact LOO refits run sequentially with one
core during the correction stage. This avoids sending Stanli's live native model
handle to parallel workers.

You can pass options to each stage:

```r
m0_loo_corrected <- loo(
  m0,
  moment_match = TRUE,
  reloo = TRUE,
  moment_match_args = list(
    max_iters = 100L
  ),
  reloo_args = list(
    refit_args = list(
      chains = 2L,
      iter = 1000L,
      warmup = 500L,
      cores = 2L,
      seed = 20260921L
    )
  )
)
```

## Storing criteria with `add_criterion()`

Store a LOO result in the fit object:

```r
m0 <- add_criterion(
  m0,
  criterion = "loo",
  moment_match = TRUE
)

m0$criteria$loo
```

A subsequent plain call reuses the stored result:

```r
loo(m0)
```

To replace an already stored criterion with a new configuration, use
`overwrite = TRUE`:

```r
m0 <- add_criterion(
  m0,
  criterion = "loo",
  moment_match = TRUE,
  reloo = TRUE,
  overwrite = TRUE
)
```

The Stanli wrapper currently supports these criteria through `add_criterion()`:

```r
add_criterion(
  m0,
  criterion = c(
    "loo",
    "waic",
    "bayes_R2",
    "loo_R2"
  )
)
```

The following `brms` criteria are not implemented for `brm_stanli_fit`:

```r
add_criterion(m0, criterion = "kfold")
add_criterion(m0, criterion = "loo_subsample")
add_criterion(m0, criterion = "marglik")
```

## Updating models

Use `update()` normally for ordinary formula, sampling, prior, family, autocorrelation,
and data updates. The method refits through `brm_stanli()` and returns another
`"brm_stanli_fit"` object.

```r
m1 <- update(
  m0,
  formula. = ~ . + I(Days^2),
  iter = 4000L,
  warmup = 2000L,
  seed = 1234L
)
```

Update the fitting data with `newdata`:

```r
m2 <- update(
  m0,
  newdata = lme4::sleepstudy[1:120, ],
  iter = 2000L,
  warmup = 1000L,
  seed = 1234L
)
```

`update.brm_stanli_fit()` is intended for normal single-model workflows. It is
not yet complete feature parity with `brms::update.brmsfit()` for every advanced
`brms` configuration, especially complex multivariate, nonlinear, or `data2`
workflows.

## Caching fitted models

Use `file` to save a fitted object:

```r
m0 <- brm_stanli(
  Reaction ~ 1 + Days + (1 + Days | Subject),
  data = lme4::sleepstudy,
  family = student(),
  chains = 4,
  cores = 4,
  seed = 1234,
  file = "fits/sleepstudy_student",
  file_refit = "on_change"
)
```

When loading a cached fit, `brm_stanli()` rebuilds and reattaches Stanli's live
model automatically. This is necessary for moment matching, exact LOO refits,
and native density evaluation.

Supported `file_refit` values:

- `"never"`: load an existing cached fit when available;
- `"always"`: always refit and overwrite the cache;
- `"on_change"`: refit only when generated Stan code or Stan data changes.

## Running tests

Source the implementation and run the test file:

```r
source("brm_stanli.R")
testthat::test_file("test-brm-stanli.R")
```

Slow integration tests for moment matching and exact LOO are disabled by default.
Enable them explicitly:

```r
Sys.setenv(BRM_STANLI_RUN_SLOW_TESTS = "true")

source("brm_stanli.R")
testthat::test_file("test-brm-stanli.R")
```

The test file covers:

- creation of a live Stanli-backed `brmsfit`;
- posterior summaries and predictions;
- ordinary PSIS-LOO;
- stored LOO criteria;
- `update()`;
- cache restoration;
- moment matching;
- sequential moment matching plus exact reloo;
- explicit rejection of unsupported operations.

## Current limitations

This is experimental code, not an official `brms` or `stanli` integration.

### Unsupported operations

Do **not** use the following with `brm_stanli_fit` objects:

```r
brms::combine_models(...)
brms::kfold(...)
brms::bridge_sampler(...)
brms::bayes_factor(...)
brms::post_prob(...)
```

`combine_models_stanli()` exists only to fail explicitly:

```r
combine_models_stanli(m0)
```

The stock `brms::combine_models()` implementation must not be used with Stanli
fits because it combines underlying fits through RStan and loses the live Stanli
model adapter.

Bridge sampling, Bayes factors, and posterior model probabilities are blocked
because the stock `brms` bridge-sampling path rebuilds RStan internals and may
replace Stanli's live native model handle.

### Other current restrictions

- Exact `reloo()` does not support models using `data2`.
- Moment matching and exact LOO refits are sequential.
- Native RStan sampling functions must not be used on Stanli fits:

```r
rstan::sampling(m0$fit)
rstan::stan(fit = m0$fit)
```

- `recompile = TRUE` is not supported for Stanli moment matching or reloo.
- Exact LOO refits are potentially expensive because one model is refitted for
  each unresolved problematic observation.

## Background

The `loo` package uses Pareto-$k$ diagnostics to assess the reliability of
PSIS-LOO estimates. Moment matching can repair some problematic importance
sampling approximations. If problematic observations remain, exact leave-one-out
refits replace only their pointwise contributions. See the official
[`loo` moment-matching documentation](https://mc-stan.org/loo/reference/loo_moment_match.html)
and the [`brms` LOO documentation](https://paulbuerkner.com/brms/reference/loo.brmsfit.html).
