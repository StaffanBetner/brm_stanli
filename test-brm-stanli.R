# Integration tests for brm_stanli.R
#
# Run manually after sourcing the implementation:
#
# source("brm_stanli.R")
# testthat::test_file("test-brm-stanli.R")
#
# The slow tests are disabled by default. Enable them with:
#
# Sys.setenv(BRM_STANLI_RUN_SLOW_TESTS = "true")
# testthat::test_file("test-brm-stanli.R")

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop(
    "Package 'testthat' must be installed to run this test file.",
    call. = FALSE
  )
}

if (!exists("brm_stanli", mode = "function", inherits = TRUE)) {
  stop(
    paste0(
      "Could not find brm_stanli(). Source brm_stanli.R before running ",
      "this test file."
    ),
    call. = FALSE
  )
}

brm_stanli_test_requirements <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("stanli")
  testthat::skip_if_not_installed("rstan")
  testthat::skip_if_not_installed("loo")
  testthat::skip_if_not_installed("lme4")

  testthat::skip_if(
    !stanli::stanli_available(),
    "Stanli is not available. Run stanli::stanli_install() before testing."
  )
}

brm_stanli_test_slow_enabled <- function() {
  identical(
    tolower(Sys.getenv("BRM_STANLI_RUN_SLOW_TESTS", unset = "false")),
    "true"
  )
}

brm_stanli_skip_if_slow_disabled <- function() {
  testthat::skip_if_not(
    brm_stanli_test_slow_enabled(),
    paste0(
      "Slow Stanli integration test skipped. Set ",
      "BRM_STANLI_RUN_SLOW_TESTS=true to run it."
    )
  )
}

brm_stanli_make_fast_fit <- function(seed) {
  brm_stanli_test_requirements()

  brm_stanli(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy,
    family = brms::gaussian(),
    chains = 2L,
    iter = 200L,
    warmup = 100L,
    cores = 2L,
    seed = seed,
    refresh = 0L
  )
}

brm_stanli_make_slow_fit <- function(seed) {
  brm_stanli_test_requirements()

  brm_stanli(
    Reaction ~ Days + (1 + Days | Subject),
    data = lme4::sleepstudy,
    family = brms::student(),
    chains = 4L,
    iter = 1000L,
    warmup = 500L,
    cores = 4L,
    seed = seed,
    refresh = 0L
  )
}

testthat::test_that("brm_stanli creates a live Stanli-backed brmsfit", {
  fit <- brm_stanli_make_fast_fit(seed = 1001L)

  testthat::expect_s3_class(fit, "brm_stanli_fit")
  testthat::expect_s3_class(fit, "brmsfit")
  testthat::expect_s4_class(fit$fit, "stanli_stanfit")

  testthat::expect_silent(
    brm_stanli_require_live_model(fit)
  )

  unconstrained_zero <- rep(
    0,
    rstan::get_num_upars(fit$fit)
  )

  log_density <- rstan::log_prob(
    fit$fit,
    upars = unconstrained_zero,
    adjust_transform = TRUE
  )

  testthat::expect_length(log_density, 1L)
  testthat::expect_true(is.finite(log_density))
})

testthat::test_that("ordinary posterior post-processing works", {
  fit <- brm_stanli_make_fast_fit(seed = 1002L)

  fixed_effects <- brms::fixef(fit)
  random_effects <- brms::ranef(fit)

  posterior_expectation <- brms::posterior_epred(
    fit,
    newdata = lme4::sleepstudy[1:5, , drop = FALSE]
  )

  posterior_prediction <- brms::posterior_predict(
    fit,
    newdata = lme4::sleepstudy[1:5, , drop = FALSE]
  )

  log_likelihood <- brms::log_lik(fit)

  testthat::expect_true(is.matrix(fixed_effects))
  testthat::expect_true(is.list(random_effects))

  testthat::expect_equal(
    dim(posterior_expectation)[2],
    5L
  )

  testthat::expect_equal(
    dim(posterior_prediction)[2],
    5L
  )

  testthat::expect_equal(
    dim(log_likelihood)[2],
    nrow(lme4::sleepstudy)
  )
})

testthat::test_that("ordinary PSIS-LOO uses the Stanli-aware method", {
  fit <- brm_stanli_make_fast_fit(seed = 1003L)

  loo_result <- suppressWarnings(
    brms::loo(
      fit,
      cores = 1L
    )
  )

  testthat::expect_s3_class(loo_result, "loo")

  testthat::expect_identical(
    attr(loo_result, "model_name"),
    "fit"
  )

  testthat::expect_equal(
    nrow(loo_result$pointwise),
    nrow(lme4::sleepstudy)
  )

  testthat::expect_true(
    all(
      c("elpd_loo", "p_loo", "looic") %in%
        colnames(loo_result$pointwise)
    )
  )
})

testthat::test_that("add_criterion stores and reuses ordinary LOO", {
  fit <- brm_stanli_make_fast_fit(seed = 1004L)

  fit_with_loo <- suppressWarnings(
    brms::add_criterion(
      fit,
      criterion = "loo"
    )
  )

  stored_loo <- fit_with_loo$criteria$loo

  reused_loo <- suppressWarnings(
    brms::loo(fit_with_loo)
  )

  testthat::expect_s3_class(fit_with_loo, "brm_stanli_fit")
  testthat::expect_s3_class(stored_loo, "loo")
  testthat::expect_s3_class(reused_loo, "loo")

  testthat::expect_equal(
    stored_loo$estimates,
    reused_loo$estimates,
    tolerance = 1e-12
  )

  testthat::expect_equal(
    stored_loo$pointwise,
    reused_loo$pointwise,
    tolerance = 1e-12
  )
})

testthat::test_that("add_criterion supports WAIC and Bayesian R-squared", {
  fit <- brm_stanli_make_fast_fit(seed = 1005L)

  fit_with_criteria <- suppressWarnings(
    brms::add_criterion(
      fit,
      criterion = c("waic", "bayes_R2")
    )
  )

  testthat::expect_s3_class(
    fit_with_criteria,
    "brm_stanli_fit"
  )

  testthat::expect_s3_class(
    fit_with_criteria$criteria$waic,
    "waic"
  )

  testthat::expect_true(
    is.matrix(fit_with_criteria$criteria$bayes_R2)
  )

  testthat::expect_true(
    "R2" %in% colnames(fit_with_criteria$criteria$bayes_R2)
  )
})

testthat::test_that("add_criterion rejects unsupported criteria", {
  fit <- brm_stanli_make_fast_fit(seed = 1006L)

  testthat::expect_error(
    brms::add_criterion(
      fit,
      criterion = "kfold"
    ),
    "Unsupported criterion"
  )

  testthat::expect_error(
    brms::add_criterion(
      fit,
      criterion = "loo_subsample"
    ),
    "Unsupported criterion"
  )

  testthat::expect_error(
    brms::add_criterion(
      fit,
      criterion = "marglik"
    ),
    "Unsupported criterion"
  )
})

testthat::test_that("update preserves the Stanli class and refits through brm_stanli", {
  fit <- brm_stanli_make_fast_fit(seed = 1007L)

  updated_fit <- update(
    fit,
    formula. = ~ . + I(Days^2),
    iter = 240L,
    warmup = 120L,
    seed = 1008L,
    refresh = 0L
  )

  fixed_effect_names <- rownames(
    brms::fixef(updated_fit)
  )

  testthat::expect_s3_class(
    updated_fit,
    "brm_stanli_fit"
  )

  testthat::expect_s3_class(
    updated_fit,
    "brmsfit"
  )

  testthat::expect_s4_class(
    updated_fit$fit,
    "stanli_stanfit"
  )

  testthat::expect_silent(
    brm_stanli_require_live_model(updated_fit)
  )

  testthat::expect_true(
    any(grepl("I\\(Days\\^2\\)", fixed_effect_names))
  )

  testthat::expect_equal(
    updated_fit$fit@sim$iter,
    240L
  )

  testthat::expect_equal(
    updated_fit$fit@sim$warmup,
    120L
  )
})

testthat::test_that("update supports replacing the fitting data", {
  fit <- brm_stanli_make_fast_fit(seed = 1009L)

  smaller_data <- lme4::sleepstudy[
    1:120,
    ,
    drop = FALSE
  ]

  updated_fit <- update(
    fit,
    newdata = smaller_data,
    iter = 240L,
    warmup = 120L,
    seed = 1010L,
    refresh = 0L
  )

  updated_log_likelihood <- brms::log_lik(updated_fit)

  testthat::expect_s3_class(
    updated_fit,
    "brm_stanli_fit"
  )

  testthat::expect_equal(
    nrow(updated_fit$data),
    nrow(smaller_data)
  )

  testthat::expect_equal(
    dim(updated_log_likelihood)[2],
    nrow(smaller_data)
  )
})

testthat::test_that("cached fits restore a live Stanli model", {
  brm_stanli_test_requirements()

  cache_file <- tempfile(
    pattern = "brm_stanli_cache_",
    fileext = ".rds"
  )

  on.exit(
    unlink(cache_file),
    add = TRUE
  )

  fitted_model <- brm_stanli(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy,
    family = brms::gaussian(),
    chains = 2L,
    iter = 200L,
    warmup = 100L,
    cores = 2L,
    seed = 1011L,
    refresh = 0L,
    file = cache_file,
    file_refit = "always"
  )

  restored_model <- brm_stanli(
    Reaction ~ Days + (1 | Subject),
    data = lme4::sleepstudy,
    family = brms::gaussian(),
    chains = 2L,
    iter = 200L,
    warmup = 100L,
    cores = 2L,
    seed = 1011L,
    refresh = 0L,
    file = cache_file,
    file_refit = "never"
  )

  unconstrained_zero <- rep(
    0,
    rstan::get_num_upars(restored_model$fit)
  )

  restored_log_density <- rstan::log_prob(
    restored_model$fit,
    upars = unconstrained_zero,
    adjust_transform = TRUE
  )

  testthat::expect_s3_class(
    fitted_model,
    "brm_stanli_fit"
  )

  testthat::expect_s3_class(
    restored_model,
    "brm_stanli_fit"
  )

  testthat::expect_s4_class(
    restored_model$fit,
    "stanli_stanfit"
  )

  testthat::expect_true(
    is.finite(restored_log_density)
  )

  testthat::expect_silent(
    brm_stanli_require_live_model(restored_model)
  )
})

testthat::test_that("unsupported Stanli operations fail explicitly", {
  fit <- brm_stanli_make_fast_fit(seed = 1012L)

  testthat::expect_error(
    combine_models_stanli(fit),
    "not implemented yet"
  )

  testthat::expect_error(
    brms::kfold(fit),
    "not implemented for brm_stanli_fit"
  )

  testthat::expect_error(
    brms::bridge_sampler(fit),
    "not implemented for brm_stanli_fit"
  )

  testthat::expect_error(
    brms::bayes_factor(fit, fit),
    "not implemented for brm_stanli_fit"
  )

  testthat::expect_error(
    brms::post_prob(fit),
    "not implemented for brm_stanli_fit"
  )
})

testthat::test_that("moment matching corrects or improves Pareto-k diagnostics", {
  brm_stanli_skip_if_slow_disabled()

  fit <- brm_stanli_make_slow_fit(seed = 1013L)

  ordinary_loo <- suppressWarnings(
    brms::loo(
      fit,
      cores = 1L
    )
  )

  problematic_rows <- loo::pareto_k_ids(
    ordinary_loo,
    threshold = 0.7
  )

  testthat::skip_if(
    length(problematic_rows) == 0L,
    "This stochastic fitted model has no problematic Pareto-k values."
  )

  corrected_loo <- brms::loo(
    fit,
    moment_match = TRUE,
    moment_match_args = list(
      max_iters = 100L
    )
  )

  corrected_rows <- loo::pareto_k_ids(
    corrected_loo,
    threshold = 0.7
  )

  testthat::expect_s3_class(
    corrected_loo,
    "loo"
  )

  testthat::expect_lte(
    length(corrected_rows),
    length(problematic_rows)
  )
})

testthat::test_that("moment matching and reloo run sequentially", {
  brm_stanli_skip_if_slow_disabled()

  fit <- brm_stanli_make_slow_fit(seed = 1014L)

  corrected_loo <- brms::loo(
    fit,
    moment_match = TRUE,
    reloo = TRUE,
    moment_match_args = list(
      max_iters = 100L
    ),
    reloo_args = list(
      refit_args = list(
        chains = 2L,
        iter = 500L,
        warmup = 250L,
        cores = 2L,
        seed = 1015L
      )
    )
  )

  remaining_problematic_rows <- loo::pareto_k_ids(
    corrected_loo,
    threshold = 0.7
  )

  testthat::expect_s3_class(
    corrected_loo,
    "loo"
  )

  testthat::expect_length(
    remaining_problematic_rows,
    0L
  )

  testthat::expect_true(
    !is.null(attr(corrected_loo, "reloo_stanli"))
  )
})
