# Experimental brms + stanli integration.
#
# This is an unofficial compatibility layer. It uses brms to generate
# Stan code/data and stanli to sample, then exposes a brmsfit-compatible
# object with class c("brm_stanli_fit", "brmsfit").
#
# Source this file before fitting:
#
# source("brm_stanli.R")
#
# Required packages:
#
# install.packages(c("brms", "loo", "rstan"))
# install.packages("stanli", repos = "https://seantalts.r-universe.dev")
# stanli::stanli_install()
#
# See README for supported workflows and current limitations.

brm_stanli_require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' must be installed.", package),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

brm_stanli_validate_integer <- function(x, name, lower = NULL) {
  if (length(x) != 1L ||
      is.na(x) ||
      !is.numeric(x) ||
      x != as.integer(x)) {
    stop(sprintf("'%s' must be an integer.", name), call. = FALSE)
  }

  x <- as.integer(x)

  if (!is.null(lower) && x < lower) {
    stop(
      sprintf("'%s' must be greater than or equal to %d.", name, lower),
      call. = FALSE
    )
  }

  x
}

brm_stanli_validate_named_list <- function(x, name) {
  if (!is.list(x)) {
    stop(sprintf("'%s' must be a list.", name), call. = FALSE)
  }

  if (length(x) > 0L &&
      (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop(sprintf("'%s' must be a named list.", name), call. = FALSE)
  }

  x
}

brm_stanli_validate_cache_path <- function(file) {
  if (is.null(file)) {
    return(NULL)
  }

  if (!is.character(file) ||
      length(file) != 1L ||
      is.na(file) ||
      !nzchar(file)) {
    stop("'file' must be NULL or a non-empty character string.", call. = FALSE)
  }

  if (!identical(tolower(tools::file_ext(file)), "rds")) {
    file <- paste0(file, ".rds")
  }

  directory <- dirname(file)

  if (!dir.exists(directory)) {
    stop(
      sprintf(
        "Directory '%s' does not exist. Create it before saving the model.",
        directory
      ),
      call. = FALSE
    )
  }

  file
}

brm_stanli_threads_per_chain <- function(threads) {
  if (is.null(threads)) {
    return(1L)
  }

  if (is.numeric(threads) && length(threads) == 1L) {
    return(
      brm_stanli_validate_integer(
        threads,
        "threads",
        lower = 1L
      )
    )
  }

  if (inherits(threads, "brmsthreads")) {
    threads_per_chain <- threads[["threads"]]

    if (is.null(threads_per_chain)) {
      stop(
        "Could not read the thread count from the 'threads' object.",
        call. = FALSE
      )
    }

    return(
      brm_stanli_validate_integer(
        threads_per_chain,
        "threads$threads",
        lower = 1L
      )
    )
  }

  stop(
    paste0(
      "'threads' must be NULL, a positive integer, or ",
      "a brms::threading() object."
    ),
    call. = FALSE
  )
}

brm_stanli_mark_fit <- function(x) {
  if (!inherits(x, "brmsfit")) {
    stop("'x' must be a brmsfit object.", call. = FALSE)
  }

  class(x) <- unique(c("brm_stanli_fit", class(x)))

  x
}

brm_stanli_unmark_fit <- function(x) {
  class(x) <- setdiff(class(x), "brm_stanli_fit")

  x
}

brm_stanli_require_live_model <- function(x) {
  if (!inherits(x, "brm_stanli_fit")) {
    stop("'x' must inherit from 'brm_stanli_fit'.", call. = FALSE)
  }

  if (!inherits(x$fit, "stanli_stanfit")) {
    stop(
      "The fit does not contain a Stanli-compatible stanfit object.",
      call. = FALSE
    )
  }

  if (!inherits(x$fit@.MISC$stanli_model, "stanli_model")) {
    stop(
      paste0(
        "This fit has no live Stanli model. If it was loaded from disk, ",
        "reload it through brm_stanli(..., file_refit = 'never') so that ",
        "the native Stanli handle is restored."
      ),
      call. = FALSE
    )
  }

  invisible(x)
}

brm_stanli_set_data_name <- function(x, data_name) {
  attr(x$data, "data_name") <- data_name

  x
}

brm_stanli_data_name <- function(data, fallback = "data") {
  data_name <- attr(data, "data_name", exact = TRUE)

  if (is.null(data_name) ||
      !is.character(data_name) ||
      length(data_name) != 1L ||
      is.na(data_name) ||
      !nzchar(data_name)) {
    return(fallback)
  }

  data_name
}

brm_stanli_model_signature <- function(x) {
  list(
    code = brms::make_stancode(x),
    data = brms::make_standata(x),
    algorithm = "sampling"
  )
}

brm_stanli_log_mean_exp <- function(x) {
  x <- as.numeric(x)

  if (length(x) == 0L) {
    stop("'x' must contain at least one value.", call. = FALSE)
  }

  max_x <- max(x)

  if (is.infinite(max_x) && max_x < 0) {
    return(-Inf)
  }

  max_x + log(mean(exp(x - max_x)))
}

brm_stanli_restore_model <- function(x) {
  if (!inherits(x, "brmsfit")) {
    stop(
      "The object read from the cache file is not a brmsfit object.",
      call. = FALSE
    )
  }

  if (!inherits(x$fit, "stanli_stanfit")) {
    return(brm_stanli_mark_fit(x))
  }

  stanli_metadata <- attr(x, "stanli_metadata")

  if (is.null(stanli_metadata) ||
      is.null(stanli_metadata$seed) ||
      is.null(stanli_metadata$threads_per_chain)) {
    stop(
      paste0(
        "The cached fit does not contain Stanli metadata required to restore ",
        "a live model. Refit with file_refit = 'always'."
      ),
      call. = FALSE
    )
  }

  restored_model <- stanli::stanli_model(
    code = brms::make_stancode(x),
    data = brms::make_standata(x),
    seed = stanli_metadata$seed,
    threads_per_chain = stanli_metadata$threads_per_chain
  )

  saved_fnames_oi <- x$fit@sim$fnames_oi
  if (!is.null(x$fit@sim$fnames_oi_old)) {
    x$fit@sim$fnames_oi <- x$fit@sim$fnames_oi_old
  }

  x$fit <- stanli::as_stanfit(
    x$fit,
    model = restored_model
  )

  x$fit@sim$fnames_oi <- saved_fnames_oi

  brm_stanli_mark_fit(x)
}

brm_stanli_read_cached_fit <- function(path, restore_model = TRUE) {
  if (!file.exists(path)) {
    return(NULL)
  }

  cached_fit <- tryCatch(
    readRDS(path),
    error = function(error) {
      stop(
        sprintf(
          "Could not read cached fit from '%s': %s",
          path,
          conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  if (!inherits(cached_fit, "brmsfit")) {
    stop(
      sprintf(
        "The object stored in '%s' is not a brmsfit object.",
        path
      ),
      call. = FALSE
    )
  }

  cached_fit$file <- path

  if (restore_model) {
    return(brm_stanli_restore_model(cached_fit))
  }

  brm_stanli_mark_fit(cached_fit)
}

brm_stanli_write_cached_fit <- function(x, path, compress) {
  x$file <- path

  saveRDS(
    x,
    file = path,
    compress = compress
  )

  x
}

brm_stanli <- function(
  formula,
  data,
  ...,
  chains = 4L,
  iter = getOption("brms.iter", 2000L),
  warmup = floor(iter / 2),
  thin = 1L,
  cores = getOption("mc.cores", 1L),
  threads = getOption("brms.threads", NULL),
  control = NULL,
  seed = NA_integer_,
  init = NULL,
  init_radius = 2,
  stanli_init = NULL,
  pathfinder_init = NULL,
  save_warmup = FALSE,
  refresh = 100L,
  file = NULL,
  file_compress = TRUE,
  file_refit = getOption("brms.file_refit", "never"),
  empty = TRUE
) {
  brm_stanli_require_namespace("brms")
  brm_stanli_require_namespace("stanli")
  brm_stanli_require_namespace("rstan")

  data_name <- paste(
    deparse(substitute(data), width.cutoff = 500L),
    collapse = ""
  )

  if (!stanli::stanli_available()) {
    stop(
      paste0(
        "Stanli is not available. Run stanli::stanli_install() once, ",
        "then restart R."
      ),
      call. = FALSE
    )
  }

  if (!isTRUE(empty)) {
    stop(
      paste0(
        "brm_stanli() always creates an empty brmsfit object before ",
        "sampling with Stanli, so 'empty' must be TRUE."
      ),
      call. = FALSE
    )
  }

  chains <- brm_stanli_validate_integer(chains, "chains", lower = 1L)
  iter <- brm_stanli_validate_integer(iter, "iter", lower = 1L)
  warmup <- brm_stanli_validate_integer(warmup, "warmup", lower = 0L)
  thin <- brm_stanli_validate_integer(thin, "thin", lower = 1L)
  cores <- brm_stanli_validate_integer(cores, "cores", lower = 1L)
  refresh <- brm_stanli_validate_integer(refresh, "refresh", lower = 0L)

  if (warmup >= iter) {
    stop(
      paste0(
        "'warmup' must be smaller than 'iter', because at least one ",
        "post-warmup draw is required."
      ),
      call. = FALSE
    )
  }

  if (length(init_radius) != 1L ||
      is.na(init_radius) ||
      !is.numeric(init_radius) ||
      init_radius < 0) {
    stop(
      "'init_radius' must be a numeric value greater than or equal to 0.",
      call. = FALSE
    )
  }

  if (is.null(control)) {
    control <- list()
  }

  control <- brm_stanli_validate_named_list(control, "control")

  supported_control <- c(
    "adapt_delta",
    "max_treedepth"
  )

  unsupported_control <- setdiff(
    names(control),
    supported_control
  )

  if (length(unsupported_control) > 0L) {
    warning(
      paste0(
        "The following control arguments are ignored because ",
        "stanli::sample_model() does not expose them: ",
        paste(unsupported_control, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  delta <- if (!is.null(control$adapt_delta)) {
    control$adapt_delta
  } else {
    0.8
  }

  max_depth <- if (!is.null(control$max_treedepth)) {
    control$max_treedepth
  } else {
    10L
  }

  if (length(delta) != 1L ||
      is.na(delta) ||
      !is.numeric(delta) ||
      delta <= 0 ||
      delta >= 1) {
    stop(
      "'control$adapt_delta' must be a numeric scalar strictly between 0 and 1.",
      call. = FALSE
    )
  }

  max_depth <- brm_stanli_validate_integer(
    max_depth,
    "control$max_treedepth",
    lower = 1L
  )

  if (length(seed) != 1L || !is.numeric(seed)) {
    stop("'seed' must be NA or an integer.", call. = FALSE)
  }

  stanli_seed <- if (is.na(seed)) {
    sample.int(.Machine$integer.max, size = 1L)
  } else {
    brm_stanli_validate_integer(seed, "seed", lower = 1L)
  }

  threads_per_chain <- brm_stanli_threads_per_chain(threads)
  parallel_chains <- min(chains, cores)
  samples <- iter - warmup

  if (!is.null(stanli_init) && !is.null(init)) {
    stop(
      "Specify only one of 'init' and 'stanli_init'.",
      call. = FALSE
    )
  }

  sampler_init <- NULL

  if (!is.null(stanli_init)) {
    sampler_init <- stanli_init
  } else if (is.null(init) || identical(init, "random")) {
    sampler_init <- NULL
  } else if (is.numeric(init) && length(init) == 1L && init == 0) {
    sampler_init <- NULL
    init_radius <- 0
  } else {
    stop(
      paste0(
        "Lists and functions supplied through 'init' are not supported because ",
        "brm() uses constrained initial values while stanli::sample_model() ",
        "requires unconstrained vectors or matrices. Use 'stanli_init' instead."
      ),
      call. = FALSE
    )
  }

  if (!is.null(pathfinder_init) && !is.null(sampler_init)) {
    stop(
      "'pathfinder_init' cannot be combined with 'stanli_init'.",
      call. = FALSE
    )
  }

  file_refit <- match.arg(
    file_refit,
    choices = c(
      "never",
      "always",
      "on_change"
    )
  )

  cache_file <- brm_stanli_validate_cache_path(file)

  if (!is.null(cache_file) && identical(file_refit, "never")) {
    cached_fit <- brm_stanli_read_cached_fit(
      cache_file,
      restore_model = TRUE
    )

    if (!is.null(cached_fit)) {
      return(cached_fit)
    }
  }

  dots <- list(...)

  prohibited_dots <- c(
    "fit",
    "init",
    "inits",
    "chains",
    "iter",
    "warmup",
    "thin",
    "cores",
    "threads",
    "control",
    "algorithm",
    "backend",
    "future",
    "seed",
    "empty",
    "rename",
    "opencl",
    "file",
    "file_compress",
    "file_refit"
  )

  duplicated_arguments <- intersect(
    names(dots),
    prohibited_dots
  )

  if (length(duplicated_arguments) > 0L) {
    stop(
      paste0(
        "Do not supply these arguments through '...': ",
        paste(duplicated_arguments, collapse = ", "),
        ". Pass them directly to brm_stanli() instead."
      ),
      call. = FALSE
    )
  }

  brm_arguments <- c(
    list(
      formula = formula,
      data = data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      cores = cores,
      threads = threads,
      control = control,
      seed = stanli_seed,
      algorithm = "sampling",
      backend = "rstan",
      empty = TRUE,
      rename = FALSE
    ),
    dots
  )

  stanli_call_arguments <- c(
    list(
      formula = formula,
      data = data,
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      cores = cores,
      threads = threads,
      control = control,
      seed = stanli_seed,
      init = init,
      init_radius = init_radius,
      stanli_init = stanli_init,
      pathfinder_init = pathfinder_init,
      save_warmup = save_warmup,
      refresh = refresh,
      file = NULL,
      file_compress = file_compress,
      file_refit = "always",
      empty = TRUE
    ),
    dots
  )

  brms_fit <- do.call(
    brms::brm,
    brm_arguments
  )

  brms_fit <- brm_stanli_set_data_name(
    brms_fit,
    data_name
  )

  current_signature <- brm_stanli_model_signature(brms_fit)

  if (!is.null(cache_file) && identical(file_refit, "on_change")) {
    cached_fit <- brm_stanli_read_cached_fit(
      cache_file,
      restore_model = FALSE
    )

    if (!is.null(cached_fit)) {
      cached_metadata <- attr(cached_fit, "stanli_metadata")
      cached_signature <- cached_metadata$model_signature

      if (!is.null(cached_signature) &&
          identical(cached_signature, current_signature)) {
        return(brm_stanli_restore_model(cached_fit))
      }
    }
  }

  stanli_model <- stanli::stanli_model(
    code = current_signature$code,
    data = current_signature$data,
    seed = stanli_seed,
    threads_per_chain = threads_per_chain
  )

  stanli_fit <- stanli::sample_model(
    model = stanli_model,
    chains = chains,
    seed = stanli_seed,
    warmup = warmup,
    samples = samples,
    thin = thin,
    delta = delta,
    max_depth = max_depth,
    save_warmup = save_warmup,
    init = sampler_init,
    init_radius = init_radius,
    pathfinder_init = pathfinder_init,
    parallel_chains = parallel_chains,
    refresh = refresh,
    threads_per_chain = threads_per_chain
  )

  brms_fit$fit <- stanli::as_stanfit(stanli_fit)
  brms_fit <- brms::rename_pars(brms_fit)
  brms_fit <- brm_stanli_mark_fit(brms_fit)

  attr(brms_fit, "stanli_metadata") <- list(
    seed = stanli_seed,
    threads_per_chain = threads_per_chain,
    report = stanli_fit$report,
    model_signature = current_signature
  )

  attr(brms_fit, "stanli_call_arguments") <- stanli_call_arguments

  if (!is.null(cache_file)) {
    brms_fit <- brm_stanli_write_cached_fit(
      brms_fit,
      path = cache_file,
      compress = file_compress
    )
  }

  brms_fit
}

loo_moment_match_stanli <- function(
  x,
  loo = NULL,
  k_threshold = 0.7,
  newdata = NULL,
  resp = NULL,
  check = TRUE,
  max_iters = 30L,
  split = TRUE,
  cov = TRUE,
  cores = 1L
) {
  brm_stanli_require_namespace("brms")
  brm_stanli_require_namespace("loo")
  brm_stanli_require_namespace("rstan")
  brm_stanli_require_namespace("stanli")

  brm_stanli_require_live_model(x)

  if (length(cores) != 1L ||
      is.na(cores) ||
      !is.numeric(cores) ||
      cores != as.integer(cores) ||
      as.integer(cores) != 1L) {
    stop(
      "'cores' must be 1 for Stanli moment matching.",
      call. = FALSE
    )
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("'check' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(split) || length(split) != 1L || is.na(split)) {
    stop("'split' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(cov) || length(cov) != 1L || is.na(cov)) {
    stop("'cov' must be TRUE or FALSE.", call. = FALSE)
  }

  if (length(k_threshold) != 1L ||
      is.na(k_threshold) ||
      !is.numeric(k_threshold) ||
      k_threshold <= 0) {
    stop("'k_threshold' must be a positive numeric scalar.", call. = FALSE)
  }

  if (length(max_iters) != 1L ||
      is.na(max_iters) ||
      !is.numeric(max_iters) ||
      max_iters < 1L ||
      max_iters != as.integer(max_iters)) {
    stop("'max_iters' must be a positive integer.", call. = FALSE)
  }

  if (is.null(loo)) {
    loo <- brms::loo(
      x,
      resp = resp,
      cores = 1L
    )
  }

  if (!inherits(loo, "loo")) {
    stop("'loo' must inherit from class 'loo'.", call. = FALSE)
  }

  if (is.null(newdata)) {
    newdata <- stats::model.frame(x)
  } else {
    newdata <- as.data.frame(newdata)
  }

  if (check) {
    yhash_loo <- attr(loo, "yhash")
    yhash_fit <- brms:::hash_response(
      x,
      newdata = newdata,
      resp = resp
    )

    if (!brms:::is_equal(yhash_loo, yhash_fit)) {
      stop(
        paste0(
          "Response values used in 'loo' and 'x' do not match. ",
          "If this is a false positive, use check = FALSE."
        ),
        call. = FALSE
      )
    }
  }

  required_old_order_fields <- c(
    "pars_oi_old",
    "dims_oi_old",
    "fnames_oi_old"
  )

  missing_old_order_fields <- required_old_order_fields[
    !vapply(
      required_old_order_fields,
      function(name) !is.null(x$fit@sim[[name]]),
      logical(1)
    )
  ]

  if (length(missing_old_order_fields) > 0L) {
    stop(
      paste0(
        "The fit lacks original Stan parameter-order metadata required for ",
        "moment matching: ",
        paste(missing_old_order_fields, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  log_lik_i <- function(x, i, newdata, ...) {
    as.vector(
      brms::log_lik(
        x,
        newdata = newdata[i, , drop = FALSE],
        ...
      )
    )
  }

  unconstrain_pars <- function(x, pars, ...) {
    brms:::unconstrain_pars_stanfit(
      x$fit,
      pars = pars,
      ...
    )
  }

  log_prob_upars <- function(x, upars, ...) {
    brms:::log_prob_upars_stanfit(
      x$fit,
      upars = upars,
      ...
    )
  }

  log_lik_i_upars <- function(
    x,
    upars,
    i,
    ndraws = NULL,
    draw_ids = NULL,
    ...
  ) {
    x <- brms:::.update_pars(
      x,
      upars = upars,
      ...
    )

    log_lik_i(
      x,
      i = i,
      newdata = newdata,
      resp = resp
    )
  }

  out <- try(
    loo::loo_moment_match.default(
      x = x,
      loo = loo,
      post_draws = as.matrix,
      log_lik_i = log_lik_i,
      unconstrain_pars = unconstrain_pars,
      log_prob_upars = log_prob_upars,
      log_lik_i_upars = log_lik_i_upars,
      max_iters = as.integer(max_iters),
      k_threshold = k_threshold,
      split = split,
      cov = cov,
      cores = 1L,
      newdata = newdata,
      resp = resp
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    condition <- attr(out, "condition")

    underlying_message <- if (inherits(condition, "condition")) {
      conditionMessage(condition)
    } else {
      as.character(out)
    }

    stop(
      paste0(
        "Stanli moment matching failed.\n\n",
        "Underlying error:\n",
        underlying_message
      ),
      call. = FALSE
    )
  }

  out
}

reloo_stanli <- function(
  x,
  loo = NULL,
  k_threshold = 0.7,
  resp = NULL,
  check = TRUE,
  refit_args = list(),
  log_lik_args = list()
) {
  brm_stanli_require_namespace("brms")
  brm_stanli_require_namespace("loo")
  brm_stanli_require_namespace("stanli")

  brm_stanli_require_live_model(x)

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("'check' must be TRUE or FALSE.", call. = FALSE)
  }

  if (length(k_threshold) != 1L ||
      is.na(k_threshold) ||
      !is.numeric(k_threshold) ||
      k_threshold <= 0) {
    stop("'k_threshold' must be a positive numeric scalar.", call. = FALSE)
  }

  refit_args <- brm_stanli_validate_named_list(
    refit_args,
    "refit_args"
  )

  log_lik_args <- brm_stanli_validate_named_list(
    log_lik_args,
    "log_lik_args"
  )

  prohibited_refit_args <- c(
    "formula",
    "data",
    "file",
    "file_refit",
    "file_compress",
    "empty",
    "refresh"
  )

  invalid_refit_args <- intersect(
    names(refit_args),
    prohibited_refit_args
  )

  if (length(invalid_refit_args) > 0L) {
    stop(
      paste0(
        "Do not supply these arguments through 'refit_args': ",
        paste(invalid_refit_args, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (!is.null(x$data2) && length(x$data2) > 0L) {
    stop(
      "reloo_stanli() does not currently support brms models using 'data2'.",
      call. = FALSE
    )
  }

  if (is.null(loo)) {
    loo <- brms::loo(
      x,
      resp = resp,
      cores = 1L
    )
  }

  if (!inherits(loo, "loo")) {
    stop("'loo' must inherit from class 'loo'.", call. = FALSE)
  }

  if (is.null(loo$diagnostics$pareto_k)) {
    stop(
      "The supplied loo object does not contain Pareto-k diagnostics.",
      call. = FALSE
    )
  }

  data <- x$data

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (nrow(data) != nrow(loo$pointwise)) {
    stop(
      paste0(
        "The number of rows in x$data does not match the number of ",
        "observations in loo$pointwise."
      ),
      call. = FALSE
    )
  }

  if (check) {
    yhash_loo <- attr(loo, "yhash")
    yhash_fit <- brms:::hash_response(
      x,
      resp = resp
    )

    if (!brms:::is_equal(yhash_loo, yhash_fit)) {
      stop(
        paste0(
          "Response values used in 'loo' and 'x' do not match. ",
          "If this is a false positive, use check = FALSE."
        ),
        call. = FALSE
      )
    }
  }

  omitted_rows <- loo::pareto_k_ids(
    loo,
    threshold = k_threshold
  )

  if (length(omitted_rows) == 0L) {
    message(
      "No problematic observations found. Returning the original loo object."
    )

    attr(loo, "reloo_stanli") <- list(
      omitted_rows = omitted_rows,
      k_threshold = k_threshold,
      refit_args = refit_args
    )

    return(loo)
  }

  stored_call_arguments <- attr(x, "stanli_call_arguments")

  if (is.null(stored_call_arguments) || !is.list(stored_call_arguments)) {
    stop(
      paste0(
        "This fit lacks the stored brm_stanli() call specification required ",
        "for exact Stanli refits. Refit the original model with the current ",
        "brm_stanli() function."
      ),
      call. = FALSE
    )
  }

  stored_call_arguments$file <- NULL
  stored_call_arguments$file_refit <- "always"
  stored_call_arguments$refresh <- 0L
  stored_call_arguments$empty <- TRUE
  stored_call_arguments$pathfinder_init <- NULL

  message(
    length(omitted_rows),
    " problematic observation(s) found.",
    "\nThe model will be refit ",
    length(omitted_rows),
    " time(s) sequentially with Stanli."
  )

  held_out_log_lik <- vector(
    mode = "list",
    length = length(omitted_rows)
  )

  for (j in seq_along(omitted_rows)) {
    omitted_row <- omitted_rows[[j]]

    message(
      "\nFitting model ",
      j,
      " of ",
      length(omitted_rows),
      " (leaving out observation ",
      omitted_row,
      ")."
    )

    training_data <- data[-omitted_row, , drop = FALSE]
    held_out_data <- data[omitted_row, , drop = FALSE]

    fit_arguments <- utils::modifyList(
      stored_call_arguments,
      refit_args
    )

    fit_arguments$data <- training_data

    refit <- do.call(
      brm_stanli,
      fit_arguments
    )

    held_out_log_lik[[j]] <- do.call(
      brms::log_lik,
      utils::modifyList(
        list(
          object = refit,
          newdata = held_out_data,
          resp = resp,
          allow_new_levels = TRUE,
          sample_new_levels = "gaussian",
          combine = TRUE
        ),
        log_lik_args
      )
    )
  }

  elpd_loo <- vapply(
    held_out_log_lik,
    brm_stanli_log_mean_exp,
    numeric(1)
  )

  full_data_problem_rows <- data[
    omitted_rows,
    ,
    drop = FALSE
  ]

  full_log_lik <- do.call(
    brms::log_lik,
    utils::modifyList(
      list(
        object = x,
        newdata = full_data_problem_rows,
        resp = resp,
        allow_new_levels = TRUE,
        sample_new_levels = "gaussian",
        combine = TRUE
      ),
      log_lik_args
    )
  )

  hat_lpd <- apply(
    full_log_lik,
    2L,
    brm_stanli_log_mean_exp
  )

  p_loo <- hat_lpd - elpd_loo

  replacement_columns <- c(
    "elpd_loo",
    "p_loo",
    "looic"
  )

  loo$pointwise[
    omitted_rows,
    replacement_columns
  ] <- cbind(
    elpd_loo,
    p_loo,
    -2 * elpd_loo
  )

  updated_pointwise <- loo$pointwise[
    ,
    replacement_columns,
    drop = FALSE
  ]

  loo$estimates[
    replacement_columns,
    "Estimate"
  ] <- colSums(updated_pointwise)

  loo$estimates[
    replacement_columns,
    "SE"
  ] <- sqrt(
    nrow(loo$pointwise) *
      apply(updated_pointwise, 2L, stats::var)
  )

  loo$diagnostics$pareto_k[omitted_rows] <- 0

  attr(loo, "reloo_stanli") <- list(
    omitted_rows = omitted_rows,
    k_threshold = k_threshold,
    refit_args = refit_args
  )

  loo
}

brm_stanli_muffle_intermediate_pareto_warnings <- function(expr, muffle) {
  if (!muffle) {
    return(expr)
  }

  withCallingHandlers(
    expr,
    warning = function(warning) {
      if (grepl(
        "pareto_k >|Pareto k diagnostic values are too high",
        conditionMessage(warning),
        ignore.case = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

loo.brm_stanli_fit <- function(
  x,
  ...,
  compare = TRUE,
  resp = NULL,
  pointwise = FALSE,
  moment_match = FALSE,
  reloo = FALSE,
  k_threshold = 0.7,
  save_psis = FALSE,
  moment_match_args = list(),
  reloo_args = list(),
  model_names = NULL
) {
  call_name <- paste(
    deparse(substitute(x), width.cutoff = 500L),
    collapse = ""
  )

  dots <- list(...)

  if (!is.list(moment_match_args)) {
    stop("'moment_match_args' must be a list.", call. = FALSE)
  }

  if (!is.list(reloo_args)) {
    stop("'reloo_args' must be a list.", call. = FALSE)
  }

  if (is.null(model_names)) {
    model_names <- call_name
  }

  if (!is.character(model_names) ||
      length(model_names) != 1L ||
      is.na(model_names) ||
      !nzchar(model_names)) {
    stop(
      "'model_names' must be a single non-empty character string.",
      call. = FALSE
    )
  }

  can_use_stored_loo <-
    !isTRUE(moment_match) &&
    !isTRUE(reloo) &&
    is.null(resp) &&
    !isTRUE(pointwise) &&
    !isTRUE(save_psis) &&
    identical(k_threshold, 0.7) &&
    length(dots) == 0L

  ordinary_fit <- brm_stanli_unmark_fit(x)

  ordinary_loo_call <- c(
    list(
      x = ordinary_fit,
      compare = compare,
      resp = resp,
      pointwise = pointwise,
      moment_match = FALSE,
      reloo = FALSE,
      k_threshold = k_threshold,
      save_psis = save_psis,
      use_stored = can_use_stored_loo,
      model_names = model_names
    ),
    dots
  )

  compute_ordinary_loo <- function() {
    do.call(
      brms::loo,
      ordinary_loo_call
    )
  }

  loo_result <- brm_stanli_muffle_intermediate_pareto_warnings(
    compute_ordinary_loo(),
    muffle = isTRUE(moment_match) || isTRUE(reloo)
  )

  attr(loo_result, "model_name") <- model_names

  if (isTRUE(moment_match)) {
    brm_stanli_require_live_model(x)

    moment_match_call <- utils::modifyList(
      list(
        x = x,
        loo = loo_result,
        k_threshold = k_threshold,
        resp = resp,
        cores = 1L
      ),
      moment_match_args
    )

    moment_match_call$cores <- 1L

    loo_result <- brm_stanli_muffle_intermediate_pareto_warnings(
      do.call(loo_moment_match_stanli, moment_match_call),
      muffle = isTRUE(reloo)
    )

    attr(loo_result, "model_name") <- model_names
  }

  if (isTRUE(reloo)) {
    brm_stanli_require_live_model(x)

    reloo_call <- utils::modifyList(
      list(
        x = x,
        loo = loo_result,
        k_threshold = k_threshold,
        resp = resp
      ),
      reloo_args
    )

    loo_result <- do.call(
      reloo_stanli,
      reloo_call
    )

    attr(loo_result, "model_name") <- model_names
  }

  loo_result
}

loo_moment_match.brm_stanli_fit <- function(
  x,
  loo = NULL,
  k_threshold = 0.7,
  newdata = NULL,
  resp = NULL,
  check = TRUE,
  recompile = FALSE,
  ...
) {
  brm_stanli_require_live_model(x)

  if (isTRUE(recompile)) {
    stop(
      paste0(
        "'recompile = TRUE' is not supported for brm_stanli_fit objects. ",
        "RStan recompilation would replace the live Stanli handle."
      ),
      call. = FALSE
    )
  }

  dots <- list(...)

  if (!is.null(dots$cores) &&
      (!is.numeric(dots$cores) ||
       length(dots$cores) != 1L ||
       as.integer(dots$cores) != 1L)) {
    stop(
      "'cores' must be 1 for Stanli moment matching.",
      call. = FALSE
    )
  }

  if (is.null(loo)) {
    loo <- brms::loo(
      x,
      resp = resp,
      cores = 1L
    )
  }

  do.call(
    loo_moment_match_stanli,
    utils::modifyList(
      list(
        x = x,
        loo = loo,
        k_threshold = k_threshold,
        newdata = newdata,
        resp = resp,
        check = check,
        cores = 1L
      ),
      dots
    )
  )
}

reloo.brm_stanli_fit <- function(
  x,
  loo = NULL,
  k_threshold = 0.7,
  newdata = NULL,
  resp = NULL,
  check = TRUE,
  recompile = NULL,
  future_args = list(),
  ...
) {
  brm_stanli_require_live_model(x)

  if (!is.null(newdata)) {
    stop(
      paste0(
        "'newdata' is not currently supported by reloo.brm_stanli_fit(). ",
        "Use the original fitting data."
      ),
      call. = FALSE
    )
  }

  if (isTRUE(recompile)) {
    stop(
      paste0(
        "'recompile = TRUE' is not supported for brm_stanli_fit objects. ",
        "Stanli refits are created directly with brm_stanli()."
      ),
      call. = FALSE
    )
  }

  if (!is.list(future_args)) {
    stop("'future_args' must be a list.", call. = FALSE)
  }

  if (length(future_args) > 0L) {
    warning(
      paste0(
        "'future_args' is ignored. Stanli exact LOO refits run sequentially ",
        "so live native handles are not sent to parallel workers."
      ),
      call. = FALSE
    )
  }

  if (is.null(loo)) {
    loo <- brms::loo(
      x,
      resp = resp,
      cores = 1L
    )
  }

  dots <- list(...)

  refit_args <- if (!is.null(dots$refit_args)) {
    dots$refit_args
  } else {
    list()
  }

  log_lik_args <- if (!is.null(dots$log_lik_args)) {
    dots$log_lik_args
  } else {
    list()
  }

  unsupported_arguments <- setdiff(
    names(dots),
    c(
      "refit_args",
      "log_lik_args"
    )
  )

  if (length(unsupported_arguments) > 0L) {
    stop(
      paste0(
        "Unsupported arguments for reloo.brm_stanli_fit(): ",
        paste(unsupported_arguments, collapse = ", "),
        ". Use 'refit_args' or 'log_lik_args' instead."
      ),
      call. = FALSE
    )
  }

  reloo_stanli(
    x = x,
    loo = loo,
    k_threshold = k_threshold,
    resp = resp,
    check = check,
    refit_args = refit_args,
    log_lik_args = log_lik_args
  )
}

add_criterion.brm_stanli_fit <- function(
  x,
  criterion,
  model_name = NULL,
  overwrite = FALSE,
  file = NULL,
  force_save = FALSE,
  ...
) {
  brm_stanli_require_namespace("brms")

  if (!inherits(x, "brm_stanli_fit")) {
    stop("'x' must inherit from 'brm_stanli_fit'.", call. = FALSE)
  }

  call_name <- paste(
    deparse(substitute(x), width.cutoff = 500L),
    collapse = ""
  )

  if (is.null(model_name)) {
    model_name <- call_name
  }

  if (!is.character(model_name) ||
      length(model_name) != 1L ||
      is.na(model_name) ||
      !nzchar(model_name)) {
    stop(
      "'model_name' must be a single non-empty character string.",
      call. = FALSE
    )
  }

  criterion <- unique(as.character(criterion))

  if (any(criterion == "R2")) {
    warning(
      "Criterion 'R2' is deprecated. Using 'bayes_R2' instead.",
      call. = FALSE
    )

    criterion[criterion == "R2"] <- "bayes_R2"
  }

  supported_criteria <- c(
    "loo",
    "waic",
    "bayes_R2",
    "loo_R2"
  )

  unsupported_criteria <- setdiff(
    criterion,
    supported_criteria
  )

  if (length(unsupported_criteria) > 0L) {
    stop(
      paste0(
        "Unsupported criterion/criteria for brm_stanli_fit: ",
        paste(unsupported_criteria, collapse = ", "),
        ". Supported criteria are: ",
        paste(supported_criteria, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (!is.logical(overwrite) ||
      length(overwrite) != 1L ||
      is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(force_save) ||
      length(force_save) != 1L ||
      is.na(force_save)) {
    stop("'force_save' must be TRUE or FALSE.", call. = FALSE)
  }

  auto_save <- FALSE

  if (!is.null(file)) {
    file <- brm_stanli_validate_cache_path(file)
  } else {
    file <- x$file

    if (!is.null(file)) {
      auto_save <- TRUE
    }
  }

  if (isTRUE(overwrite)) {
    new_criteria <- criterion
  } else {
    new_criteria <- criterion[
      vapply(
        criterion,
        function(name) is.null(x$criteria[[name]]),
        logical(1)
      )
    ]
  }

  x$criteria[new_criteria] <- NULL

  dots <- list(...)

  if ("loo" %in% new_criteria) {
    x$criteria$loo <- do.call(
      brms::loo,
      c(
        list(
          x = x,
          model_names = model_name
        ),
        dots
      )
    )
  }

  if ("waic" %in% new_criteria) {
    x$criteria$waic <- do.call(
      brms::waic,
      c(
        list(
          x = brm_stanli_unmark_fit(x),
          model_names = model_name
        ),
        dots
      )
    )
  }

  if ("bayes_R2" %in% new_criteria) {
    x$criteria$bayes_R2 <- do.call(
      brms::bayes_R2,
      c(
        list(
          object = x,
          summary = FALSE
        ),
        dots
      )
    )
  }

  if ("loo_R2" %in% new_criteria) {
    x$criteria$loo_R2 <- do.call(
      brms::loo_R2,
      c(
        list(
          object = x,
          summary = FALSE
        ),
        dots
      )
    )
  }

  if (!is.null(file) &&
      (isTRUE(force_save) || length(new_criteria) > 0L)) {
    if (isTRUE(auto_save)) {
      message(
        "Automatically saving the model object in '",
        file,
        "'."
      )
    }

    x$file <- file

    saveRDS(
      x,
      file = file
    )
  }

  x
}

update.brm_stanli_fit <- function(
  object,
  formula.,
  newdata = NULL,
  recompile = NULL,
  ...
) {
  brm_stanli_require_namespace("brms")
  brm_stanli_require_live_model(object)

  dots <- list(...)

  if ("data" %in% names(dots)) {
    stop(
      "Use 'newdata' rather than 'data' when updating a brm_stanli_fit object.",
      call. = FALSE
    )
  }

  if (!is.null(recompile)) {
    if (!is.logical(recompile) ||
        length(recompile) != 1L ||
        is.na(recompile)) {
      stop("'recompile' must be NULL, TRUE, or FALSE.", call. = FALSE)
    }

    warning(
      paste0(
        "'recompile' is ignored for brm_stanli_fit objects. Stanli always ",
        "constructs a fresh native model from the updated Stan code and data."
      ),
      call. = FALSE
    )
  }

  stored_arguments <- attr(object, "stanli_call_arguments")

  if (is.null(stored_arguments) || !is.list(stored_arguments)) {
    stop(
      paste0(
        "This fit does not contain the stored brm_stanli() call specification ",
        "required for update(). Refit it with the current brm_stanli() wrapper."
      ),
      call. = FALSE
    )
  }

  if (missing(formula.) || is.null(formula.)) {
    updated_formula <- object$formula
  } else {
    updated_formula <- stats::update(
      object$formula,
      formula.
    )
  }

  if ("family" %in% names(dots) || "autocor" %in% names(dots)) {
    updated_family <- if ("family" %in% names(dots)) {
      dots$family
    } else {
      object$family
    }

    updated_autocor <- if ("autocor" %in% names(dots)) {
      dots$autocor
    } else {
      object$autocor
    }

    updated_formula <- brms::bf(
      updated_formula,
      family = updated_family,
      autocor = updated_autocor
    )

    dots$family <- NULL
    dots$autocor <- NULL
  }

  if (is.null(newdata)) {
    updated_data <- object$data
    updated_data_name <- brm_stanli_data_name(updated_data)
  } else {
    updated_data <- newdata
    updated_data_name <- paste(
      deparse(substitute(newdata), width.cutoff = 500L),
      collapse = ""
    )
  }

  prohibited_arguments <- c(
    "formula",
    "data",
    "fit",
    "backend",
    "algorithm",
    "empty",
    "rename",
    "file",
    "file_refit"
  )

  invalid_arguments <- intersect(
    names(dots),
    prohibited_arguments
  )

  if (length(invalid_arguments) > 0L) {
    stop(
      paste0(
        "Do not supply these arguments to update.brm_stanli_fit(): ",
        paste(invalid_arguments, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  updated_arguments <- stored_arguments
  updated_arguments$formula <- updated_formula
  updated_arguments$data <- updated_data

  if (length(dots) > 0L) {
    updated_arguments[names(dots)] <- dots
  }

  updated_arguments$file <- NULL
  updated_arguments$file_refit <- "always"
  updated_arguments$empty <- TRUE

  updated_fit <- do.call(
    brm_stanli,
    updated_arguments
  )

  attr(updated_fit$data, "data_name") <- updated_data_name

  updated_fit
}

combine_models_stanli <- function(...) {
  stop(
    paste0(
      "combine_models_stanli() is not implemented yet. Do not use ",
      "brms::combine_models() with brm_stanli_fit objects because the stock ",
      "brms implementation combines fits through rstan::sflist2stanfit(), ",
      "which loses Stanli's live model adapter."
    ),
    call. = FALSE
  )
}

kfold.brm_stanli_fit <- function(...) {
  stop(
    paste0(
      "kfold() is not implemented for brm_stanli_fit objects yet. ",
      "Use loo(..., moment_match = TRUE, reloo = TRUE) instead when ",
      "appropriate."
    ),
    call. = FALSE
  )
}

bridge_sampler.brm_stanli_fit <- function(samples, ...) {
  stop(
    paste0(
      "bridge_sampler() is not implemented for brm_stanli_fit objects yet. ",
      "The stock brms method calls brms:::update_misc_env(), which can replace ",
      "the live Stanli model handle with RStan internals."
    ),
    call. = FALSE
  )
}

bayes_factor.brm_stanli_fit <- function(x1, x2, ...) {
  stop(
    paste0(
      "bayes_factor() is not implemented for brm_stanli_fit objects yet ",
      "because it depends on bridge sampling."
    ),
    call. = FALSE
  )
}

post_prob.brm_stanli_fit <- function(x, ...) {
  stop(
    paste0(
      "post_prob() is not implemented for brm_stanli_fit objects yet ",
      "because it depends on bridge sampling."
    ),
    call. = FALSE
  )
}
