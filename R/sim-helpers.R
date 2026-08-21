################################################################################
# Simulation utilities for the marked Hawkes process, built on top of the
# existing sim_mhp() thinning simulator. Used to perform a simulation study of
# the developed sampler.
################################################################################

#' Evaluate the ground intensity of a marked Hawkes process
#'
#' Computes the conditional ground intensity at time \code{t} for use with
#' \code{\link{sim_mhp}}.
#'
#' @param t Numeric scalar; current time.
#' @param times Numeric vector; previous event times.
#' @param marks Numeric vector; marks corresponding to \code{times}.
#' @param lambda0 Numeric scalar; baseline intensity.
#' @param A Numeric scalar; mark productivity parameter.
#' @param beta Numeric scalar; exponential mark-productivity parameter.
#'   Ignored for \code{"linear"} productivity.
#' @param theta Numeric vector; kernel atom locations.
#' @param w Numeric vector; raw stick-breaking weights.
#' @param C Numeric scalar; kernel normalising constant.
#' @param kernel Character; kernel type, \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; mark-productivity form,
#'   \code{"linear"} or \code{"exponential"}.
#'
#' @return Numeric scalar; conditional ground intensity at \code{t}.
#'
#' @export
mhp_intensity <- function(
  t,
  times,
  marks,
  lambda0,
  A,
  beta = NULL,
  theta,
  w,
  C,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential")
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  if (length(times) == 0) {
    return(lambda0)
  }

  dt <- t - times # all > 0, since sim_mhp() only ever passes accepted history

  if (kernel == "step") {
    basis <- outer(dt, theta, function(d, th) as.numeric(d < th))
  } else {
    basis <- outer(dt, theta, function(d, th) pmax(th - d, 0))
  }
  f_vals <- as.numeric(basis %*% w) / C # f(dt_j) for each past event j

  if (mark_productivity == "linear") {
    lambda0 + A * sum(marks * f_vals)
  } else if (mark_productivity == "exponential") {
    lambda0 + A * sum(exp(beta * marks) * f_vals)
  }
}


#' Simulate a marked Hawkes process from model parameters
#'
#' Wraps \code{\link{sim_mhp}} using the supplied parameters, with
#' stick-breaking weights and kernel normalisation computed internally.
#' Marks are generated independently as \eqn{\mathrm{Exp}(\gamma)}.
#'
#' @param T_max Numeric scalar; simulation window length.
#' @param params List containing \code{lambda0}, \code{A}, \code{beta},
#'   \code{theta}, \code{v}, and \code{gamma}.
#' @param kernel Character; kernel type, \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; mark-productivity form,
#'   \code{"linear"} or \code{"exponential"}.
#'
#' @return An object of class \code{c("marked_pp_sim", "point_process_sim")}
#'   containing the simulated event times and marks.
#'
#' @examples
#' \dontrun{
#' true_params <- list(
#'   lambda0 = 0.07, A = 0.08, beta = 1.8,
#'   theta = sort(rexp(10, 0.03)), v = rbeta(9, 1, 0.5),
#'   gamma = 2.2
#' )
#' sim <- simulate_mhp(T_max = 500, params = true_params, kernel = "pwlin")
#' }
#'
#' @export
simulate_mhp <- function(
  T_max,
  params,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential")
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  K <- length(params$theta)
  stopifnot(length(params$v) == K - 1)

  remaining <- cumprod(c(1, 1 - params$v))
  w <- c(params$v, 1) * remaining

  if (kernel == "step") {
    C <- sum(w * params$theta)
  } else {
    C <- 0.5 * sum(w * params$theta^2)
  }

  sim_mhp(
    T_max = T_max,
    intensity_func = mhp_intensity,
    rmark_func = function() stats::rexp(1, params$gamma),
    lambda0 = params$lambda0,
    A = params$A,
    beta = params$beta,
    theta = params$theta,
    w = w,
    C = C,
    kernel = kernel,
    mark_productivity = mark_productivity
  )
}


#' Construct the kernel function from model parameters
#'
#' Returns the normalised kernel implied by the supplied atom locations and
#' stick-breaking weights.
#'
#' @param params List containing \code{theta} and \code{v}.
#' @param kernel Character; kernel type, \code{"step"} or \code{"pwlin"}.
#'
#' @return A vectorised function of \code{x} returning the kernel value.
#'
#' @export
true_kernel_fn <- function(params, kernel = c("step", "pwlin")) {
  kernel <- match.arg(kernel)

  remaining <- cumprod(c(1, 1 - params$v))
  w_raw <- c(params$v, 1) * remaining

  if (kernel == "step") {
    C <- sum(w_raw * params$theta)
  } else {
    C <- 0.5 * sum(w_raw * params$theta^2)
  }
  w <- w_raw / C

  function(x) {
    if (kernel == "step") {
      basis <- outer(x, params$theta, function(xx, th) as.numeric(xx < th))
    } else {
      basis <- outer(x, params$theta, function(xx, th) pmax(th - xx, 0))
    }
    as.numeric(basis %*% w)
  }
}


#' Simulate and refit one marked Hawkes process replicate
#'
#' Simulates one dataset from \code{true_params} and fits the model using
#' \code{\link{run_mcmc}}.
#'
#' @param T_max Numeric scalar; simulation and fitting window length.
#' @param true_params List containing the true simulation parameters:
#'   \code{lambda0}, \code{A}, \code{beta}, \code{theta}, \code{v}, and
#'   \code{gamma}.
#' @param kernel Character; kernel type, \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; mark-productivity form,
#'   \code{"linear"} or \code{"exponential"}.
#' @param K Integer; number of kernel components used when fitting.
#' @param n_chains Integer; number of MCMC chains.
#' @param n_iter Integer; number of iterations per chain.
#' @param fit_seed Integer or \code{NULL}; seed for MCMC fitting.
#' @param sim_seed Integer or \code{NULL}; seed for data simulation.
#' @param prior_params List of prior parameters passed to \code{\link{run_mcmc}}.
#' @param proposal_sds List of proposal standard deviations passed to
#'   \code{\link{run_mcmc}}.
#' @param progress Logical; whether to display MCMC progress.
#'
#' @return A list containing the simulated data in \code{$sim} and the fitted
#'   model in \code{$fit}.
#'
#' @export
run_mhp_replicate <- function(
  T_max,
  true_params,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential"),
  K,
  n_chains = 2,
  n_iter = 4000,
  fit_seed = NULL,
  sim_seed = NULL,
  prior_params = default_prior_params(),
  proposal_sds = default_proposal_sds(),
  progress = FALSE
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  if (!is.null(sim_seed)) {
    set.seed(sim_seed)
  }
  sim <- simulate_mhp(
    T_max = T_max,
    params = true_params,
    kernel = kernel,
    mark_productivity = mark_productivity
  )

  fit <- run_mcmc(
    times = sim$events,
    marks = sim$marks,
    T_max = T_max,
    kernel = kernel,
    mark_productivity = mark_productivity,
    K = K,
    n_chains = n_chains,
    n_iter = n_iter,
    seed = fit_seed,
    prior_params = prior_params,
    proposal_sds = proposal_sds,
    progress = progress
  )

  list(sim = sim, fit = fit)
}


#' Run a marked Hawkes process simulation study
#'
#' Simulates and fits \code{R} replicates under a fixed parameter scenario,
#' then summarises bias, RMSE, and credible-interval coverage for the scalar
#' model parameters.
#'
#' Failed replicates are omitted with a warning.
#'
#' @param T_max Numeric scalar; simulation and fitting window length.
#' @param true_params List of true simulation parameters.
#' @param kernel Character; kernel type, \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; mark-productivity form,
#'   \code{"linear"} or \code{"exponential"}.
#' @param K Integer; number of kernel components used when fitting.
#' @param R Integer; number of simulation replicates.
#' @param n_chains Integer; number of MCMC chains per replicate.
#' @param n_iter Integer; number of MCMC iterations per chain.
#' @param n_burn Integer or \code{NULL}; number of initial samples to discard
#'   from each chain before summarising.
#' @param base_seed Integer; base seed used to initialise each replicate.
#' @param prior_params List of prior parameters passed to \code{\link{run_mcmc}}.
#' @param proposal_sds List of proposal standard deviations passed to
#'   \code{\link{run_mcmc}}.
#' @param ci_level Numeric in \code{(0, 1)}; credible interval level.
#' @param progress Logical; whether to display MCMC progress.
#' @param verbose Logical; whether to print replicate progress.
#'
#' @return A list containing:
#'   \describe{
#'     \item{summary}{Data frame containing the true value, mean estimate,
#'       bias, RMSE, and coverage for each scalar parameter.}
#'     \item{mean_mat}{Matrix of posterior means by replicate.}
#'     \item{lower_mat}{Matrix of lower credible bounds by replicate.}
#'     \item{upper_mat}{Matrix of upper credible bounds by replicate.}
#'     \item{results}{List of successful replicate results.}
#'     \item{true_params}{The supplied true parameter set.}
#'     \item{kernel}{The kernel type used in the study.}
#'     \item{ci_level}{The credible interval level used for coverage.}
#'   }
#'
#' @export
run_simulation_study <- function(
  T_max,
  true_params,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential"),
  K,
  R = 20,
  n_chains = 2,
  n_iter = 4000,
  n_burn = NULL,
  base_seed = 1,
  prior_params = default_prior_params(),
  proposal_sds = default_proposal_sds(),
  ci_level = 0.95,
  progress = FALSE,
  verbose = TRUE
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  if (mark_productivity == "linear") {
    scalar_names <- c("lambda0", "A", "gamma")
  } else if (mark_productivity == "exponential") {
    scalar_names <- c("lambda0", "A", "beta", "gamma")
  }
  true_vals <- unlist(true_params[scalar_names])
  tail_p <- (1 - ci_level) / 2

  results <- vector("list", R)

  for (r in seq_len(R)) {
    if (verbose) {
      cat(sprintf("Replicate %d/%d...\n", r, R))
    }

    rep_out <- tryCatch(
      run_mhp_replicate(
        T_max = T_max,
        true_params = true_params,
        kernel = kernel,
        mark_productivity = mark_productivity,
        K = K,
        n_chains = n_chains,
        n_iter = n_iter,
        fit_seed = base_seed + r,
        sim_seed = base_seed + r,
        prior_params = prior_params,
        proposal_sds = proposal_sds,
        progress = progress
      ),
      error = function(e) {
        warning(sprintf("Replicate %d failed: %s", r, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(rep_out)) {
      next
    }

    if (!is.null(n_burn)) {
      for (j in seq_along(rep_out$fit$chains)) {
        rep_out$fit$chains[[j]]$samples <-
          rep_out$fit$chains[[j]]$samples[
            (n_burn + 1):nrow(rep_out$fit$chains[[j]]$samples),
            ,
            drop = FALSE
          ]
      }
    }

    samples <- do.call(
      rbind,
      lapply(rep_out$fit$chains, function(ch) ch$samples)
    )

    est <- sapply(scalar_names, function(nm) {
      x <- samples[, nm]
      c(
        mean = mean(x),
        lower = as.numeric(quantile(x, tail_p)),
        upper = as.numeric(quantile(x, 1 - tail_p))
      )
    })

    results[[r]] <- list(
      n_events = rep_out$sim$n,
      events = rep_out$sim$events,
      marks = rep_out$sim$marks,
      est = est,
      fit = rep_out$fit
    )
  }

  ok <- !sapply(results, is.null)
  results <- results[ok]
  R_ok <- length(results)
  if (verbose) {
    cat(sprintf("%d/%d replicates completed successfully.\n", R_ok, R))
  }
  if (R_ok == 0) {
    stop("All replicates failed -- check true_params/settings.")
  }

  mean_mat <- t(sapply(results, function(res) res$est["mean", ]))
  lower_mat <- t(sapply(results, function(res) res$est["lower", ]))
  upper_mat <- t(sapply(results, function(res) res$est["upper", ]))
  colnames(mean_mat) <- colnames(lower_mat) <- colnames(
    upper_mat
  ) <- scalar_names

  covered <- sapply(scalar_names, function(nm) {
    (lower_mat[, nm] <= true_vals[nm]) & (true_vals[nm] <= upper_mat[, nm])
  })

  true_vals_rep <- matrix(
    true_vals,
    nrow(mean_mat),
    length(true_vals),
    byrow = TRUE
  )

  summary_tbl <- data.frame(
    parameter = scalar_names,
    true_value = true_vals,
    mean_est = colMeans(mean_mat),
    bias = colMeans(mean_mat) - true_vals,
    rmse = sqrt(colMeans((mean_mat - true_vals_rep)^2)),
    coverage = colMeans(covered),
    row.names = NULL
  )

  list(
    summary = summary_tbl,
    mean_mat = mean_mat,
    lower_mat = lower_mat,
    upper_mat = upper_mat,
    results = results,
    true_params = true_params,
    kernel = kernel,
    ci_level = ci_level
  )
}

#' Diagnose MCMC convergence for a simulation-study replicate
#'
#' Produces standard MCMC summaries and, optionally, trace plots. If
#' \code{lags} is supplied, also computes Gelman-Rubin statistics for the
#' posterior kernel evaluated at those lags.
#'
#' @param study List; output from \code{\link{run_simulation_study}}.
#' @param R Integer; replicate number to diagnose.
#' @param start Integer; first iteration included in the diagnostics.
#' @param thin Integer; thinning interval.
#' @param params Character vector; parameters to include in the standard
#'   MCMC diagnostics, or \code{"all"}.
#' @param plot Logical; whether to produce trace plots.
#' @param lags Numeric vector or \code{NULL}; lag values at which to evaluate
#'   the kernel for convergence diagnostics.
#' @param ... Additional arguments passed to
#'   \code{MCMCvis::MCMCtrace}.
#'
#' @return If \code{lags = NULL}, a parameter summary from
#'   \code{MCMCvis::MCMCsummary}. Otherwise, a list containing the parameter
#'   summary, kernel-specific Gelman-Rubin statistics, and the corresponding
#'   MCMC object.
#'
#' @export
diagnose_sim_study <- function(
  study,
  R,
  start = 1,
  thin = 1,
  params = "all",
  plot = FALSE,
  lags = NULL,
  ...
) {
  # Extract the MCMC chains object
  chains <- study$results[[R]]$fit$chains

  mcmc_chains <- lapply(chains, function(chain) {
    coda::mcmc(chain$samples)
  })

  mcmc_chains <- window(
    coda::mcmc.list(mcmc_chains),
    start = start,
    thin = thin
  )

  # Standard parameter diagnostics
  if (plot) {
    MCMCvis::MCMCtrace(mcmc_chains, params = params, ...)
  }

  summary <- MCMCvis::MCMCsummary(
    mcmc_chains,
    round = 4,
    params = params
  )

  # Compute Gelman-Rubin convergence statistics for the kernel values atspecified lags
  if (!is.null(lags)) {
    kernel <- study$kernel

    kernel_chains <- lapply(mcmc_chains, function(chain) {
      samples <- as.matrix(chain)

      apply(samples, 1, function(s) {
        theta <- s[grep("^theta", names(s))]
        v <- s[grep("^v", names(s))]

        # Stick-breaking weights
        remaining <- cumprod(c(1, 1 - v))
        w <- c(v, 1) * remaining

        # Kernel normalising constant
        if (kernel == "step") {
          C <- sum(w * theta)
          basis <- outer(
            lags,
            theta,
            function(x, th) as.numeric(x < th)
          )
        } else if (kernel == "pwlin") {
          C <- 0.5 * sum(w * theta^2)
          basis <- outer(
            lags,
            theta,
            function(x, th) pmax(th - x, 0)
          )
        } else {
          stop("Unknown kernel: ", kernel)
        }
        as.numeric(basis %*% w / C)
      })
    })

    # Convert each chain to an MCMC object
    kernel_chains <- lapply(kernel_chains, function(x) {
      x <- t(x)
      colnames(x) <- paste0("f_", lags)
      coda::mcmc(x)
    })

    kernel_mcmc <- coda::mcmc.list(kernel_chains)

    # Gelman-Rubin diagnostics
    gelman <- coda::gelman.diag(
      kernel_mcmc,
      autoburnin = FALSE,
      multivariate = FALSE
    )

    gelman_stats <- gelman$psrf[, "Point est."]

    names(gelman_stats) <- paste0("f(", lags, ")")

    kernel_summary <- data.frame(
      lag = lags,
      Rhat = as.numeric(gelman_stats)
    )

    return(list(
      parameter_summary = summary,
      kernel_Rhat = kernel_summary,
      kernel_mcmc = kernel_mcmc
    ))
  }

  return(summary)
}
