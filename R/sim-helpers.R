################################################################################
# Simulation utilities for the marked Hawkes process, built on top of the
# existing sim_mhp() thinning simulator. Used to perform a simulation study of 
# the developed sampler. 
################################################################################

#' Ground intensity of the marked Hawkes process, for use with sim_mhp()
#'
#' Computes \eqn{\lambda_0 + \sum_j A e^{\beta M_j} f(t - T_j)}, the ground
#' intensity at time \code{t} given event history \code{(times, marks)}.
#' Matches the same kernel-evaluation logic as \code{\link{update_zs}} and
#' \code{\link{compensator}} (basis functions evaluated via \code{outer()}
#' over history x components, then combined via the raw stick-breaking
#' weights \code{w} and normalising constant \code{C}), just returning the
#' intensity itself rather than an integrated/compensator form.
#'
#' Intended to be passed as \code{sim_mhp()}'s \code{intensity_func}
#' argument, not called directly in most workflows -- see
#' \code{\link{simulate_mhp}}.
#'
#' @param t Numeric scalar; current time.
#' @param times Numeric vector; event times strictly before \code{t}
#'   (supplied by \code{sim_mhp()}'s thinning loop).
#' @param marks Numeric vector, same length as \code{times}; marks of those
#'   events.
#' @param lambda0 Numeric scalar; baseline intensity.
#' @param A Numeric scalar; mark productivity parameter.
#' @param beta Numeric scalar; optional mark productivity parameter, not 
#'   required if \code{mark_productivity = "linear"}
#' @param theta Numeric vector of length K; kernel atom locations.
#' @param w Numeric vector of length K; raw stick-breaking weights (summing
#'   to 1) -- not yet divided by \code{C}.
#' @param C Numeric scalar; kernel normalising constant (as in
#'   \code{\link{run_sampler}}'s \code{C_fun}).
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}
#'
#' @return Numeric scalar; the ground intensity at \code{t}.
#'
#' @export
mhp_intensity <- function(t, times, marks, lambda0, A, beta = NULL, theta, w, C,
                          kernel = c("step", "pwlin"), 
                          mark_productivity = c("linear", "exponential")) {
  
  kernel <- match.arg(kernel)
  mark_productivity = match.arg(mark_productivity)
  
  if (length(times) == 0) {
    return(lambda0)
  }
  
  dt <- t - times   # all > 0, since sim_mhp() only ever passes accepted history
  
  if (kernel == "step") {
    basis <- outer(dt, theta, function(d, th) as.numeric(d < th))
  } else {
    basis <- outer(dt, theta, function(d, th) pmax(th - d, 0))
  }
  f_vals <- as.numeric(basis %*% w) / C   # f(dt_j) for each past event j
  
  if (mark_productivity == "linear") {
    lambda0 + A * sum(exp(beta * marks) * f_vals)
  } else if (mark_productivity == "exponential") {
    lambda0 + A * sum(marks * f_vals)
  }
}


#' Simulate a marked Hawkes process from a parameter set
#'
#' Thin convenience wrapper around \code{sim_mhp()}: unpacks a parameter
#' list (in the same shape used elsewhere in this codebase --
#' \code{\link{run_sampler}}'s \code{init}, or one row of a posterior
#' samples matrix), builds the stick-breaking weights and normalising
#' constant once, and calls \code{sim_mhp()} with \code{\link{mhp_intensity}}
#' and an iid \code{Exp(gamma)} mark generator -- matching the "marks are
#' unpredictable" assumption in the generative model.
#'
#' @param T_max Numeric scalar; length of the simulation window.
#' @param params A list with elements \code{lambda0}, \code{A}, \code{beta},
#'   \code{theta} (length K), \code{v} (length K-1), \code{gamma}. Matches
#'   the shape of \code{\link{run_sampler}}'s \code{init} argument (minus
#'   \code{alpha}/\code{phi}, which govern the priors, not the generative
#'   process itself).
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#'
#' @return An object of class \code{c("marked_pp_sim", "point_process_sim")},
#'   as returned by \code{sim_mhp()} -- \code{$events} (times) and
#'   \code{$marks} are what you'll feed into \code{\link{run_mcmc}} for
#'   refitting.
#'
#' @examples
#' \dontrun{
#' true_params <- list(
#'   lambda0 = 0.07, A = 0.08, beta = 1.8,
#'   theta = sort(rexp(10, 0.03)), v = rbeta(9, 1, 0.5),
#'   gamma = 2.2
#' )
#' sim <- simulate_mhp(T_max = 500, params = true_params, kernel = "pwlin")
#' fit <- run_mcmc(sim$events, sim$marks, T_max = 500, kernel = "pwlin", K = 10)
#' }
#'
#' @export
simulate_mhp <- function(T_max, params, kernel = c("step", "pwlin"), 
                         mark_productivity = c("linear", "exponential")) {
  
  kernel <- match.arg(kernel)
  mark_productivity = match.arg(mark_productivity)
  
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


#' Build the "true" kernel function from a set of parameters
#'
#' Constructs \eqn{f(x) = \sum_k w_k g_k(x)} from a fixed parameter set -- for 
#' use as the \code{true_kernel} argument to \code{\link{plot_hawkes_kernel}} 
#' when checking kernel-shape recovery in a simulation study.
#'
#' @param params A list with elements \code{theta} (length K), 
#'   \code{v} (length K-1) -- as in \code{\link{simulate_mhp}}.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#'
#' @return A function of \code{x} (vectorised), giving the true kernel value.
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


#' Run one simulate-and-refit replicate
#'
#' Simulates one synthetic catalog from \code{true_params} via
#' \code{\link{simulate_mhp}}, then refits it with \code{\link{run_mcmc}}.
#' A single building block for \code{\link{run_simulation_study}}; call
#' directly if you just want to inspect one replicate (e.g. for the kernel
#' recovery plot).
#'
#' @param T_max Numeric scalar; simulation/fit window length.
#' @param true_params A list as expected by \code{\link{simulate_mhp}}
#'   (\code{lambda0, A, beta (if "exponential" productivity), theta, v, gamma}).
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#' @param K Integer; truncation level used when *refitting* (need not equal
#'   \code{length(true_params$theta)} -- e.g. you can simulate from a
#'   simpler true kernel and fit with more components, or vice versa).
#' @param n_chains,n_iter As in \code{\link{run_mcmc}}. Defaults here are
#'   lighter than \code{run_mcmc}'s own defaults (2 chains, 4000 iterations)
#'   since a simulation study cares about calibration across many
#'   replicates rather than a single maximally-precise fit -- increase if
#'   individual replicate Rhats look poor.
#' @param fit_seed,sim_seed Optional integer seeds for the MCMC chains and
#'   for the data simulation respectively -- kept separate since they seed
#'   different RNG mechanisms (\code{parallel::clusterSetRNGStream} vs the
#'   global RNG).
#' @param prior_params,proposal_sds As in \code{\link{run_mcmc}}.
#' @param progress Logical; passed to \code{\link{run_mcmc}}.
#'
#' @return A list with \code{sim} (the simulated \code{point_process_sim}
#'   object) and \code{fit} (the \code{\link{run_mcmc}} output).
#'
#' @export
run_mhp_replicate <- function(T_max, true_params, kernel = c("step", "pwlin"), 
                              mark_productivity = c("linear", "exponential"), K,
                              n_chains = 2, n_iter = 4000,
                              fit_seed = NULL, sim_seed = NULL,
                              prior_params = default_prior_params(),
                              proposal_sds = default_proposal_sds(),
                              progress = FALSE) {
  
  kernel <- match.arg(kernel)
  mark_productivity = match.arg(mark_productivity)
  
  if (!is.null(sim_seed)) set.seed(sim_seed)
  sim <- simulate_mhp(T_max = T_max, params = true_params, kernel = kernel, 
                      mark_productivity = mark_productivity)
  
  fit <- run_mcmc(
    times = sim$events, marks = sim$marks, T_max = T_max, kernel = kernel, 
    mark_productivity = mark_productivity, K = K, n_chains = n_chains, 
    n_iter = n_iter, seed = fit_seed, prior_params = prior_params, 
    proposal_sds = proposal_sds, progress = progress
  )
  
  list(sim = sim, fit = fit)
}


#' Run a full simulation-study scenario (many replicates)
#'
#' Repeats \code{\link{run_mhp_replicate}} \code{R} times under a fixed
#' \code{true_params} scenario, and summarises bias, RMSE, and credible
#' interval coverage for \code{lambda0}, \code{A}, \code{beta}, \code{gamma}
#' -- the scalars with an unambiguous "true value" (kernel-shape recovery,
#' for \code{theta}/\code{v}, is a separate check -- see the example below
#' using \code{\link{true_kernel_fn}} with \code{\link{plot_hawkes_kernel}}
#' on one replicate's fit, since that check is already robust to
#' mixture-component label switching, whereas a per-component
#' theta_k/v_k coverage table would not be).
#'
#' Replicates that error (simulation or fitting) are dropped with a warning
#' rather than aborting the whole study.
#'
#' @param T_max Numeric scalar; simulation/fit window length, shared across
#'   replicates.
#' @param true_params A list as in \code{\link{simulate_mhp}}.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#' @param K Integer; truncation level for refitting.
#' @param R Integer; number of replicates. Default 20.
#' @param n_chains,n_iter As in \code{\link{run_mhp_replicate}}.
#' @param base_seed Integer; replicate \code{r} uses seed
#'   \code{base_seed + r} for both simulation and fitting.
#' @param prior_params,proposal_sds As in \code{\link{run_mcmc}}.
#' @param ci_level Numeric in (0, 1); credible interval width used for the
#'   coverage check. Default 0.90.
#' @param progress Logical; passed through to each replicate's
#'   \code{\link{run_mcmc}} call.
#' @param verbose Logical; print replicate progress to the console.
#'
#' @return A list with:
#'   \describe{
#'     \item{summary}{A data frame: one row per parameter, with
#'       \code{true_value}, \code{mean_est} (mean of posterior means across
#'       replicates), \code{bias}, \code{rmse}, and \code{coverage} (empirical
#'       fraction of replicates whose \code{ci_level} credible interval
#'       contained the true value).}
#'     \item{mean_mat, lower_mat, upper_mat}{\code{R x 4} matrices of
#'       per-replicate posterior mean / lower / upper credible bound, for
#'       \code{lambda0, A, beta, gamma}.}
#'     \item{results}{List of length \code{R} (successful replicates only),
#'       each with \code{$fit} -- useful for e.g. picking one replicate for
#'       a kernel recovery plot.}
#'     \item{true_params, kernel, ci_level}{Echoed back for convenience.}
#'   }
#'
#' @examples
#' \dontrun{
#' true_params <- list(
#'   lambda0 = 0.07, A = 0.08, beta = 1.5,       # comfortably subcritical
#'   theta = sort(rexp(10, 0.03)), v = rbeta(9, 1, 0.5),
#'   gamma = 2.2
#' )
#' study <- run_simulation_study(T_max = 500, true_params = true_params,
#'                                kernel = "pwlin", K = 10, R = 20)
#' study$summary
#'
#' # kernel-shape recovery, on one representative replicate
#' plot_hawkes_kernel(study$results[[1]]$fit$chains, kernel = "pwlin",
#'                     true_kernel = true_kernel_fn(true_params, "pwlin"))
#' }
#'
#' @export
run_simulation_study <- function(T_max, true_params, kernel = c("step", "pwlin"), 
                                 mark_productivity = c("linear", "exponential"), K,
                                 R = 20, n_chains = 2, n_iter = 4000, n_burn=NULL,
                                 base_seed = 1,
                                 prior_params = default_prior_params(),
                                 proposal_sds = default_proposal_sds(),
                                 ci_level = 0.95, progress = FALSE, verbose = TRUE) {
  
  kernel <- match.arg(kernel)
  mark_productivity = match.arg(mark_productivity)
  
  if (mark_productivity == "linear") {
    scalar_names <- c("lambda0", "A", "gamma")
  } else if (mark_productivity == "exponential") {
    scalar_names <- c("lambda0", "A", "beta", "gamma")
  }
  true_vals <- unlist(true_params[scalar_names])
  tail_p <- (1 - ci_level) / 2
  
  results <- vector("list", R)
  
  for (r in seq_len(R)) {
    if (verbose) cat(sprintf("Replicate %d/%d...\n", r, R))
    
    rep_out <- tryCatch(
      run_mhp_replicate(
        T_max = T_max, true_params = true_params, kernel = kernel, 
        mark_productivity = mark_productivity, K = K,
        n_chains = n_chains, n_iter = n_iter,
        fit_seed = base_seed + r, sim_seed = base_seed + r,
        prior_params = prior_params, proposal_sds = proposal_sds, progress = progress
      ),
      error = function(e) {
        warning(sprintf("Replicate %d failed: %s", r, conditionMessage(e)))
        NULL
      }
    )
    
    if !is.null(n_burn) {
      # drop first n_burn samples if provided a non null value
      for (j in seq_along(rep_out$fit$chains)) {
        rep_out$fit$chains[[j]]$samples <-
          rep_out$fit$chains[[j]]$samples[
            (n_burn + 1):nrow(rep_out$fit$chains[[j]]$samples), , drop = FALSE
          ]
      }
    }
    
    if (is.null(rep_out)) next
    
    samples <- do.call(
      rbind, 
      lapply(rep_out$fit$chains, function(ch) ch$samples)
    )
    
    est <- sapply(scalar_names, function(nm) {
      x <- samples[, nm]
      c(mean  = mean(x),
        lower = as.numeric(quantile(x, tail_p)),
        upper = as.numeric(quantile(x, 1 - tail_p)))
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
  if (verbose) cat(sprintf("%d/%d replicates completed successfully.\n", R_ok, R))
  if (R_ok == 0) stop("All replicates failed -- check true_params/settings.")
  
  mean_mat  <- t(sapply(results, function(res) res$est["mean", ]))
  lower_mat <- t(sapply(results, function(res) res$est["lower", ]))
  upper_mat <- t(sapply(results, function(res) res$est["upper", ]))
  colnames(mean_mat) <- colnames(lower_mat) <- colnames(upper_mat) <- scalar_names
  
  covered <- sapply(scalar_names, function(nm) {
    (lower_mat[, nm] <= true_vals[nm]) & (true_vals[nm] <= upper_mat[, nm])
  })
  
  true_vals_rep <- matrix(true_vals, nrow(mean_mat), length(true_vals), byrow = TRUE)
  
  summary_tbl <- data.frame(
    parameter  = scalar_names,
    true_value = true_vals,
    mean_est   = colMeans(mean_mat),
    bias       = colMeans(mean_mat) - true_vals,
    rmse       = sqrt(colMeans((mean_mat - true_vals_rep)^2)),
    coverage   = colMeans(covered),
    row.names  = NULL
  )
  
  list(
    summary     = summary_tbl,
    mean_mat    = mean_mat,
    lower_mat   = lower_mat,
    upper_mat   = upper_mat,
    results     = results,
    true_params = true_params,
    kernel      = kernel,
    ci_level    = ci_level
  )
}

#' Plot MCMC diagnostics for a simulation study replicate
#'
#' @param study List; output from \code{run_simulation_study}.
#' @param R Integer; number of individual replicate to diagnose.
#' @param start Integer; first iteration of interest.
#' @param thin Integer; the required interval between successive samples.
#' @param ... optional arguments to pass to MCMCvis::MCMCtrace
#' 
#' @export
diagnose_sim_study <- function(study, R, start = 1, thin = 1, params = "all", ...) {
  
  # Extract the MCMC chains object
  chains <- study$results[[R]]$fit$chains
  
  mcmc_chains <- lapply(chains, function(chain) {
    coda::mcmc(chain$samples)
  })
  mcmc_chains <- window(coda::mcmc.list(mcmc_chains), start = start, thin = thin)
  # mcmc_chains <- window(mcmc_chains, start = start, thin = thin)
  
  MCMCvis::MCMCtrace(mcmc_chains, params = params, ...)
  
  return(MCMCvis::MCMCsummary(mcmc_chains, round = 4, params = params))
}