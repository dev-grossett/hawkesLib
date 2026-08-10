################################################################################
# Marked Hawkes process: MCMC sampler + posterior kernel plotting
#
# Model: a Hawkes process with a K-component truncated stick-breaking
# (Dirichlet process) mixture excitation kernel, extended to incorporate
# event marks m_i: a parent event's offspring intensity scales with its own
# via two possible parameterisations:
#   - Linear productivity      - A * mark
#   - Exponential productivity - A * exp(beta * mark) 
# In both cases, marks themselves are modelled as iid Exp(gamma), giving a 
# conjugate Gamma full conditional
#
# Public/exported entry points: run_mcmc(), plot_hawkes_kernel()
# Everything else is an internal helper (@noRd).
################################################################################

#' Sample latent branching structure and mixture-component allocations
#'
#' For each event, samples which prior event (if any) is its parent, and
#' which mixture component of the excitation kernel produced it, via a
#' single categorical draw over {immigrant, every valid (parent, component)
#' pair}. 
#'
#' @param lambda0 Numeric scalar; baseline (immigrant) intensity.
#' @param A Numeric scalar; mark productivity parameter 
#' @param beta Numeric scalar; optional mark productivity parameter, not 
#'   required if \code{mark_productivity = "linear"}
#' @param w Numeric vector of length K; stick-breaking mixture weights
#'   (raw weights that sum to 1, as constructed by the caller).
#' @param C Numeric scalar; kernel normalising constant, e.g.
#'   \code{sum(w * theta)} (step) or \code{0.5 * sum(w * theta^2)} (pwlin).
#' @param theta Numeric vector of length K; kernel atom locations.
#' @param times Numeric vector of event times, sorted ascending.
#' @param marks Numeric vector, same length as \code{times}; the mark
#'   (e.g. magnitude) of each event.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}
#'
#' @return A list with:
#'   \describe{
#'     \item{z}{Integer vector, same length as \code{times}. \code{z[i] == 0}
#'       marks event \code{i} as an immigrant; otherwise \code{z[i]} is the
#'       index of event \code{i}'s parent.}
#'     \item{s}{Integer vector, same length as \code{times}. The mixture
#'       component that produced event \code{i}, or \code{NA} for
#'       immigrants.}
#'   }
#'
#' @keywords internal
update_zs <- function(lambda0, A, beta = NULL, w, C, theta, times, marks, 
                      kernel = c("step", "pwlin"), 
                      mark_productivity = c("linear", "exponential")) {
  
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  stopifnot(
    "Need 'beta' param when mark_productivity = 'exponential'" = 
    mark_productivity == "exponential" & !is.null(beta)
  )
  
  n <- length(times)
  z <- integer(n)
  s <- integer(n)
  s[] <- NA_integer_
  
  if (n <= 1) return(list(z = z, s = s))
  
  for (i in 2:n) {
    # Time since each possible parent event
    dt <- times[i] - times[1:(i - 1)]
    
    # Matrix of basis-function values:
    # rows = possible parents j = i - 1
    # columns = mixture components k
    if (kernel == "step") {
      g <- dt < matrix(theta, nrow = i - 1, ncol = length(theta), byrow = TRUE)
    } else if (kernel == "pwlin") {
      g <- pmax(outer(dt, theta, FUN = function(dt, theta) {theta - dt}), 0)
    }
    
    # Find valid (parent, component) combinations
    active_flat_idx <- which(g > 0, arr.ind = TRUE)
    
    # Handle no valid (j, k)
    if (nrow(active_flat_idx) == 0) {
      z[i] <- 0
      s[i] <- NA_integer_
      next
    }
    
    # Map matrix coordinates back to parent and component
    parent_idx <- active_flat_idx[, 1]
    comp_idx   <- active_flat_idx[, 2]
    
    # Offspring probabilities
    if (mark_productivity == "linear") {
      offspring_probs <- A * marks[parent_idx] * w[comp_idx] * 
        g[active_flat_idx] / C
    } else if (mark_productivity == "exponential") {
      offspring_probs <- A * exp(beta * marks[parent_idx]) * w[comp_idx] *
        g[active_flat_idx] / C
    }    
    
    # Include immigrant probability
    # (clip at 0 as a defensive guard against floating-point cancellation
    # elsewhere producing a tiny negative weight; should not trigger once
    # the stick-breaking weights are computed without residual subtraction)
    probs <- pmax(c(lambda0, offspring_probs), 0)
    draw <- sample.int(length(probs), 1, prob = probs)
    
    if (draw > 1) {
      idx <- draw - 1
      z[i] <- parent_idx[idx]
      s[i] <- comp_idx[idx]
    }
  }
  
  return(list(z = z, s = s))
}

#' Full conditional log-posterior for a single kernel atom theta_k
#'
#' Evaluates the (unnormalised) log full conditional of \code{theta_k} on
#' the log scale, for use in a Metropolis-within-Gibbs update. Only offspring
#' events allocated to component \code{k} depend on \code{theta_k} through
#' the kernel basis function; offspring allocated to other components are
#' constant with respect to \code{theta_k} and are deliberately excluded from
#' both the validity check and the log-density sum, so that an unrelated
#' component's (in)validity can never leak a spurious \code{-Inf} into this
#' component's Metropolis ratio. The normalising constant \code{C}, by
#' contrast, genuinely depends on every component and is always recomputed
#' in full. (Unaffected by marks/A/beta -- those only enter through
#' \code{\link{update_zs}} and their own full conditionals, not through the
#' kernel shape.)
#' 
#' @param log_theta_k Numeric scalar; proposed value of \code{theta_k} on
#'   the log scale (this function samples on the log scale so a Gamma-type
#'   positive parameter can use a symmetric random-walk proposal).
#' @param k Integer; index of the component being updated.
#' @param w Numeric vector of length K; current stick-breaking weights.
#' @param theta Numeric vector of length K; current atom locations (entry
#'   \code{k} is overwritten internally by the proposed value).
#' @param phi Numeric scalar; rate parameter of the Gamma prior on
#'   \code{theta_k}.
#' @param times Numeric vector of all event times.
#' @param O_idx Integer vector; indices of offspring events (i.e.
#'   \code{which(z != 0)}).
#' @param z Integer vector; parent index for each event (0 for immigrants).
#' @param s Integer vector; mixture component allocation for each event
#'   (\code{NA} for immigrants).
#' @param C_fun Function of \code{(w, theta)} returning the kernel
#'   normalising constant.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#'
#' @return Numeric scalar: the log full conditional density (posterior plus
#'   log-Jacobian of the log transform) at \code{log_theta_k}, or
#'   \code{-Inf} if the proposal is incompatible with an offspring event's
#'   observed waiting time.
#'
#' @keywords internal
logpost_theta_k <- function(log_theta_k, k, w, theta, phi, times, O_idx,
                            z, s, C_fun, kernel = c("step", "pwlin")) {
  
  kernel <- match.arg(kernel)
  
  theta_k <- exp(log_theta_k)
  
  # Proposed atom locations
  theta_prop <- theta
  theta_prop[k] <- theta_k
  
  # Recalculate normalising constant (genuinely depends on every component)
  C_prop <- C_fun(w, theta_prop)
  
  # Latent mixture components for offspring events
  s_O <- s[O_idx]
  
  # Time since parent event
  dt_O <- times[O_idx] - times[z[O_idx]]
  
  # Restrict to offspring actually allocated to component k
  own    <- (s_O == k)
  dt_own <- dt_O[own]
  
  if (kernel == "step") {
    if (any(theta_k <= dt_own)) return(-Inf)
    log_g_own_sum <- 0  # log(1) for every satisfied step-indicator
  } else if (kernel == "pwlin") {
    g_own <- theta_k - dt_own
    if (any(g_own <= 0)) return(-Inf)
    log_g_own_sum <- sum(log(g_own))
  }
  
  # sum(log(w[s_O])) does not depend on theta_k and cannot be -Inf
  # (stick-breaking weights are always > 0), so it's safe to include in full
  ll_sum <- sum(log(w[s_O])) + log_g_own_sum - length(O_idx)*log(C_prop)
  
  # Return posterior + log-Jacobian
  return(ll_sum - phi*theta_k + log_theta_k)
}


#' Full conditional log-posterior for a single stick-breaking variable v_k
#'
#' Evaluates the (unnormalised) log full conditional of \code{v_k} on the
#' logit scale, for use in a Metropolis-within-Gibbs update. Stick-breaking
#' weight \code{w_m} depends on \code{v_k} through a factor of \code{v_k}
#' when \code{m == k}, through a common factor of \code{(1 - v_k)} when
#' \code{m > k}, and not at all when \code{m < k}. Rather than summing
#' \code{log(w_prop[s_O])} over every offspring event (which risks a
#' \code{-Inf / -Inf -> NaN} collision if some unrelated component has
#' already saturated to ~0 weight elsewhere in the same sweep), only the
#' part that truly depends on \code{v_k} is computed directly, using the
#' precomputed per-component offspring counts \code{n_comp} (allocated to
#' exactly component \code{k}) and \code{m_comp} (allocated to any component
#' above \code{k}). Unaffected by marks, same reasoning as
#' \code{\link{logpost_theta_k}}.
#'
#' @param logit_v_k Numeric scalar; proposed value of \code{v_k} on the
#'   logit scale.
#' @param k Integer; index of the stick-breaking variable being updated
#'   (\code{1, ..., K-1}).
#' @param v Numeric vector of length K-1; current stick-breaking variables
#'   (entry \code{k} is overwritten internally by the proposed value).
#' @param theta Numeric vector of length K; current atom locations.
#' @param alpha Numeric scalar; concentration parameter of the Beta(1,
#'   alpha) stick-breaking prior.
#' @param n_comp Integer vector of length K; number of offspring events
#'   allocated to exactly each component.
#' @param m_comp Integer vector of length K; number of offspring events
#'   allocated to any component with a strictly greater index.
#' @param n_off Integer scalar; total number of offspring events (i.e.
#'   \code{sum(n_comp)}).
#' @param C_fun Function of \code{(w, theta)} returning the kernel
#'   normalising constant.
#'
#' @return Numeric scalar: the log full conditional density (posterior plus
#'   log-Jacobian of the logit transform) at \code{logit_v_k}.
#'
#' @keywords internal
logpost_v_k <- function(logit_v_k, k, v, theta, alpha, 
                        n_comp, m_comp, n_off, C_fun) {
  v_k <- plogis(logit_v_k)
  
  # Proposed stick-breaking variables
  v_prop <- v
  v_prop[k] <- v_k
  
  # Construct stick-breaking weights (trust the product directly --
  # it is guaranteed >= 0 term-by-term; no residual subtraction needed)
  remaining <- cumprod(c(1, 1 - v_prop))
  w_prop <- c(v_prop, 1) * remaining
  
  # Normalising constant (safe to compute in full: a zero weight just
  # contributes zero here, not NaN, since this is not a log-sum)
  C_prop <- C_fun(w_prop, theta)
  
  # Only the v_k-dependent part of the offspring log-likelihood
  ll_sum <- n_comp[k]*log(v_k) + m_comp[k]*log(1 - v_k) - n_off*log(C_prop)
  
  # Add prior contribution and the log-Jacobian for the logit transformation
  # Prior contribution: (alpha - 1)*log(1 - v_k)
  # Jacobian contribution: log(v_k) + log(1 - v_k)
  return(ll_sum + (alpha - 1)*log(1 - v_k) + log(v_k) + log(1 - v_k))
}


#' Full conditional log-posterior for the mark-productivity slope beta
#'
#' Evaluates the (unnormalised) log full conditional of \code{beta}, for use
#' in a Metropolis-within-Gibbs update. Unlike \code{theta_k}/\code{v_k},
#' \code{beta} is already unconstrained (any sign is meaningful -- the data
#' determines whether productivity increases or decreases with mark), so
#' this uses a plain symmetric random-walk proposal on \code{beta} itself,
#' with no transform or Jacobian needed. 
#'
#' @param beta Numeric scalar; proposed value of the mark-productivity
#'   slope.
#' @param A Numeric scalar; current productivity scale.
#' @param marks Numeric vector of all event marks (every event, not just
#'   offspring -- every event can be a candidate parent).
#' @param sum_M_off Numeric scalar; sum of the parent's mark over every
#'   offspring event, i.e. \code{sum(marks[z[O_idx]])}.
#' @param mu_beta Numeric scalar; mean of the Normal prior on \code{beta}.
#' @param sd_beta Numeric scalar; standard deviation of the Normal prior on
#'   \code{beta}.
#'
#' @return Numeric scalar: the log full conditional density at \code{beta}.
#'
#' @keywords internal
logpost_beta <- function(beta, A, marks, sum_M_off, mu_beta, sd_beta) {
  # log-likelihood: beta * sum_{i in O} M_{z_i} - A * sum_j exp(beta * M_j)
  ll <- beta * sum_M_off - A * sum(exp(beta * marks))
  
  # Normal prior contribution (symmetric RW proposal -> no Jacobian needed)
  lp <- -(beta - mu_beta)^2 / (2 * sd_beta^2)
  
  return(ll + lp)
}


#' Run a single Metropolis-within-Gibbs marked-Hawkes MCMC chain
#'
#' Runs one chain of the sampler: at each iteration, resamples the latent
#' branching structure and component allocations
#' (\code{\link{update_zs}}), updates \code{lambda0}, \code{A}, \code{alpha},
#' \code{phi}, and \code{gamma} via their conjugate Gamma full conditionals,
#' and updates every \code{theta_k} and \code{v_k} via a
#' Metropolis-within-Gibbs sweep using \code{\link{logpost_theta_k}} and
#' \code{\link{logpost_v_k}}. This is an internal worker called once per
#' chain (typically in parallel) by \code{\link{run_mcmc}}; end users should
#' call \code{run_mcmc} directly rather than this function.
#'
#' @param times Numeric vector of event times, sorted ascending (already
#'   time-scaled if scaling is in use -- scaling itself is handled by
#'   \code{\link{run_mcmc}}, not here).
#' @param marks Numeric vector, same length as \code{times}; the mark of
#'   each event (not time-scaled -- marks are a separate physical quantity).
#' @param T_max Numeric scalar; observation window length (already
#'   time-scaled, consistent with \code{times}).
#' @param n_iter Integer; number of MCMC iterations to run.
#' @param init A list giving the initial state, with elements
#'   \code{lambda0}, \code{A}, \code{theta} (length K), \code{v} (length
#'   K-1), \code{alpha}, \code{phi}, and \code{gamma}.
#' @param prior_params A list of prior hyperparameters; see
#'   \code{\link{default_prior_params}} for the expected names.
#' @param proposal_sds A list of Metropolis random-walk proposal standard
#'   deviations; see \code{\link{default_proposal_sds}} for the expected
#'   names.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param progress Logical; if \code{TRUE}, display a text progress bar.
#'
#' @return A list with:
#'   \describe{
#'     \item{samples}{An \code{n_iter x (5 + K + (K-1))} matrix of posterior
#'       draws, with columns \code{lambda0}, \code{A}, \code{theta1...K},
#'       \code{v1...(K-1)}, \code{alpha}, \code{phi}, \code{gamma}.}
#'     \item{acceptance_rates}{Named numeric vector; overall Metropolis
#'       acceptance rate for the \code{theta} and \code{v} sweeps.}
#'     \item{n_immigrant}{Numeric vector of length \code{n_iter}; number of
#'       immigrant events sampled at each iteration.}
#'     \item{n_offspring}{Numeric vector of length \code{n_iter}; number of
#'       offspring events sampled at each iteration.}
#'   }
#'
#' @keywords internal
run_sampler <- function(times, marks, T_max, n_iter, init, prior_params, 
                        proposal_sds, kernel = c("step", "pwlin"), 
                        mark_productivity = c("linear", "exponential"), 
                        progress = TRUE) {
  
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)
  
  samples <- matrix(NA, n_iter, sum(lengths(init)))

  if (mark_productivity == "linear") {
    acceptance <- matrix(NA, n_iter, 2)
    proposals <- matrix(c(length(init$theta), length(init$v)), 
                        n_iter, 2, byrow = TRUE)
    colnames(samples) <- c("lambda0", "A", 
                           paste0("theta", 1:(length(init$theta))), 
                           paste0("v", 1:(length(init$v))),
                           "alpha", "phi", "gamma")
    colnames(acceptance) <- c("theta", "v")
    colnames(proposals) <- c("theta", "v")
  } else if (mark_productivity == "exponential") {
    acceptance <- matrix(NA, n_iter, 3)
    proposals <- matrix(c(length(init$theta), length(init$v), 1), 
                        n_iter, 3, byrow = TRUE)
    colnames(samples) <- c("lambda0", "A", "beta",
                           paste0("theta", 1:(length(init$theta))), 
                           paste0("v", 1:(length(init$v))),
                           "alpha", "phi", "gamma")
    colnames(acceptance) <- c("theta", "v", "beta")
    colnames(proposals) <- c("theta", "v", "beta")
  }
  n_immigrant <- numeric(n_iter)
  n_offspring <- numeric(n_iter)
  
  # initialise
  lambda0 <- init$lambda0
  A       <- init$A
  theta   <- init$theta
  v       <- init$v
  alpha   <- init$alpha
  phi     <- init$phi
  gamma   <- init$gamma
  if (mark_productivity == "exponential") beta <- init$beta 
  
  # some preliminary calculations
  N_T <- length(times)
  K <- length(init$theta)
  
  # Stick-break weights
  remaining <- cumprod(c(1, 1 - v))
  w <- c(v, 1)*remaining
  # Scaling required to make mu(t)/eta a probability density
  if (kernel == "step") {
    C_fun <- function(w, theta) sum(w*theta)
  } else if (kernel == "pwlin") {
    C_fun <- function(w, theta) 0.5*sum(w*theta^2)
  }
  C <- C_fun(w, theta)
  
  
  if (progress) {
    pb <- txtProgressBar(min = 0,      
                         max = n_iter, 
                         style = 3,    
                         width = 50,   
                         char = "=")   
  }
  
  for (iter in 1:n_iter) {
    
    # resample latent parameters
    zs <- update_zs(lambda0, A, 
                    ifelse(mark_productivity == "linear", NULL, beta),
                    w, C, theta, times, marks, kernel, mark_productivity)
    z <- zs$z
    s <- zs$s
    
    # calculate variables used in various posteriors
    I_idx <- which(z == 0)           # immigrant set
    O_idx <- which(z != 0)           # offspring set
    n_imm <- sum(z == 0)             # immigrant count
    n_off <- sum(z != 0)             # offspring count
    parent_idx <- unique(z[z != 0])  # indices of unique parent point
    n_comp <- integer(K)
    m_comp <- integer(K)
    for (k in 1:K) {
      n_comp[k] <- sum(s == k & !is.na(s))
      m_comp[k] <- sum(s > k & !is.na(s))
    }
    # sum of the parent's mark over every offspring event (needed for beta)
    sum_M_off <- sum(marks[z[O_idx]])
    
    n_immigrant[iter] = n_imm
    n_offspring[iter] = n_off
    
    # lambda0
    lambda0 <- rgamma(
      1,
      shape = prior_params$a_l0 + n_imm,
      rate = prior_params$b_l0 + T_max
    )
    
    #A 
    A <- rgamma(
      1, 
      shape = prior_params$a_A + n_off, 
      rate = prior_params$b_A + 
        ifelse(mark_productivity == "linear", sum(marks), sum(exp(beta * marks)))
    )

    # beta (mark-productivity slope; Metropolis, no transform needed since
    # beta is already on an unconstrained scale)
    if (mark_productivity == "exponential") {
      beta_prop <- rnorm(1, beta, proposal_sds$beta)
      lp_prop <- logpost_beta(beta_prop, A, marks, sum_M_off, 
                              prior_params$mu_beta, prior_params$sd_beta)
      lp_curr <- logpost_beta(beta, A, marks, sum_M_off, 
                              prior_params$mu_beta, prior_params$sd_beta)
      if (!is.finite(lp_curr)) {
        accept_beta <- is.finite(lp_prop)
      } else {
        accept_beta <- is.finite(lp_prop) && (log(runif(1)) < (lp_prop - lp_curr))
      }
      if (accept_beta) {
        beta <- beta_prop
      }
      acceptance[iter, "beta"] <- as.integer(accept_beta)
    }
      
    # Theta_k (Metropolis-within-Gibbs)
    accepted <- 0
    for (k_theta in 1:K) {
      current <- log(theta[k_theta])
      proposal <- rnorm(1, current, proposal_sds$theta_k)
      lp_prop <- logpost_theta_k(proposal, k_theta, w, theta, phi, 
                                 times, O_idx, z, s, C_fun, kernel)
      lp_curr <- logpost_theta_k(current, k_theta, w, theta, phi, 
                                 times, O_idx, z, s, C_fun, kernel)
      
      # Guard against non-finite log-posteriors so log_acc is never NaN:
      # if the current state is somehow non-finite, accept any finite proposal;
      # otherwise use the ordinary MH ratio.
      if (!is.finite(lp_curr)) {
        accept <- is.finite(lp_prop)
      } else {
        accept <- is.finite(lp_prop) && (log(runif(1)) < (lp_prop - lp_curr))
      }
      
      if (accept) {
        theta[k_theta] <- exp(proposal)
        accepted <- accepted + 1
      }
    }
    acceptance[iter, "theta"] <- accepted
    
    # v_k (Metropolis-within-Gibbs, systematic sweep over all K-1 sticks)
    accepted <- 0
    for (k_v in 1:(K - 1)) {
      current <- qlogis(v[k_v])
      proposal <- rnorm(1, current, proposal_sds$v_k)
      
      lp_prop <- logpost_v_k(proposal, k_v, v, theta, alpha, n_comp, m_comp, n_off, C_fun)
      lp_curr <- logpost_v_k(current,  k_v, v, theta, alpha, n_comp, m_comp, n_off, C_fun)
      
      if (!is.finite(lp_curr)) {
        accept <- is.finite(lp_prop)
      } else {
        accept <- is.finite(lp_prop) && (log(runif(1)) < (lp_prop - lp_curr))
      }
      
      if (accept) {
        v[k_v] <- plogis(proposal)
        accepted <- accepted + 1
      }
    }
    acceptance[iter, "v"] <- accepted
    
    # recalculate stick-break weights after sampling v_k
    remaining <- cumprod(c(1, 1 - v))
    w <- c(v, 1)*remaining
    # Scaled to make mu(t)/eta a probability density
    C <- C_fun(w, theta)
    
    # alpha 
    alpha <- rgamma(
      1, 
      shape = prior_params$a_alpha + K - 1, 
      rate  = prior_params$b_alpha - sum(log(1 - v))
    )
    
    # phi 
    phi <- rgamma(
      1, 
      shape = prior_params$a_phi + length(theta), 
      rate = prior_params$b_phi + sum(theta)
    )
    
    # gamma (rate of the Exponential mark distribution)
    gamma <- rgamma(
      1, 
      shape = prior_params$a_gamma + N_T, 
      rate = prior_params$b_gamma + sum(marks)
    )
    
    if (mark_productivity == "linear") {
      samples[iter, ] <- c(lambda0, A, theta, v, alpha, phi, gamma)
    } else if (mark_productivity == "exponential") {
      samples[iter, ] <- c(lambda0, A, beta, theta, v, alpha, phi, gamma)
    }
    
    if (progress) {
      setTxtProgressBar(pb, iter)
    }
  }
  
  if (progress) {
    close(pb)
  }
  
  acceptance_rates <- colSums(acceptance)/colSums(proposals)
  
  return(list(
    samples = samples, 
    acceptance_rates = acceptance_rates,
    n_immigrant = n_immigrant,
    n_offspring = n_offspring
  ))
}

#' Fit a marked Hawkes process with a Dirichlet process mixture kernel
#'
#' Runs \code{n_chains} independent Metropolis-within-Gibbs MCMC chains (in
#' parallel, via a \code{parallel::makeCluster} PSOCK cluster) for a marked
#' Hawkes process whose excitation kernel is a K-component truncated
#' stick-breaking (Dirichlet process) mixture of step or piecewise-linear
#' basis functions, and whose offspring productivity scales either linearly or 
#' exponentially with each parent's mark. This is the main user-facing entry 
#' point; \code{\link{run_sampler}} and friends are internal workers called once
#' per chain.
#'
#' Event times are optionally rescaled before fitting (\code{scale_time =
#' TRUE}, the default) so that the baseline intensity is roughly of order 1,
#' which helps the default proposal standard deviations behave sensibly
#' across data sets with very different time units or event rates. Marks are
#' never rescaled -- they are a separate physical quantity unrelated to
#' time. Posterior samples are transformed back to the original time scale
#' automatically before being returned (via \code{\link{unscale_chain}}), so
#' downstream code never needs to know whether scaling was used.
#'
#' @param times Numeric vector of event times, sorted ascending, on their
#'   original (unscaled) time unit.
#' @param marks Numeric vector, same length as \code{times}; the mark
#'   (e.g. magnitude) of each event.
#' @param T_max Numeric scalar; length of the observation window, on the
#'   same original time unit as \code{times}.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}. Partially
#'   matched via \code{match.arg}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}
#' @param K Integer; truncation level of the stick-breaking mixture (number
#'   of mixture components). Default 20.
#' @param n_chains Integer; number of independent MCMC chains to run in
#'   parallel. Default 4.
#' @param n_iter Integer; number of MCMC iterations per chain. Default 8000.
#' @param seed Integer; seed for \code{parallel::clusterSetRNGStream}, giving
#'   reproducible, independent RNG streams across chains. Default 123.
#' @param prior_params A list of prior hyperparameters; see
#'   \code{\link{default_prior_params}}. Defaults to
#'   \code{default_prior_params()}.
#' @param proposal_sds A list of Metropolis random-walk proposal standard
#'   deviations; see \code{\link{default_proposal_sds}}. Defaults to
#'   \code{default_proposal_sds()}.
#' @param scale_time Logical; if \code{TRUE} (the default), rescale time so
#'   the baseline intensity is roughly 1 before fitting, then transform
#'   posterior samples back to the original scale.
#' @param save_path Optional file path; if supplied, the fitted chains are
#'   saved to this path via \code{saveRDS} as a side effect.
#' @param progress Logical; if \code{TRUE}, display a text progress bar for
#'   each chain. Default \code{FALSE} (recommended when running many chains
#'   in parallel, since progress bars from parallel workers do not display
#'   in the calling session).
#'
#' @return A list with:
#'   \describe{
#'     \item{chains}{A list of length \code{n_chains}, each element being
#'       the list returned by \code{\link{run_sampler}} for that chain
#'       (already transformed back to the original time scale if
#'       \code{scale_time = TRUE}).}
#'     \item{time_scale}{Numeric scalar; the scaling factor applied to time
#'       (1 if \code{scale_time = FALSE}).}
#'     \item{settings}{A list recording the call's configuration
#'       (\code{T_max}, \code{n_events}, \code{kernel}, \code{K},
#'       \code{n_chains}, \code{n_iter}, \code{seed}, \code{scale_time}).}
#'     \item{prior_params, proposal_sds, init_list}{The prior
#'       hyperparameters, proposal standard deviations, and per-chain
#'       initial values used to fit the model.}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- run_mcmc(times, marks, T_max, kernel = "pwlin", 
#'                    mark_productivity = "linear", K = 10, n_chains = 4, 
#'                    n_iter = 8000)
#' result$chains[[1]]$acceptance_rates
#' }
#'
#' @export
run_mcmc <- function(
    times, 
    marks,
    T_max, 
    kernel = c("step", "pwlin"), 
    mark_productivity = c("linear", "exponential"),
    K = 20, 
    n_chains = 4,
    n_iter = 8000, 
    seed = 123, 
    prior_params = default_prior_params(),  
    proposal_sds = default_proposal_sds(), 
    scale_time = TRUE,
    save_path = NULL, 
    progress = FALSE
) {
  
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)
  
  if (scale_time) {
    # Approximate total event rate
    rate_total <- length(times)/T_max
    
    # Assume E(eta) approximately 0.5
    # Choose scaling so that the background intensity is approximately 1
    time_scale <- 2/rate_total
    
  } else {
    time_scale <- 1
  }
  
  times_scaled <- times / time_scale
  T_max_scaled <- T_max / time_scale
  
  # Random initialisation of chains
  if (mark_productivity == "linear") {
    init_list <- lapply(seq_len(n_chains), function(i) {
      list(
        lambda0 = runif(1, 0.25, 2),
        A       = runif(1, 0.2, 2),
        theta   = sort(exp(runif(K, log(0.05), log(5)))),
        v       = rbeta(K - 1, 1, 1),
        alpha   = runif(1, 0.5, 2),
        phi     = runif(1, 0.2, 1),
        gamma   = runif(1, 0.2, 2)
      )
    })
  } else if (mark_productivity == "exponential") {
    init_list <- lapply(seq_len(n_chains), function(i) {
      list(
        lambda0 = runif(1, 0.25, 2),
        A       = runif(1, 0.2, 2),
        beta    = rnorm(1, 0, 0.5),
        theta   = sort(exp(runif(K, log(0.05), log(5)))),
        v       = rbeta(K - 1, 1, 1),
        alpha   = runif(1, 0.5, 2),
        phi     = runif(1, 0.2, 1),
        gamma   = runif(1, 0.2, 2)
      )
    })
  }
  
  # Create parallel workers
  cl <- parallel::makeCluster(n_chains)
  
  # Ensure the cluster is stopped even if run_sampler errors
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Reproducible independent RNG streams
  parallel::clusterSetRNGStream(cl, iseed = seed)
  
  # Export required objects and functions
  parallel::clusterExport(
    cl,
    varlist = c(
      "times_scaled",
      "marks",
      "T_max_scaled", 
      "n_iter", 
      "init_list", 
      "prior_params",
      "proposal_sds", 
      "kernel", 
      "mark_productivity",
      "progress",
      "run_sampler", 
      "logpost_theta_k", 
      "logpost_v_k",
      "logpost_beta",
      "update_zs"
    ),
    envir = environment()
  )
  
  # Run MCMC chains
  chains <- parallel::parLapply(
    cl,
    seq_len(n_chains),
    function(i) {
      run_sampler(
        times = times_scaled,
        marks = marks,
        T_max = T_max_scaled,
        n_iter = n_iter,
        init = init_list[[i]],
        prior_params = prior_params,
        proposal_sds = proposal_sds,
        kernel = kernel,
        mark_productivity = mark_productivity,
        progress = progress
      )
    }
  )
  
  # Transform posterior samples back if scaled
  if (scale_time) {
    chains <- lapply(
      chains,
      unscale_chain,
      time_scale = time_scale
    )
  }
  
  # Save chains if requested
  if (!is.null(save_path)) {
    saveRDS(chains, save_path)
  }
  
  # Return results
  list(
    chains = chains,
    time_scale = time_scale,
    settings = list(
      T_max = T_max,
      n_events = length(times),
      kernel = kernel,
      K = K,
      n_chains = n_chains,
      n_iter = n_iter,
      seed = seed,
      scale_time = scale_time
    ),
    prior_params = prior_params,
    proposal_sds = proposal_sds,
    init_list = init_list
  )
}

#' Default prior hyperparameters
#'
#' @keywords internal
default_prior_params <- function() {
  list(
    a_l0 = 1,
    b_l0 = 1,
    a_alpha = 1,
    b_alpha = 1,
    a_phi = 1,
    b_phi = 1,
    a_A = 1,
    b_A = 2,
    mu_beta = 0,
    sd_beta = 1,
    a_gamma = 1, 
    b_gamma = 1
  )
}

#' Default Metropolis proposal standard deviations
#'
#' @keywords internal
default_proposal_sds <- function() {
  list(theta_k = 1, v_k = 1, beta = 0.2)
}

#' Transform a chain's posterior samples back to the original time scale
#'
#' Undoes the rescaling applied by \code{\link{run_mcmc}} when
#' \code{scale_time = TRUE}: \code{lambda0} and \code{phi} are rates (events
#' per unit time / inverse-time), so they are divided by \code{time_scale};
#' the kernel atoms \code{theta1, ..., thetaK} are themselves time
#' durations, so they are multiplied by \code{time_scale}. \code{A}, 
#' \code{beta}, \code{v1, ..., v(K-1)}, \code{alpha}, and \code{gamma} are 
#' dimensionless or relate to marks rather than time, and are left unchanged.
#'
#' @param chain A list as returned by \code{\link{run_sampler}}, containing
#'   a \code{$samples} matrix with columns \code{lambda0}, \code{A},
#'   \code{theta1...K}, \code{v1...(K-1)}, \code{alpha}, \code{phi},
#'   \code{gamma}.
#' @param time_scale Numeric scalar; the scaling factor originally applied
#'   to time (i.e. \code{times_scaled = times / time_scale}).
#'
#' @return The input \code{chain} list, with \code{chain$samples}
#'   transformed back to the original time scale.
#'
#' @keywords internal
unscale_chain <- function(chain, time_scale) {
  samples <- chain$samples
  
  # Baseline intensity
  samples[, "lambda0"] <- samples[, "lambda0"]/time_scale
  
  # Atom locations
  theta_cols <- grep("^theta", colnames(samples))
  samples[, theta_cols] <- samples[, theta_cols] * time_scale
  
  # Exponential rate parameter
  samples[, "phi"] <- samples[, "phi"] / time_scale
  
  chain$samples <- samples
  
  return(chain)
}

################################################################################
# Posterior kernel plotting
################################################################################

#' Plot the posterior marked-Hawkes excitation kernel
#'
#' Produces a posterior-mean-and-credible-band plot and/or a "spaghetti" plot
#' of individual posterior draws for the Hawkes excitation kernel with either 
#' the step or piecewise-linear basis function. 
#'
#' Posterior draws are combined from either the raw list of chain objects
#' returned by \code{run_mcmc()} or an already row-bound matrix of draws.
#' Computation is vectorised over posterior draws: the function loops over
#' grid points in \code{x_grid} (typically a few hundred) rather than over
#' posterior draws (which can number in the tens of thousands across chains
#' and iterations), since each grid-point iteration is then a single
#' vectorised \code{n_draws x K} matrix operation.
#'
#' @param chains Either (a) the list returned in \code{run_mcmc()$chains},
#'   i.e. a list of chain objects each containing a \code{$samples} matrix,
#'   or (b) an already row-bound matrix of posterior draws with columns
#'   \code{theta1, ..., thetaK}, \code{v1, ..., v(K-1)}, and \code{A}.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}. Partially
#'   matched via \code{match.arg}.
#' @param x_grid Optional numeric vector of x-values at which to evaluate the
#'   kernel. If \code{NULL} (the default), a regular grid of length
#'   \code{n_grid} from 0 to \code{x_max} is constructed.
#' @param x_max Upper limit for the automatically constructed x grid, used
#'   only when \code{x_grid} is \code{NULL}. Defaults to the 99th percentile
#'   of all posterior \code{theta} draws.
#' @param n_grid Integer; number of points in the automatically constructed
#'   x grid, used only when \code{x_grid} is \code{NULL}. Default 200.
#' @param ci_level Numeric in (0, 1); width of the posterior credible band
#'   shown in the summary panel, e.g. \code{0.90} gives the 5th/95th
#'   percentiles. Default 0.90.
#' @param n_spaghetti Integer; number of individual posterior draws to
#'   overlay in the spaghetti panel, subsampled (without replacement) from
#'   all available draws for readability. Default 500.
#' @param true_kernel Optional function of \code{x} giving a reference or
#'   true kernel to overlay in red on the relevant panel(s), e.g.
#'   \code{function(x) 0.6 * dunif(x)}. Default \code{NULL} (no overlay).
#' @param panel Character; which panel(s) to draw: \code{"both"} (default),
#'   \code{"summary"} (posterior mean + credible band only), or
#'   \code{"spaghetti"} (individual draws only).
#' @param seed Optional integer seed for the spaghetti-panel subsampling, for
#'   reproducible figures. Default \code{NULL} (no seed set).
#'
#' @return Invisibly, a list with components:
#'   \describe{
#'     \item{x_grid}{The grid of x-values the kernel was evaluated at.}
#'     \item{kernel_draws}{An \code{n_draws x length(x_grid)} matrix of the
#'       kernel evaluated at every posterior draw and grid point.}
#'     \item{kernel_mean}{Posterior mean of the kernel at each grid point.}
#'     \item{kernel_lower, kernel_upper}{Posterior credible band bounds at
#'       each grid point, per \code{ci_level}.}
#'   }
#'   The plot(s) are drawn as a side effect.
#'
#' @examples
#' \dontrun{
#' result <- run_mcmc(times, marks, T_max, kernel = "pwlin", K = 10)
#'
#' # directly from run_mcmc()'s output list of chains:
#' plot_hawkes_kernel(result$chains, kernel = "pwlin",
#'                     true_kernel = function(x) 0.6 * dunif(x))
#'
#' # or from an already row-bound matrix:
#' chains_matrix <- do.call(rbind, lapply(result$chains, `[[`, "samples"))
#' plot_hawkes_kernel(chains_matrix, kernel = "pwlin", n_spaghetti = 500,
#'                     seed = 1, true_kernel = function(x) 0.6 * dunif(x))
#'
#' # summary panel only, custom grid:
#' plot_hawkes_kernel(chains_matrix, kernel = "pwlin", panel = "summary",
#'                     x_grid = seq(0, 5, length.out = 300))
#' }
#'
#' @export
plot_hawkes_kernel <- function(
  chains,
  kernel = c("step", "pwlin"),
  x_grid = NULL,
  x_max = NULL,
  n_grid = 200,
  ci_level = 0.90,
  n_spaghetti = 500,
  true_kernel = NULL,
  panel = c("both", "summary", "spaghetti"),
  seed = NULL
) {
  kernel <- match.arg(kernel)
  panel  <- match.arg(panel)
  
  # ---- assemble posterior draws into a single matrix -----------------------
  if (is.matrix(chains)) {
    samples <- chains
  } else if (is.list(chains) && !is.null(chains[[1]]$samples)) {
    samples <- do.call(rbind, lapply(chains, function(ch) ch$samples))
  } else {
    stop("`chains` must be a matrix of posterior draws, or a list of chain ",
         "objects each containing a $samples matrix (e.g. run_mcmc()$chains).")
  }
  
  # ---- infer K and pull out the relevant columns ----------------------------
  theta_cols <- grep("^theta[0-9]+$", colnames(samples), value = TRUE)
  v_cols     <- grep("^v[0-9]+$",     colnames(samples), value = TRUE)
  if (length(theta_cols) == 0 || length(v_cols) == 0) {
    stop("Could not find theta*/v* columns in the samples matrix.")
  }
  # order by numeric suffix (column order in the matrix isn't guaranteed)
  theta_cols <- theta_cols[order(as.integer(sub("theta", "", theta_cols)))]
  v_cols     <- v_cols[order(as.integer(sub("v", "", v_cols)))]
  K <- length(theta_cols)
  stopifnot(length(v_cols) == K - 1)
  
  Theta   <- samples[, theta_cols, drop = FALSE]   # n_draws x K
  V       <- samples[, v_cols,     drop = FALSE]   # n_draws x (K-1)
  A       <- samples[, "A"]                        # n_draws
  n_draws <- nrow(samples)
  
  # ---- stick-breaking weights + kernel normalisation, vectorised over draws -
  one_minus_v <- 1 - V
  remaining   <- t(apply(cbind(1, one_minus_v), 1, cumprod))  # n_draws x K
  W_raw       <- cbind(V, 1) * remaining                       # rows sum to 1
  
  if (kernel == "step") {
    norm_const <- rowSums(W_raw * Theta)
  } else {
    norm_const <- 0.5 * rowSums(W_raw * Theta^2)
  }
  W <- W_raw / norm_const  
  
  # ---- build the x grid ------------------------------------------------------
  if (is.null(x_grid)) {
    if (is.null(x_max)) {
      x_max <- as.numeric(quantile(Theta, 0.99))
    }
    x_grid <- seq(0, x_max, length.out = n_grid)
  }
  
  # ---- evaluate the kernel at each grid point ------------------------------
  kernel_draws <- matrix(0, n_draws, length(x_grid))
  for (j in seq_along(x_grid)) {
    xj <- x_grid[j]
    if (kernel == "step") {
      basis <- (Theta > xj) * 1
    } else {
      basis <- pmax(Theta - xj, 0)
    }
    kernel_draws[, j] <- rowSums(W * basis)
  }
  
  # ---- posterior summaries ---------------------------------------------------
  tail_p       <- (1 - ci_level) / 2
  kernel_mean  <- colMeans(kernel_draws)
  kernel_lower <- apply(kernel_draws, 2, quantile, probs = tail_p)
  kernel_upper <- apply(kernel_draws, 2, quantile, probs = 1 - tail_p)
  
  # ---- plotting ---------------------------------------------------------------
  if (panel == "both") {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)
  }
  
  if (panel %in% c("both", "summary")) {
    plot(x_grid, kernel_mean, type = "l", lwd = 2,
         ylim = range(kernel_lower, kernel_upper),
         xlab = "x", ylab = "f(x)",
         main = paste0("Posterior kernel (", kernel, "), ",
                       round(ci_level * 100), "% CI"))
    lines(x_grid, kernel_lower, lty = 2)
    lines(x_grid, kernel_upper, lty = 2)
    if (!is.null(true_kernel)) {
      lines(x_grid, true_kernel(x_grid), lwd = 2, col = "red")
    }
  }
  
  if (panel %in% c("both", "spaghetti")) {
    if (!is.null(seed)) set.seed(seed)
    n_show   <- min(n_spaghetti, n_draws)
    draw_idx <- sample.int(n_draws, n_show)
    matplot(x_grid, t(kernel_draws[draw_idx, , drop = FALSE]),
            type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
            xlab = "x", ylab = "f(x)",
            main = paste0("Posterior kernel draws (", kernel, ")"))
    if (!is.null(true_kernel)) {
      lines(x_grid, true_kernel(x_grid), lwd = 3, col = "red")
    }
  }
  
  invisible(list(
    x_grid       = x_grid,
    kernel_draws = kernel_draws,
    kernel_mean  = kernel_mean,
    kernel_lower = kernel_lower,
    kernel_upper = kernel_upper
  ))
}

################################################################################
# Compensator and residual diagnostics for the marked Hawkes process
# linear mark-productivity: A * mark, or 
# exponential mark-productivity: A * exp(beta * mark)
################################################################################

#' Ground-process compensator for the marked Hawkes process
#'
#' Computes \eqn{\Lambda(t) = \int_0^t \lambda_g^*(s)\,ds}, the integrated
#' ground intensity up to time \code{t}, for use in residual/goodness-of-fit
#' diagnostics via the random time change theorem.
#'
#' @param t Numeric scalar; time at which to evaluate the compensator.
#' @param times Numeric vector of event times, sorted ascending.
#' @param marks Numeric vector, same length as \code{times}; event marks.
#' @param params Numeric vector \code{c(lambda0, A, theta_1, ..., theta_K, 
#'   v_1, ..., v_{K-1})} -- e.g. one row/summary of the posterior
#'   samples matrix restricted to these columns. Should be named (as
#'   \code{colMeans()} or a single row of a samples matrix would be) so
#'   \code{theta}/\code{v} are matched by name rather than position.
#' @param K Integer; number of mixture components.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}
#' @param n_before Optional integer; number of events strictly before
#'   \code{t}. If \code{NULL} (the default), computed as
#'   \code{sum(times < t)}. Passing this in directly (as
#'   \code{\link{hawkes_resid}} does, since it always evaluates at
#'   \code{t = times[i]} and so knows \code{n_before = i - 1} without a
#'   scan) avoids an O(n) scan per call -- useful when \code{compensator()}
#'   is called many times, e.g. across many posterior draws.
#'
#' @return Numeric scalar, \eqn{\Lambda(t)}.
#'
#' @export
compensator <- function(t, times, marks, params, K, kernel = c("step", "pwlin"),
                        mark_productivity = c("linear", "exponential"), 
                        n_before = NULL) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)
  
  # Number of events strictly before t. times is sorted ascending, so this
  # is a direct count rather than max(which(...)), which returns -Inf (not
  # 0) when nothing precedes t and then breaks any 1:(idx - 1) slicing.
  if (is.null(n_before)) {
    n_before <- sum(times < t)
  }
  
  # Extract by name 
  stopifnot(!is.null(names(params)))

  lambda_0 <- params[["lambda0"]]
  A <- params[["A"]]
  if (mark_productivity == "exponential") beta <- params[["beta"]]
  
  theta_names <- grep("^theta[0-9]+$", names(params), value = TRUE)
  theta_names <- theta_names[order(as.integer(sub("theta", "", theta_names)))]
  theta <- unname(params[theta_names])
  
  v_names <- grep("^v[0-9]+$", names(params), value = TRUE)
  v_names <- v_names[order(as.integer(sub("v", "", v_names)))]
  v <- unname(params[v_names])
  
  stopifnot(length(theta) == K, length(v) == K - 1)
  
  remaining <- cumprod(c(1, 1 - v))
  w <- c(v, 1) * remaining
  
  if (kernel == "step") {
    C_fun <- function(w, theta) sum(w * theta)
  } else if (kernel == "pwlin") {
    C_fun <- function(w, theta) 0.5 * sum(w * theta^2)
  }
  C <- C_fun(w, theta)
  
  # No events precede t: only the baseline term contributes
  if (n_before == 0) {
    return(lambda_0 * t)
  }
  
  parent_times <- times[1:n_before]
  parent_marks <- marks[1:n_before]
  dt <- t - parent_times   # elapsed time since each parent (length n_before)
  
  # min(dt_i, theta_k) for every (parent, component) pair -- n_before x K
  dt_theta_min <- outer(dt, theta, pmin)
  
  if (kernel == "step") {
    # F_step(dt_i) = sum_k w_k * min(dt_i, theta_k)
    F_vals <- as.numeric(dt_theta_min %*% w)
  } else if (kernel == "pwlin") {
    # F_pwlin(dt_i) = sum_k w_k * (theta_k*min(dt_i,theta_k) - 0.5*min(dt_i,theta_k)^2)
    theta_mat <- matrix(theta, nrow = n_before, ncol = K, byrow = TRUE)
    F_vals <- as.numeric((theta_mat * dt_theta_min - 0.5 * dt_theta_min^2) %*% w)
  }
  
  # productivity of each parent: (A * mark) or (A * exp(beta * mark))
  if (mark_productivity == "linear") {
    productivity <- A * parent_marks
  } else if (mark_productivity == "exponential") {
    productivity <- A * exp(beta * parent_marks)
  }
  
  lambda_0 * t + (1 / C) * sum(productivity * F_vals)
}


#' Compensator values at each observed event time
#'
#' Evaluates \code{\link{compensator}} at every event time \code{times[i]},
#' returning the sequence \eqn{\Lambda(T_1), \ldots, \Lambda(T_n)}. Under a
#' correctly specified model, these are the arrival times of a unit-rate
#' Poisson process, so \code{diff(c(0, tau))} should look like iid
#' \eqn{\text{Exp}(1)} residuals -- this function returns the raw
#' compensator values, not the differenced residuals, so take
#' \code{diff()} yourself if that's what you need next.
#'
#' @param times Numeric vector of event times, sorted ascending.
#' @param marks Numeric vector, same length as \code{times}; event marks.
#' @param params Numeric vector, as in \code{\link{compensator}}.
#' @param K Integer; number of mixture components.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#'
#' @return Numeric vector, same length as \code{times}: \eqn{\Lambda(T_i)}
#'   for each \code{i}.
#'
#' @export
hawkes_resid <- function(times, marks, params, K, kernel = c("step", "pwlin"), 
                         mark_productivity = c("linear", "exponential")) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)
  
  tau <- sapply(seq_along(times), function(i) {
    # n_before = i - 1 always, since t = times[i] and times is sorted --
    # passing it in skips compensator()'s O(n) sum(times < t) scan, turning
    # an O(n^2) call sequence into O(n) scans total.
    compensator(times[i], times, marks, params, K, kernel, 
                mark_productivity, n_before = i - 1)
  })
  
  return(tau)
}


#' Posterior-predictive compensator, averaged across posterior draws
#'
#' Rather than plugging \code{colMeans()} of the posterior samples into a
#' single \code{\link{compensator}} evaluation, this evaluates the
#' compensator separately for each of a subsample of posterior draws, then
#' averages the resulting \eqn{\Lambda(T_i)} sequences pointwise across
#' draws. This matters specifically for mixture-component parameters
#' (\code{theta_k}, \code{v_k}): they are only weakly identified
#' individually (exchangeable across the K components, prone to label
#' switching across MCMC iterations), so averaging them directly in
#' parameter space can reconstruct a kernel that doesn't correspond to any
#' actual posterior draw. Averaging the compensator *values* instead
#' (function-space averaging) is the same principle already used by
#' \code{\link{plot_hawkes_kernel}}, and is robust to this.
#'
#' @param times Numeric vector of event times, sorted ascending.
#' @param marks Numeric vector, same length as \code{times}; event marks.
#' @param chains Either the list returned in \code{run_mcmc()$chains}, or
#'   an already row-bound matrix of posterior draws (as in
#'   \code{\link{plot_hawkes_kernel}}).
#' @param K Integer; number of mixture components.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param n_draws Integer; number of posterior draws to subsample (without
#'   replacement) for the average. Default 200 -- each draw costs a full
#'   \code{\link{hawkes_resid}} evaluation (O(n) compensator calls, each
#'   O(K) after the \code{n_before} fix above), so this trades precision
#'   for runtime on large event catalogs.
#' @param seed Optional integer seed for reproducible subsampling.
#'
#' @return A list with:
#'   \describe{
#'     \item{tau_mean}{Numeric vector, same length as \code{times}: the
#'       posterior-predictive-averaged \eqn{\Lambda(T_i)} at each event
#'       time. Feed \code{diff(c(0, tau_mean))} into your existing Q-Q /
#'       K-S code exactly as you would \code{\link{hawkes_resid}}'s
#'       output.}
#'     \item{tau_draws}{\code{n_draws x length(times)} matrix; each row is
#'       one draw's own (self-consistent, not averaged) \eqn{\Lambda(T_i)}
#'       sequence -- useful for a spaghetti-style residual plot showing
#'       posterior uncertainty in the diagnostic itself.}
#'     \item{draw_idx}{Integer vector; which rows of the combined samples
#'       matrix were used.}
#'   }
#'
#' @export
hawkes_resid_posterior <- function(
  times, 
  marks, 
  chains, 
  K,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential"),
  n_draws = 200, 
  seed = NULL
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)
  
  if (is.matrix(chains)) {
    samples <- chains
  } else if (is.list(chains) && !is.null(chains[[1]]$samples)) {
    samples <- do.call(rbind, lapply(chains, function(ch) ch$samples))
  } else {
    stop("`chains` must be a matrix of posterior draws, or a list of chain ",
         "objects each containing a $samples matrix (e.g. run_mcmc()$chains).")
  }
  
  n_total <- nrow(samples)
  if (!is.null(seed)) set.seed(seed)
  n_draws <- min(n_draws, n_total)
  draw_idx <- sample.int(n_total, n_draws)
  
  n_events <- length(times)
  tau_draws <- matrix(NA_real_, n_draws, n_events)
  for (d in seq_len(n_draws)) {
    tau_draws[d, ] <- hawkes_resid(times, marks, samples[draw_idx[d], ], K, 
                                   kernel, mark_productivity)
  }
  
  list(
    tau_mean  = colMeans(tau_draws),
    tau_draws = tau_draws,
    draw_idx  = draw_idx
  )
}

#' Plot posterior-predictive residual Q-Q diagnostic
#'
#' Plots a residual Q-Q diagnostic from the output of
#' \code{\link{hawkes_resid_posterior}}, showing how posterior parameter
#' uncertainty (evaluated on the fixed, observed event sequence) propagates
#' into the residual diagnostic itself. Three display styles are available: 
#' a shaded credible ribbon (\code{"ribbon"}, a dashed-line credible band 
#' (\code{"band"}), or raw per-draw spaghetti lines (\code{"spaghetti"} -- shows 
#' the full draw-by-draw variation.
#'
#' Note that the resulting band/ribbon reflects posterior *parameter*
#' uncertainty propagated onto the observed catalog (\code{times}/
#' \code{marks} are fixed; only the parameter draw varies).
#'
#' @param tau_draws \code{n_draws x n_events} matrix, as returned by
#'   \code{hawkes_resid_posterior()$tau_draws}.
#' @param tau_mean Optional numeric vector of length \code{n_events} (e.g.
#'   \code{hawkes_resid_posterior()$tau_mean}); if supplied, its residuals
#'   are overlaid as points on top.
#' @param style Character; \code{"ribbon"} (default), \code{"band"}, or
#'   \code{"spaghetti"}.
#' @param band_level Numeric in (0, 1); credible band/ribbon width for
#'   \code{style \%in\% c("ribbon", "band")}, e.g. \code{0.90} gives the
#'   5th/95th percentiles. Default 0.90.
#' @param ribbon_col Colour for the shaded ribbon when
#'   \code{style = "ribbon"}. Default a light grey, chosen to stay legible
#'   in greyscale printing.
#' @param main Character; plot title.
#'
#' @return Invisibly, a list with \code{theoretical_q} (the shared
#'   theoretical quantiles) and \code{gaps_matrix} (the \code{n_events x
#'   n_draws} matrix of sorted per-draw gaps used for plotting).
#'
#' @export
plot_hawkes_resid_qq <- function(tau_draws, tau_mean = NULL,
                                 style = c("ribbon", "band", "spaghetti"),
                                 band_level = 0.90,
                                 ribbon_col = grey(0.85),
                                 main = "Posterior Predictive Residual Q-Q Plot") {
  style <- match.arg(style)
  n_events <- ncol(tau_draws)
  theoretical_q <- qexp(ppoints(n_events))
  
  # sorted gaps for every draw -> n_events x n_draws matrix (one column per draw)
  gaps_matrix <- apply(tau_draws, 1, function(tau_d) sort(diff(c(0, tau_d))))
  
  gaps_mean <- if (!is.null(tau_mean)) sort(diff(c(0, tau_mean))) else NULL
  
  max_val <- max(theoretical_q, gaps_matrix, gaps_mean) * 1.05
  
  par(pty = "s")
  plot(NA, xlim = c(0, max_val), ylim = c(0, max_val),
       xlab = "Theoretical Exp(1) Quantiles", ylab = "Sample Residual Quantiles",
       main = main)
  
  if (style %in% c("ribbon", "band")) {
    tail_p <- (1 - band_level) / 2
    gaps_lower <- apply(gaps_matrix, 1, quantile, probs = tail_p)
    gaps_upper <- apply(gaps_matrix, 1, quantile, probs = 1 - tail_p)
    
    if (style == "ribbon") {
      # shaded polygon, drawn before the reference line so the line sits on top
      polygon(c(theoretical_q, rev(theoretical_q)),
              c(gaps_lower, rev(gaps_upper)),
              col = ribbon_col, border = NA)
    } else {
      lines(theoretical_q, gaps_lower, lty = 2)
      lines(theoretical_q, gaps_upper, lty = 2)
    }
  } else {
    matlines(theoretical_q, gaps_matrix, lty = 1, col = rgb(0, 0, 0, 0.05))
  }
  
  abline(0, 1, col = "slategrey", lty = 2, lwd = 2)
  
  if (!is.null(gaps_mean)) {
    points(theoretical_q, gaps_mean, pch = 20, cex = 0.5)
  }
  
  invisible(list(theoretical_q = theoretical_q, gaps_matrix = gaps_matrix))
}