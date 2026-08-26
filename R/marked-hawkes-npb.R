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

#' Sample latent parent and kernel-component allocations
#'
#' Samples the parent event and mixture component for each event using
#' a categorical draw over all valid immigrant and offspring assignments.
#'
#' @param lambda0 Numeric; baseline intensity.
#' @param A Numeric; mark productivity scale.
#' @param beta Numeric; exponential mark-productivity slope, if required.
#' @param w Numeric vector; length-K mixture weights.
#' @param C Numeric; kernel normalising constant.
#' @param theta Numeric vector; length-K kernel atom locations.
#' @param times Numeric vector; sorted event times.
#' @param marks Numeric vector; event marks.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#'
#' @return List with integer vectors \code{z} (parent indices, with 0 for
#'   immigrants) and \code{s} (kernel-component allocations).
#'
#' @keywords internal
update_zs <- function(
  lambda0,
  A,
  beta = NULL,
  w,
  C,
  theta,
  times,
  marks,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential")
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  if (mark_productivity == "exponential" && is.null(beta)) {
    stop("Need 'beta' param when mark_productivity = 'exponential'")
  }

  n <- length(times)
  z <- integer(n)
  s <- integer(n)
  s[] <- NA_integer_

  if (n <= 1) {
    return(list(z = z, s = s))
  }

  for (i in 2:n) {
    # Time since each possible parent event
    dt <- times[i] - times[1:(i - 1)]

    # Matrix of basis-function values:
    # rows = possible parents j = i - 1
    # columns = mixture components k
    if (kernel == "step") {
      g <- dt < matrix(theta, nrow = i - 1, ncol = length(theta), byrow = TRUE)
    } else if (kernel == "pwlin") {
      g <- pmax(
        outer(dt, theta, FUN = function(dt, theta) {
          theta - dt
        }),
        0
      )
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
    comp_idx <- active_flat_idx[, 2]

    # Offspring probabilities
    if (mark_productivity == "linear") {
      offspring_probs <- A *
        marks[parent_idx] *
        w[comp_idx] *
        g[active_flat_idx] /
        C
    } else if (mark_productivity == "exponential") {
      offspring_probs <- A *
        exp(beta * marks[parent_idx]) *
        w[comp_idx] *
        g[active_flat_idx] /
        C
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

#' Log full conditional for a kernel atom
#'
#' Evaluates the log full conditional of \code{theta_k} on the log scale,
#' including the Jacobian of the log transformation.
#'
#' @param log_theta_k Numeric; proposed log atom location.
#' @param k Integer; component index.
#' @param w Numeric vector; mixture weights.
#' @param theta Numeric vector; current atom locations.
#' @param phi Numeric; Gamma prior rate for \code{theta_k}.
#' @param times Numeric vector; event times.
#' @param O_idx Integer vector; offspring-event indices.
#' @param z Integer vector; parent indices.
#' @param s Integer vector; kernel-component allocations.
#' @param C_fun Function; computes the kernel normalising constant from
#'   \code{w} and \code{theta}.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#'
#' @return Numeric; log full conditional density, or \code{-Inf} for an
#'   invalid proposal.
#'
#' @keywords internal
logpost_theta_k <- function(
  log_theta_k,
  k,
  w,
  theta,
  phi,
  times,
  O_idx,
  z,
  s,
  C_fun,
  kernel = c("step", "pwlin")
) {
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
  own <- (s_O == k)
  dt_own <- dt_O[own]

  if (kernel == "step") {
    if (any(theta_k <= dt_own)) {
      return(-Inf)
    }
    log_g_own_sum <- 0 # log(1) for every satisfied step-indicator
  } else if (kernel == "pwlin") {
    g_own <- theta_k - dt_own
    if (any(g_own <= 0)) {
      return(-Inf)
    }
    log_g_own_sum <- sum(log(g_own))
  }

  # sum(log(w[s_O])) does not depend on theta_k and cannot be -Inf
  # (stick-breaking weights are always > 0), so it's safe to include in full
  ll_sum <- sum(log(w[s_O])) + log_g_own_sum - length(O_idx) * log(C_prop)

  # Return posterior + log-Jacobian
  return(ll_sum - phi * theta_k + log_theta_k)
}


#' Log full conditional for a stick-breaking variable
#'
#' Evaluates the log full conditional of \code{v_k} on the logit scale,
#' including the Beta prior and transformation Jacobian.
#'
#' @param logit_v_k Numeric; proposed logit stick-breaking value.
#' @param k Integer; component index, from 1 to K-1.
#' @param v Numeric vector; current stick-breaking variables.
#' @param theta Numeric vector; kernel atom locations.
#' @param alpha Numeric; Beta(1, alpha) concentration parameter.
#' @param n_comp Integer vector; offspring counts by component.
#' @param m_comp Integer vector; offspring counts allocated to components
#'   with index greater than each component.
#' @param n_off Integer; total number of offspring events.
#' @param C_fun Function; computes the kernel normalising constant.
#'
#' @return Numeric; log full conditional density at \code{logit_v_k}.
#'
#' @keywords internal
logpost_v_k <- function(
  logit_v_k,
  k,
  v,
  theta,
  alpha,
  n_comp,
  m_comp,
  n_off,
  C_fun
) {
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
  ll_sum <- n_comp[k] *
    log(v_k) +
    m_comp[k] * log(1 - v_k) -
    n_off * log(C_prop)

  # Add prior contribution and the log-Jacobian for the logit transformation
  # Prior contribution: (alpha - 1)*log(1 - v_k)
  # Jacobian contribution: log(v_k) + log(1 - v_k)
  return(ll_sum + (alpha - 1) * log(1 - v_k) + log(v_k) + log(1 - v_k))
}


#' Log full conditional for the exponential mark-productivity slope
#'
#' Evaluates the log full conditional of \code{beta}.
#'
#' @param beta Numeric; proposed mark-productivity slope.
#' @param A Numeric; productivity scale.
#' @param marks Numeric vector; event marks.
#' @param sum_M_off Numeric; sum of parent marks over offspring events.
#' @param mu_beta Numeric; Normal prior mean.
#' @param sd_beta Numeric; Normal prior standard deviation.
#'
#' @return Numeric; log full conditional density at \code{beta}.
#'
#' @keywords internal
logpost_beta <- function(beta, A, marks, sum_M_off, mu_beta, sd_beta) {
  # log-likelihood: beta * sum_{i in O} M_{z_i} - A * sum_j exp(beta * M_j)
  ll <- beta * sum_M_off - A * sum(exp(beta * marks))

  # Normal prior contribution (symmetric RW proposal -> no Jacobian needed)
  lp <- -(beta - mu_beta)^2 / (2 * sd_beta^2)

  return(ll + lp)
}


#' Run one marked-Hawkes MCMC chain
#'
#' Runs a Metropolis-within-Gibbs sampler for one marked Hawkes process chain.
#' This is an internal worker called by \code{\link{run_mcmc}}.
#'
#' @param times Numeric vector; sorted event times on the fitting time scale.
#' @param marks Numeric vector; event marks.
#' @param T_max Numeric; observation-window length on the fitting time scale.
#' @param n_iter Integer; number of MCMC iterations.
#' @param init List; initial values for \code{lambda0}, \code{A}, \code{beta}
#'   (if exponential), \code{theta}, \code{v}, \code{alpha}, \code{phi}, and
#'   \code{gamma}.
#' @param prior_params List; prior hyperparameters, as returned by
#'   \code{\link{default_prior_params}}.
#' @param proposal_sds List; Metropolis proposal standard deviations, as
#'   returned by \code{\link{default_proposal_sds}}.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#' @param progress Logical; whether to display a progress bar.
#' @param adapt Logical; whether to adapt proposal standard deviations.
#' @param adapt_start Integer; first iteration eligible for adaptation.
#' @param adapt_end Integer; last iteration eligible for adaptation.
#' @param adapt_interval Integer; number of iterations between adaptations.
#' @param target_accept Numeric; target Metropolis acceptance rate.
#'
#' @return List containing:
#'   \describe{
#'     \item{samples}{Matrix of posterior draws.}
#'     \item{acceptance_rates}{Named acceptance rates for each Metropolis-updated parameter.}
#'     \item{n_immigrant}{Immigrant count at each iteration.}
#'     \item{n_offspring}{Offspring count at each iteration.}
#'     \item{tuned_proposal_sds}{Final adapted proposal standard deviations.}
#'   }
#'
#' @keywords internal
#' @export
run_sampler <- function(
  times,
  marks,
  T_max,
  n_iter,
  init,
  prior_params,
  proposal_sds,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential"),
  progress = TRUE,
  adapt = TRUE,
  adapt_start = 1,
  adapt_end = floor(n_iter / 2),
  adapt_interval = 100,
  target_accept = 0.3
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  samples <- matrix(NA, n_iter, sum(lengths(init)))

  if (mark_productivity == "linear") {
    acceptance <- matrix(NA, n_iter, length(init$theta) + length(init$v))
    # proposals <- matrix(c(length(init$theta), length(init$v)),
    #                     n_iter, 2, byrow = TRUE)
    colnames(samples) <- c(
      "lambda0",
      "A",
      paste0("theta", 1:(length(init$theta))),
      paste0("v", 1:(length(init$v))),
      "alpha",
      "phi",
      "gamma"
    )
    colnames(acceptance) <- c(
      paste0("theta", 1:length(init$theta)),
      paste0("v", 1:length(init$v))
    )
    # colnames(proposals) <- c("theta", "v")
  } else if (mark_productivity == "exponential") {
    acceptance <- matrix(NA, n_iter, length(init$theta) + length(init$v) + 1)
    # proposals <- matrix(c(length(init$theta), length(init$v), 1),
    #                     n_iter, 3, byrow = TRUE)
    colnames(samples) <- c(
      "lambda0",
      "A",
      "beta",
      paste0("theta", 1:(length(init$theta))),
      paste0("v", 1:(length(init$v))),
      "alpha",
      "phi",
      "gamma"
    )
    colnames(acceptance) <- c(
      paste0("theta", 1:length(init$theta)),
      paste0("v", 1:length(init$v)),
      "beta"
    )
    # colnames(proposals) <- c("theta", "v", "beta")
  }
  n_immigrant <- numeric(n_iter)
  n_offspring <- numeric(n_iter)

  # initialise
  lambda0 <- init$lambda0
  A <- init$A
  theta <- init$theta
  v <- init$v
  alpha <- init$alpha
  phi <- init$phi
  gamma <- init$gamma
  if (mark_productivity == "exponential") {
    beta <- init$beta
  }

  # some preliminary calculations
  N_T <- length(times)
  K <- length(init$theta)

  # Proposal SDs can be supplied either as a single value or as a separate value
  # for each parameter.
  theta_proposal_sds <- rep(proposal_sds$theta_k, length.out = K)
  v_proposal_sds <- rep(proposal_sds$v_k, length.out = K - 1)
  if (mark_productivity == "exponential") {
    beta_proposal_sd <- proposal_sds$beta
  }

  # Stick-break weights
  remaining <- cumprod(c(1, 1 - v))
  w <- c(v, 1) * remaining
  # Scaling required to make mu(t)/eta a probability density
  if (kernel == "step") {
    C_fun <- function(w, theta) sum(w * theta)
  } else if (kernel == "pwlin") {
    C_fun <- function(w, theta) 0.5 * sum(w * theta^2)
  }
  C <- C_fun(w, theta)

  if (progress) {
    pb <- txtProgressBar(
      min = 0,
      max = n_iter,
      style = 3,
      width = 50,
      char = "="
    )
  }

  # Per-parameter acceptance counters used for adaptive tuning
  batch_accept_theta <- numeric(K)
  batch_proposals_theta <- numeric(K)
  batch_accept_v <- numeric(K - 1)
  batch_proposals_v <- numeric(K - 1)
  if (mark_productivity == "exponential") {
    batch_accept_beta <- 0
    batch_proposals_beta <- 0
  }

  for (iter in 1:n_iter) {
    # resample latent parameters
    if (mark_productivity == "linear") {
      zs <- update_zs(
        lambda0,
        A,
        NULL,
        w,
        C,
        theta,
        times,
        marks,
        kernel,
        mark_productivity
      )
    } else {
      zs <- update_zs(
        lambda0,
        A,
        beta,
        w,
        C,
        theta,
        times,
        marks,
        kernel,
        mark_productivity
      )
    }
    z <- zs$z
    s <- zs$s

    # calculate variables used in various posteriors
    I_idx <- which(z == 0) # immigrant set
    O_idx <- which(z != 0) # offspring set
    n_imm <- sum(z == 0) # immigrant count
    n_off <- sum(z != 0) # offspring count
    parent_idx <- unique(z[z != 0]) # indices of unique parent point
    n_comp <- integer(K)
    m_comp <- integer(K)
    for (k in 1:K) {
      n_comp[k] <- sum(s == k & !is.na(s))
      m_comp[k] <- sum(s > k & !is.na(s))
    }
    # sum of the parent's mark over every offspring event (needed for beta)
    sum_M_off <- sum(marks[z[O_idx]])

    n_immigrant[iter] <- n_imm
    n_offspring[iter] <- n_off

    # lambda0
    lambda0 <- rgamma(
      1,
      shape = prior_params$a_l0 + n_imm,
      rate = prior_params$b_l0 + T_max
    )

    # A
    if (mark_productivity == "linear") {
      A_rate <- prior_params$b_A + sum(marks)
    } else {
      A_rate <- prior_params$b_A + sum(exp(beta * marks))
    }

    A <- rgamma(1, shape = prior_params$a_A + n_off, rate = A_rate)

    # beta (mark-productivity slope; Metropolis, no transform needed since
    # beta is already on an unconstrained scale)
    if (mark_productivity == "exponential") {
      beta_prop <- rnorm(1, beta, beta_proposal_sd)
      lp_prop <- logpost_beta(
        beta_prop,
        A,
        marks,
        sum_M_off,
        prior_params$mu_beta,
        prior_params$sd_beta
      )
      lp_curr <- logpost_beta(
        beta,
        A,
        marks,
        sum_M_off,
        prior_params$mu_beta,
        prior_params$sd_beta
      )
      if (!is.finite(lp_curr)) {
        accept_beta <- is.finite(lp_prop)
      } else {
        accept_beta <- is.finite(lp_prop) &&
          (log(runif(1)) < (lp_prop - lp_curr))
      }
      if (accept_beta) {
        beta <- beta_prop
      }
      acceptance[iter, "beta"] <- as.integer(accept_beta)

      # Adaptation bookkeeping
      if (adapt && iter >= adapt_start && iter <= adapt_end) {
        batch_proposals_beta <- batch_proposals_beta + 1
        batch_accept_beta <- batch_accept_beta + as.integer(accept_beta)
      }
    }

    # Theta_k (Metropolis-within-Gibbs)
    for (k_theta in 1:K) {
      current <- log(theta[k_theta])
      proposal <- rnorm(1, current, theta_proposal_sds[k_theta])
      lp_prop <- logpost_theta_k(
        proposal,
        k_theta,
        w,
        theta,
        phi,
        times,
        O_idx,
        z,
        s,
        C_fun,
        kernel
      )
      lp_curr <- logpost_theta_k(
        current,
        k_theta,
        w,
        theta,
        phi,
        times,
        O_idx,
        z,
        s,
        C_fun,
        kernel
      )

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
      }
      acceptance[iter, paste0("theta", k_theta)] <- as.integer(accept)

      # Adaptation bookkeeping
      if (adapt && iter >= adapt_start && iter <= adapt_end) {
        batch_proposals_theta[k_theta] <- batch_proposals_theta[k_theta] + 1
        batch_accept_theta[k_theta] <- batch_accept_theta[k_theta] +
          as.integer(accept)
      }
    }

    # v_k (Metropolis-within-Gibbs, systematic sweep over all K-1 sticks)
    for (k_v in 1:(K - 1)) {
      current <- qlogis(v[k_v])
      proposal <- rnorm(1, current, v_proposal_sds[k_v])

      lp_prop <- logpost_v_k(
        proposal,
        k_v,
        v,
        theta,
        alpha,
        n_comp,
        m_comp,
        n_off,
        C_fun
      )
      lp_curr <- logpost_v_k(
        current,
        k_v,
        v,
        theta,
        alpha,
        n_comp,
        m_comp,
        n_off,
        C_fun
      )

      if (!is.finite(lp_curr)) {
        accept <- is.finite(lp_prop)
      } else {
        accept <- is.finite(lp_prop) && (log(runif(1)) < (lp_prop - lp_curr))
      }

      if (accept) {
        v[k_v] <- plogis(proposal)
      }
      acceptance[iter, paste0("v", k_v)] <- as.integer(accept)

      # Adaptation bookkeeping
      if (adapt && iter >= adapt_start && iter <= adapt_end) {
        batch_proposals_v[k_v] <- batch_proposals_v[k_v] + 1
        batch_accept_v[k_v] <- batch_accept_v[k_v] + as.integer(accept)
      }
    }

    # Adaptive tuning during burn-in
    if (
      adapt &&
        iter >= adapt_start &&
        iter <= adapt_end &&
        iter %% adapt_interval == 0
    ) {
      # ---- theta_k ----
      theta_rates <- batch_accept_theta / pmax(batch_proposals_theta, 1)

      for (k in 1:K) {
        if (batch_proposals_theta[k] > 0) {
          # smooth log-scale adaptation
          theta_proposal_sds[k] <- theta_proposal_sds[k] *
            exp(1 * (theta_rates[k] - target_accept))
          # apply some lower and upper bounds
          theta_proposal_sds[k] <- min(max(theta_proposal_sds[k], 0.01), 5)
        }
      }

      # ---- v_k ----
      v_rates <- batch_accept_v / pmax(batch_proposals_v, 1)
      for (k in 1:(K - 1)) {
        if (batch_proposals_v[k] > 0) {
          # smooth log-scale adaptation
          v_proposal_sds[k] <- v_proposal_sds[k] *
            exp(1 * (v_rates[k] - target_accept))
          # apply some lower and upper bounds
          v_proposal_sds[k] <- min(max(v_proposal_sds[k], 0.01), 3)
        }
      }

      # ---- beta ----
      if (mark_productivity == "exponential") {
        if (batch_proposals_beta > 0) {
          beta_rate <- batch_accept_beta / batch_proposals_beta
          # smooth log-scale adaptation
          beta_proposal_sd <- beta_proposal_sd *
            exp(1 * (beta_rate - target_accept))
          # apply some lower and upper bounds
          beta_proposal_sd <- min(max(beta_proposal_sd, 0.01), 3)
        }
      }

      # Reset batch counters
      batch_accept_theta <- numeric(K)
      batch_proposals_theta <- numeric(K)
      batch_accept_v <- numeric(K - 1)
      batch_proposals_v <- numeric(K - 1)
      if (mark_productivity == "exponential") {
        batch_accept_beta <- 0
        batch_proposals_beta <- 0
      }
    }

    # recalculate stick-break weights after sampling v_k
    remaining <- cumprod(c(1, 1 - v))
    w <- c(v, 1) * remaining
    # Scaled to make mu(t)/eta a probability density
    C <- C_fun(w, theta)

    # alpha
    alpha <- rgamma(
      1,
      shape = prior_params$a_alpha + K - 1,
      rate = prior_params$b_alpha - sum(log(1 - v))
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

  # post adaptation acceptance
  acceptance_rates <- colMeans(acceptance[adapt_end:n_iter, ])

  return(list(
    samples = samples,
    acceptance_rates = acceptance_rates,
    n_immigrant = n_immigrant,
    n_offspring = n_offspring,
    tuned_proposal_sds = list(
      theta_k = theta_proposal_sds,
      v_k = v_proposal_sds,
      beta = if (mark_productivity == "exponential") {
        beta_proposal_sd
      } else {
        NULL
      }
    )
  ))
}

#' Fit a marked Hawkes process with a DP mixture excitation kernel
#'
#' Runs independent parallel MCMC chains for a marked Hawkes process with a
#' truncated stick-breaking mixture of step or piecewise-linear kernels and
#' linear or exponential mark-dependent productivity.
#'
#' Event times can be rescaled for fitting and are automatically transformed
#' back to the original time scale in the returned posterior samples.
#'
#' @param times Numeric vector; sorted event times.
#' @param marks Numeric vector; event marks.
#' @param T_max Numeric; observation-window length.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#' @param K Integer; number of mixture components. Default 20.
#' @param n_chains Integer; number of parallel MCMC chains. Default 4.
#' @param n_iter Integer; iterations per chain. Default 8000.
#' @param seed Integer; random-number seed. Default 123.
#' @param prior_params List; prior hyperparameters. Defaults to
#'   \code{\link{default_prior_params}}().
#' @param proposal_sds List; Metropolis proposal standard deviations. Defaults
#'   to \code{\link{default_proposal_sds}}().
#' @param scale_time Logical; whether to rescale time before fitting. Default
#'   \code{TRUE}.
#' @param save_path Character or \code{NULL}; optional path for saving chains
#'   with \code{saveRDS}.
#' @param progress Logical; whether workers display progress bars. Default
#'   \code{FALSE}.
#' @param adapt Logical; whether to adapt proposal standard deviations.
#' @param adapt_start Integer; first iteration eligible for adaptation.
#' @param adapt_end Integer; last iteration eligible for adaptation.
#' @param adapt_interval Integer; iterations between adaptations.
#' @param target_accept Numeric; target Metropolis acceptance rate. Default 0.30.
#'
#' @return List containing the fitted \code{chains}, \code{time_scale},
#'   \code{settings}, \code{prior_params}, \code{proposal_sds}, and
#'   \code{init_list}.
#'
#' @examples
#' \dontrun{
#' fit <- run_mcmc(times, marks, T_max, kernel = "pwlin",
#'                 mark_productivity = "linear")
#' fit$chains[[1]]$acceptance_rates
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
  progress = FALSE,
  adapt = TRUE,
  adapt_start = 1,
  adapt_end = floor(n_iter / 2),
  adapt_interval = 100,
  target_accept = 0.30
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  if (scale_time) {
    # Approximate total event rate
    rate_total <- length(times) / T_max

    # Assume E(eta) approximately 0.5
    # Choose scaling so that the background intensity is approximately 1
    time_scale <- 2 / rate_total
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
        A = runif(1, 0.2, 2),
        theta = sort(exp(runif(K, log(0.05), log(5)))),
        v = rbeta(K - 1, 1, 1),
        alpha = runif(1, 0.5, 2),
        phi = runif(1, 0.2, 1),
        gamma = runif(1, 0.2, 2)
      )
    })
  } else if (mark_productivity == "exponential") {
    init_list <- lapply(seq_len(n_chains), function(i) {
      list(
        lambda0 = runif(1, 0.25, 2),
        A = runif(1, 0.2, 2),
        beta = rnorm(1, 0, 0.5),
        theta = sort(exp(runif(K, log(0.05), log(5)))),
        v = rbeta(K - 1, 1, 1),
        alpha = runif(1, 0.5, 2),
        phi = runif(1, 0.2, 1),
        gamma = runif(1, 0.2, 2)
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
      "adapt",
      "adapt_start",
      "adapt_end",
      "adapt_interval",
      "target_accept",
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
        progress = progress,
        adapt = adapt,
        adapt_start = adapt_start,
        adapt_end = adapt_end,
        adapt_interval = adapt_interval,
        target_accept = target_accept
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
#' Returns the default Gamma and Normal prior hyperparameters used by the
#' marked-Hawkes sampler.
#'
#' @return Named list of prior hyperparameters for \code{lambda0}, \code{alpha},
#'   \code{phi}, \code{A}, \code{beta}, and \code{gamma}.
#'
#' @keywords internal
#' @export
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

#' Default MCMC proposal standard deviations
#'
#' Returns the default random-walk proposal standard deviations for
#' \code{theta}, \code{v}, and \code{beta}.
#'
#' @return Named list with elements \code{theta_k}, \code{v_k}, and
#'   \code{beta}.
#'
#' @keywords internal
#' @export
default_proposal_sds <- function() {
  list(theta_k = 1, v_k = 1, beta = 0.2)
}

#' Transform posterior samples back to the original time scale
#'
#' Reverses the time scaling applied by \code{\link{run_mcmc}} to rates and
#' kernel atom locations.
#'
#' @param chain List; fitted chain returned by \code{\link{run_sampler}}.
#' @param time_scale Numeric; time scaling factor used during fitting.
#'
#' @return The input chain with time-dependent posterior parameters restored
#'   to their original scale.
#'
#' @keywords internal
unscale_chain <- function(chain, time_scale) {
  samples <- chain$samples

  # Baseline intensity
  samples[, "lambda0"] <- samples[, "lambda0"] / time_scale

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

#' Plot the posterior excitation kernel
#'
#' Plots the posterior mean and credible band, posterior draws, or both for
#' a step or piecewise-linear Hawkes excitation kernel.
#'
#' @param chains Posterior draws as either a list of chains from
#'   \code{\link{run_mcmc}} or a row-bound samples matrix.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param x_grid Numeric vector; evaluation grid. If \code{NULL}, constructed
#'   from \code{x_max} and \code{n_grid}.
#' @param x_max Numeric; upper grid limit when \code{x_grid = NULL}. If
#'   \code{NULL}, chosen from posterior kernel tail mass.
#' @param n_grid Integer; number of automatically generated grid points.
#' @param ci_level Numeric in (0, 1); posterior credible-band level.
#' @param n_spaghetti Integer; number of posterior draws to plot.
#' @param true_kernel Function or \code{NULL}; optional reference kernel to
#'   overlay.
#' @param panel Character; \code{"both"}, \code{"summary"}, or
#'   \code{"spaghetti"}.
#' @param legend Logical; whether to display the plot legend.
#' @param seed Integer or \code{NULL}; seed for posterior-draw subsampling.
#'
#' @return Invisibly returns a list containing \code{x_grid},
#'   \code{kernel_draws}, \code{kernel_mean}, \code{kernel_lower}, and
#'   \code{kernel_upper}.
#'
#' @examples
#' \dontrun{
#' fit <- run_mcmc(times, marks, T_max, kernel = "pwlin")
#' plot_hawkes_kernel(fit$chains, kernel = "pwlin")
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
  legend = TRUE,
  seed = NULL
) {
  kernel <- match.arg(kernel)
  panel <- match.arg(panel)

  # ---- assemble posterior draws into a single matrix -----------------------
  if (is.matrix(chains)) {
    samples <- chains
  } else if (is.list(chains) && !is.null(chains[[1]]$samples)) {
    samples <- do.call(rbind, lapply(chains, function(ch) ch$samples))
  } else {
    stop(
      "`chains` must be a matrix of posterior draws, or a list of chain ",
      "objects each containing a $samples matrix (e.g. run_mcmc()$chains)."
    )
  }

  # ---- infer K and pull out the relevant columns ----------------------------
  theta_cols <- grep("^theta[0-9]+$", colnames(samples), value = TRUE)
  v_cols <- grep("^v[0-9]+$", colnames(samples), value = TRUE)
  if (length(theta_cols) == 0 || length(v_cols) == 0) {
    stop("Could not find theta*/v* columns in the samples matrix.")
  }
  # order by numeric suffix (column order in the matrix isn't guaranteed)
  theta_cols <- theta_cols[order(as.integer(sub("theta", "", theta_cols)))]
  v_cols <- v_cols[order(as.integer(sub("v", "", v_cols)))]
  K <- length(theta_cols)
  stopifnot(length(v_cols) == K - 1)

  Theta <- samples[, theta_cols, drop = FALSE] # n_draws x K
  V <- samples[, v_cols, drop = FALSE] # n_draws x (K-1)
  A <- samples[, "A"] # n_draws
  n_draws <- nrow(samples)

  # ---- stick-breaking weights + kernel normalisation, vectorised over draws -
  one_minus_v <- 1 - V
  remaining <- t(apply(cbind(1, one_minus_v), 1, cumprod)) # n_draws x K
  W_raw <- cbind(V, 1) * remaining # rows sum to 1

  if (kernel == "step") {
    norm_const <- rowSums(W_raw * Theta)
  } else {
    norm_const <- 0.5 * rowSums(W_raw * Theta^2)
  }
  W <- W_raw / norm_const

  # ---- build the x grid ------------------------------------------------------
  # if no x_grid given, use (0, x_max) with n_grid points (200 default)
  if (is.null(x_grid)) {
    # if no x_max, calculate the point at which no more than 5% of kernel mass
    # is left out, (ci_level)% of the time. i.e. for default ci_level of 0.90,
    # the plot range will contain at least 95% of the mass 90% of the time
    if (is.null(x_max)) {
      tail_tol <- 0.05 # 5% - fraction of kernel mass willing to leave out

      # Find x_max for each posterior draw such that
      # P_kernel(T > x_max) <= tail_tol.
      x_max_draw <- numeric(n_draws)
      for (i in seq_len(n_draws)) {
        w_i <- W_raw[i, ]
        theta_i <- Theta[i, ]

        # Kernel normalising constant
        if (kernel == "step") {
          C_i <- sum(w_i * theta_i)
        } else {
          C_i <- 0.5 * sum(w_i * theta_i^2)
        }

        # Tail mass of the normalised kernel at x
        tail_mass <- function(x) {
          if (kernel == "step") {
            # Integral_x^infinity (theta_k - x)_+
            # weighted by w_k
            sum(w_i * pmax(theta_i - x, 0)) / C_i
          } else {
            # For the piecewise-linear kernel:
            # integral_x^theta (theta - t) dt
            # = 0.5 * (theta - x)^2
            sum(w_i * 0.5 * pmax(theta_i - x, 0)^2) / C_i
          }
        }

        # At x = 0 the tail mass should be 1.
        # At x = max(theta) it is 0.
        x_max_draw[i] <- uniroot(
          function(x) tail_mass(x) - tail_tol,
          interval = c(0, max(theta_i)),
          tol = 1e-8
        )$root
      }

      # Use a high posterior quantile so that the plotting range
      # covers essentially all posterior draws.
      x_max <- as.numeric(quantile(x_max_draw, ci_level))
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
  tail_p <- (1 - ci_level) / 2
  kernel_mean <- colMeans(kernel_draws)
  kernel_lower <- apply(kernel_draws, 2, quantile, probs = tail_p)
  kernel_upper <- apply(kernel_draws, 2, quantile, probs = 1 - tail_p)

  # ---- plotting ---------------------------------------------------------------
  if (panel == "both") {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par), add = TRUE)
  }

  if (!is.null(true_kernel)) {
    y_lim <- range(
      min(kernel_lower, true_kernel(x_grid)),
      max(kernel_upper, true_kernel(x_grid))
    )
  } else {
    y_lim <- range(kernel_lower, kernel_upper)
  }

  if (panel %in% c("both", "summary")) {
    plot(
      x_grid,
      kernel_mean,
      type = "l",
      lwd = 2,
      ylim = y_lim,
      xlab = "x",
      ylab = "f(x)",
      main = paste0(
        "Posterior kernel (",
        kernel,
        "), ",
        round(ci_level * 100),
        "% CI"
      )
    )
    lines(x_grid, kernel_lower, lty = 2)
    lines(x_grid, kernel_upper, lty = 2)
    if (!is.null(true_kernel)) {
      lines(x_grid, true_kernel(x_grid), lwd = 2, col = "red")
    }

    if (legend) {
      legend_labels <- c(
        "Mean",
        paste0(round(ci_level * 100), "% Credible Interval")
      )
      legend_col <- c("black", "black")
      legend_lty <- c(1, 2)
      legend_lwd <- c(2, 1)

      if (!is.null(true_kernel)) {
        legend_labels <- c(legend_labels, "Truth")
        legend_col <- c(legend_col, "red")
        legend_lty <- c(legend_lty, 1)
        legend_lwd <- c(legend_lwd, 2)
      }

      legend(
        x = "topright",
        legend = legend_labels,
        col = legend_col,
        lty = legend_lty,
        lwd = legend_lwd,
        bty = "n"
      )
    }
  }

  if (panel %in% c("both", "spaghetti")) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    n_show <- min(n_spaghetti, n_draws)
    draw_idx <- sample.int(n_draws, n_show)

    if (!is.null(true_kernel)) {
      y_lim <- range(
        min(kernel_draws[draw_idx, ], true_kernel(x_grid)),
        max(kernel_draws[draw_idx, ], true_kernel(x_grid))
      )
    } else {
      y_lim <- range(
        min(kernel_draws[draw_idx, ]),
        max(kernel_draws[draw_idx, ])
      )
    }

    matplot(
      x_grid,
      t(kernel_draws[draw_idx, , drop = FALSE]),
      type = "l",
      lty = 1,
      col = rgb(0, 0, 0, 0.1),
      xlab = "x",
      ylab = "f(x)",
      ylim = y_lim,
      main = paste0("Posterior kernel draws (", kernel, ")")
    )
    if (!is.null(true_kernel)) {
      lines(x_grid, true_kernel(x_grid), lwd = 3, col = "red")
    }

    if (legend) {
      legend_labels <- c("Posterior Draw")
      legend_col <- c(rgb(0, 0, 0, 0.1))
      legend_lty <- c(1)
      legend_lwd <- c(1)

      if (!is.null(true_kernel)) {
        legend_labels <- c(legend_labels, "Truth")
        legend_col <- c(legend_col, "red")
        legend_lty <- c(legend_lty, 1)
        legend_lwd <- c(legend_lwd, 3)
      }

      legend(
        x = "topright",
        legend = legend_labels,
        col = legend_col,
        lty = legend_lty,
        lwd = legend_lwd,
        bty = "n"
      )
    }
  }

  invisible(list(
    x_grid = x_grid,
    kernel_draws = kernel_draws,
    kernel_mean = kernel_mean,
    kernel_lower = kernel_lower,
    kernel_upper = kernel_upper
  ))
}

################################################################################
# Compensator and residual diagnostics for the marked Hawkes process
# linear mark-productivity: A * mark, or
# exponential mark-productivity: A * exp(beta * mark)
################################################################################

#' Compute the marked-Hawkes ground-process compensator
#'
#' Computes the integrated ground intensity \eqn{\Lambda(t)} up to time
#' \code{t} for residual and goodness-of-fit diagnostics.
#'
#' @param t Numeric; evaluation time.
#' @param times Numeric vector; sorted event times.
#' @param marks Numeric vector; event marks.
#' @param params Named numeric vector; posterior parameter values.
#' @param K Integer; number of mixture components.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#' @param n_before Integer or \code{NULL}; number of events before \code{t}.
#'   If \code{NULL}, computed from \code{times}.
#'
#' @return Numeric scalar; \eqn{\Lambda(t)}.
#'
#' @export
compensator <- function(
  t,
  times,
  marks,
  params,
  K,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential"),
  n_before = NULL
) {
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
  if (mark_productivity == "exponential") {
    beta <- params[["beta"]]
  }

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
  dt <- t - parent_times # elapsed time since each parent (length n_before)

  # min(dt_i, theta_k) for every (parent, component) pair -- n_before x K
  dt_theta_min <- outer(dt, theta, pmin)

  if (kernel == "step") {
    # F_step(dt_i) = sum_k w_k * min(dt_i, theta_k)
    F_vals <- as.numeric(dt_theta_min %*% w)
  } else if (kernel == "pwlin") {
    # F_pwlin(dt_i) = sum_k w_k * (theta_k*min(dt_i,theta_k) - 0.5*min(dt_i,theta_k)^2)
    theta_mat <- matrix(theta, nrow = n_before, ncol = K, byrow = TRUE)
    F_vals <- as.numeric(
      (theta_mat * dt_theta_min - 0.5 * dt_theta_min^2) %*% w
    )
  }

  # productivity of each parent: (A * mark) or (A * exp(beta * mark))
  if (mark_productivity == "linear") {
    productivity <- A * parent_marks
  } else if (mark_productivity == "exponential") {
    productivity <- A * exp(beta * parent_marks)
  }

  lambda_0 * t + (1 / C) * sum(productivity * F_vals)
}


#' Compute compensators at observed event times
#'
#' Evaluates the ground-process compensator at each observed event time.
#'
#' @param times Numeric vector; sorted event times.
#' @param marks Numeric vector; event marks.
#' @param params Named numeric vector; posterior parameter values.
#' @param K Integer; number of mixture components.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#'
#' @return Numeric vector of \eqn{\Lambda(T_i)} values.
#'
#' @export
hawkes_resid <- function(
  times,
  marks,
  params,
  K,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential")
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  tau <- sapply(seq_along(times), function(i) {
    # n_before = i - 1 always, since t = times[i] and times is sorted --
    # passing it in skips compensator()'s O(n) sum(times < t) scan, turning
    # an O(n^2) call sequence into O(n) scans total.
    compensator(
      times[i],
      times,
      marks,
      params,
      K,
      kernel,
      mark_productivity,
      n_before = i - 1
    )
  })

  return(tau)
}


#' Compute posterior-averaged Hawkes compensators
#'
#' Evaluates the compensator for a subsample of posterior draws and averages
#' the resulting values at each observed event time.
#'
#' @param times Numeric vector; sorted event times.
#' @param marks Numeric vector; event marks.
#' @param chains Posterior draws as a chain list or samples matrix.
#' @param K Integer; number of mixture components.
#' @param kernel Character; \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Character; \code{"linear"} or \code{"exponential"}.
#' @param n_draws Integer; number of posterior draws to use. Default 200.
#' @param seed Integer or \code{NULL}; seed for draw subsampling.
#'
#' @return List containing \code{tau_mean}, \code{tau_draws}, and
#'   \code{draw_idx}.
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
    stop(
      "`chains` must be a matrix of posterior draws, or a list of chain ",
      "objects each containing a $samples matrix (e.g. run_mcmc()$chains)."
    )
  }

  n_total <- nrow(samples)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_draws <- min(n_draws, n_total)
  draw_idx <- sample.int(n_total, n_draws)

  n_events <- length(times)
  tau_draws <- matrix(NA_real_, n_draws, n_events)
  for (d in seq_len(n_draws)) {
    tau_draws[d, ] <- hawkes_resid(
      times,
      marks,
      samples[draw_idx[d], ],
      K,
      kernel,
      mark_productivity
    )
  }

  list(
    tau_mean = colMeans(tau_draws),
    tau_draws = tau_draws,
    draw_idx = draw_idx
  )
}

#' Plot the posterior-predictive residual Q-Q diagnostic
#'
#' Plots Exp(1) Q-Q diagnostics from posterior compensator draws using a
#' credible ribbon, credible band, or draw-by-draw spaghetti plot.
#'
#' @param tau_draws Numeric matrix; compensator values from posterior draws.
#' @param tau_mean Numeric vector or \code{NULL}; optional posterior-mean
#'   compensator values to overlay.
#' @param style Character; \code{"ribbon"}, \code{"band"}, or \code{"spaghetti"}.
#' @param band_level Numeric in (0, 1); credible-band level. Default 0.90.
#' @param ribbon_col Character; ribbon colour when \code{style = "ribbon"}.
#' @param main Character; plot title.
#'
#' @return Invisibly returns \code{theoretical_q} and \code{gaps_matrix}.
#'
#' @export
plot_hawkes_resid_qq <- function(
  tau_draws,
  tau_mean = NULL,
  style = c("ribbon", "band", "spaghetti"),
  band_level = 0.90,
  ribbon_col = grey(0.85),
  main = "Posterior Predictive Residual Q-Q Plot"
) {
  style <- match.arg(style)
  n_events <- ncol(tau_draws)
  theoretical_q <- qexp(ppoints(n_events))

  # sorted gaps for every draw -> n_events x n_draws matrix (one column per draw)
  gaps_matrix <- apply(tau_draws, 1, function(tau_d) sort(diff(c(0, tau_d))))

  gaps_mean <- if (!is.null(tau_mean)) sort(diff(c(0, tau_mean))) else NULL

  max_val <- max(theoretical_q, gaps_matrix, gaps_mean) * 1.05

  par(pty = "s")
  plot(
    NA,
    xlim = c(0, max_val),
    ylim = c(0, max_val),
    xlab = "Theoretical Exp(1) Quantiles",
    ylab = "Sample Residual Quantiles",
    main = main
  )

  if (style %in% c("ribbon", "band")) {
    tail_p <- (1 - band_level) / 2
    gaps_lower <- apply(gaps_matrix, 1, quantile, probs = tail_p)
    gaps_upper <- apply(gaps_matrix, 1, quantile, probs = 1 - tail_p)

    if (style == "ribbon") {
      # shaded polygon, drawn before the reference line so the line sits on top
      polygon(
        c(theoretical_q, rev(theoretical_q)),
        c(gaps_lower, rev(gaps_upper)),
        col = ribbon_col,
        border = NA
      )
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

################################################################################
# Model selection using WAIC
################################################################################

#' Compute WAIC for Hawkes Process MCMC Sampler
#'
#' @param samples Matrix of posterior draws with named columns.
#' @param times Sorted numeric vector of event timestamps.
#' @param marks Numeric vector of event marks.
#' @param T_max Total observation time horizon.
#' @param kernel Kernel type: \code{"step"} or \code{"pwlin"}.
#' @param mark_productivity Productivity law: \code{"linear"} or \code{"exponential"}.
#' @param burn_in Proportion of MCMC iterations to discard (0 to 1).
#'
#' @return List containing total \code{waic}, \code{lppd}, \code{p_waic},
#'   standard error \code{se}, and a \code{pointwise} data frame.
#' @export
compute_waic_hawkes <- function(
  samples,
  times,
  marks,
  T_max,
  kernel = c("step", "pwlin"),
  mark_productivity = c("linear", "exponential"),
  burn_in = 0.5
) {
  kernel <- match.arg(kernel)
  mark_productivity <- match.arg(mark_productivity)

  # Discard burn-in samples
  n_iter <- nrow(samples)
  post_burn <- seq(floor(n_iter * burn_in) + 1, n_iter)
  S_eff <- length(post_burn)
  n <- length(times)

  # Setup clean matrix for pointwise log-likelihoods
  log_lik_matrix <- matrix(0, nrow = S_eff, ncol = n)

  # Event time differences
  dt_mat <- outer(times, times, "-")
  dt_mat[upper.tri(dt_mat, diag = TRUE)] <- NA # Only keep past events (j < i)

  message("Computing pointwise log-likelihood matrix...")

  for (s_idx in seq_len(S_eff)) {
    # Extract parameter vector for this MCMC draw
    draw <- samples[post_burn[s_idx], ]

    # Extract parameter values from parameter vector
    l0 <- draw["lambda0"]
    A <- draw["A"]

    # Broadcasting calculations using our pre-computed dt_mat matrix
    # Reconstruct stick-breaking weights for the mixture kernel
    theta_vals <- draw[grep("^theta", names(draw))]
    v_vals <- draw[grep("^v", names(draw))]
    K <- length(theta_vals)
    remaining <- cumprod(c(1, 1 - v_vals))
    w <- c(v_vals, 1) * remaining

    # Kernel Normalising Constants C
    if (kernel == "step") {
      C_const <- sum(w * theta_vals)
    } else if (kernel == "pwlin") {
      C_const <- 0.5 * sum(w * theta_vals^2)
    }

    # Compute productivity term eta(M) for all historical items
    if (mark_productivity == "linear") {
      eta <- A * marks
    } else {
      beta <- draw["beta"]
      eta <- A * exp(beta * marks)
    }

    # Evaluate Kernel density f(t - t_j) for all pairs
    f_mat <- matrix(0, nrow = n, ncol = n)

    if (kernel == "step") {
      for (k in 1:K) {
        # indicator basis I(theta > \tau)
        f_mat <- f_mat +
          (w[k] * (dt_mat < theta_vals[k] & dt_mat > 0)) / C_const
      }
    } else if (kernel == "pwlin") {
      for (k in 1:K) {
        # max(0, theta - dt) evaluation
        diff_val <- theta_vals[k] - dt_mat
        # Replace future points (NA) and expired points (< 0) with 0
        diff_val[is.na(diff_val) | diff_val < 0] <- 0
        f_mat <- f_mat + (w[k] * diff_val) / C_const
      }
    }
    f_mat[is.na(f_mat)] <- 0

    # Multiply rows by eta and sum them up using matrix multiplication.
    trigger_intensities <- as.vector(f_mat %*% eta)

    # Total instantaneous conditional intensity at point i
    lambda_i <- l0 + trigger_intensities

    # Assign pointwise log-likelihood allocations
    log_lik_matrix[s_idx, ] <- log(lambda_i) - ((l0 * T_max) / n) - eta
  }

  message("Calculating final Information Criteria...")

  # Stable log-sum-exp per column
  log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }

  lppd_i <- apply(log_lik_matrix, 2, function(col) {
    log_sum_exp(col) - log(S_eff)
  })
  lppd <- sum(lppd_i)

  # Effective parameters penalty (Sample variance per column)
  p_waic_i <- apply(log_lik_matrix, 2, var)
  p_waic <- sum(p_waic_i)

  waic_value <- -2 * (lppd - p_waic)

  # Compute Standard Error of the difference
  se_waic <- sqrt(n * var(-2 * (lppd_i - p_waic_i)))

  return(list(
    waic = waic_value,
    lppd = lppd,
    p_waic = p_waic,
    se = se_waic,
    pointwise = data.frame(lppd = lppd_i, p_waic = p_waic_i)
  ))
}
