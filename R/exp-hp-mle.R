#' @title Exponential Hawkes Process Intensity
#' @description Calculates the conditional intensity \eqn{\lambda^*(t)} at time 
#'   \eqn{t}.
#' @param t Numeric value. The time at which to evaluate the intensity.
#' @param H_t Numeric vector of event timestamps occurring before \eqn{t}.
#' @param theta Numeric vector of parameters \eqn{c(\lambda, \alpha, \beta)}.
#' @export
exp_hp_intensity <- function(t, H_t, theta) {
  # Standard calculation for a single point t
  # unpack parameters
  lambda <- theta[1]
  alpha <- theta[2]
  beta <- theta[3]
  
  lambda + alpha*sum(exp(-beta*(t - H_t)))
}

#' @keywords internal
.exp_hp_intensity_at_events <- function(theta, H_t) {
  # Internal recursive O(n) trick to get intensity at all t_i simultaneously
  lambda <- theta[1]
  alpha <- theta[2]
  beta <- theta[3]
  
  n <- length(H_t)
  # Recursive computation of the summation term A(i)
  A <- numeric(n)
  if (n > 1) {
    for (i in 2:n) {
      A[i] <- exp(-beta*(H_t[i] - H_t[i-1])) * (1 + A[i-1])
    }
  }
  lambda + alpha*A
}

#' @title Exponential Hawkes Compensator
#' @description Calculates the integrated intensity \eqn{\Lambda(t)} over the 
#'   interval \eqn{[0, t]}.
#' @param theta Numeric vector of parameters \eqn{c(\lambda, \alpha, \beta)}.
#' @param H_t Numeric vector of event timestamps.
#' @param t Numeric value. The time at which to evaluate the compensator 
#'   (typically T_max).
#' @export
exp_hp_compensator <- function(theta, H_t, t) {
  lambda <- theta[1]
  alpha <- theta[2]
  beta <- theta[3]
  
  if (t <= 0) return(0)
  lambda*t + (alpha/beta)*sum(1 - exp(-beta*(t - H_t)))
}

#' @title Hawkes Log-Likelihood (Exponential)
#' @description Calculates the log-likelihood of a Hawkes process given a 
#'   history.
#' @param theta Numeric vector of parameters \eqn{c(\lambda, \alpha, \beta)}.
#' @param H_t Numeric vector of event timestamps.
#' @param T_max Numeric value. The end of the observation window.
#' @param lik_method Character string, either "fast" (recursive) or "slow" 
#'   (direct).
#' @export
exp_hp_loglik <- function(theta, H_t, T_max, lik_method = c("fast", "slow")) {
  # Switch between .exp_hp_loglik_fast and .exp_hp_loglik_slow  
  # Match the method argument (defaults to "fast")
  lik_method <- match.arg(lik_method)
  
  # validation
  if (any(theta <= 0)) return(-Inf)
  
  if (lik_method == "fast") {
    return(.exp_hp_loglik_fast(theta, H_t, T_max))
  } else {
    return(.exp_hp_loglik_slow(theta, H_t, T_max))
  }
}

#' @keywords internal
.exp_hp_loglik_fast <- function(theta, H_t, T_max) {
  # Implementation of the recursive O(n) log-likelihood
  intensities <- .exp_hp_intensity_at_events(theta = theta, H_t = H_t)
  log_intensities <- log(intensities)
  sum(log_intensities) - exp_hp_compensator(theta = theta, t = T_max, H_t = H_t)
}

#' @keywords internal
.exp_hp_loglik_slow <- function(theta, H_t, T_max) {
  # Implementation of the direct O(n^2) log-likelihood
  log_intensities <- sapply(H_t, function(t_i) {
    past_events <- H_t[H_t < t_i] # Find history by comparing values
    log(exp_hp_intensity(theta = theta, t = t_i, H_t = past_events))
  })
  sum(log_intensities) - exp_hp_compensator(theta = theta, t = T_max, H_t = H_t)
}

#' @keywords internal
.exp_hp_negloglik <- function(theta, 
                              H_t, 
                              T_max, 
                              lik_method = c("fast", "slow")) {
  # Wrapper returning negative log-likelihood for optim()
  lik_method <- match.arg(lik_method) 
  # Return the negative of the public function
  -exp_hp_loglik(theta, H_t, T_max, lik_method = lik_method)
}

#' @title Fit Hawkes Process via MLE
#' @description Estimates parameters \eqn{\lambda, \alpha, \beta} using Maximum 
#'   Likelihood.
#' @param object Either a numeric vector of event times or a 'point_process_sim'
#'   object.
#' @param T_max Numeric value. Required if 'object' is a numeric vector.
#' @param lik_method Character string, either "fast" (recursive) or "slow" 
#'   (direct).
#' @param init Numeric vector. Initial guesses for the parameters.
#' @return An object of class 'hawkes_fit'.
#' @export
exp_hp_fit <- function(object, 
                       T_max = NULL, 
                       lik_method = c("fast", "slow"), 
                       init = c(0.1, 0.5, 1.0)) {
  
  # defaults to "fast"
  lik_method <- match.arg(lik_method)
  
  # handle dual-input ("point_process_sim" class vs. numeric vector)
  if (inherits(object, "point_process_sim")) {
    H_t   <- object$events
    T_max <- object$T_max
  } else {
    H_t <- object
    if (is.null(T_max)) {
      stop("T_max must be provided when 'object' is a numeric vector.")
    }
  }
  
  # optimisation, passing 'H_t', 'T_max', and 'method as extra arguments
  fit_out <- stats::optim(
    par = init,
    fn = .exp_hp_negloglik, 
    H_t = H_t,
    T_max = T_max,
    lik_method = lik_method,
    method  = "L-BFGS-B",
    lower   = c(1e-6, 1e-6, 1e-6), # Ensure parameters stay positive
    hessian = TRUE
  )
  
  # for the process to be stable, alpha / beta must be < 1
  is_stable <- fit_out$par[2] < fit_out$par[3]
  if (!is_stable) {
    warning("Estimated process is non-stationary (alpha >= beta). Results may be unreliable.")
  }
  
  # return the 'hawkes_fit' Object
  structure(
    list(
      par         = stats::setNames(fit_out$par, c("lambda", "alpha", "beta")),
      hessian     = fit_out$hessian,
      loglik      = -fit_out$value,
      H_t         = H_t,
      T_max       = T_max,
      convergence = fit_out$convergence,
      n           = length(H_t),
      stable      = is_stable
    ),
    class = "hawkes_fit"
  )
}

#' @export
logLik.hawkes_fit <- function(object, ...) {
  structure(
    object$loglik,
    df = length(object$par),
    nobs = object$n,
    class = "logLik"
  )
}

#' @export
summary.hawkes_fit <- function(object, ...) {
  v_cov <- solve(object$hessian)
  se <- sqrt(diag(v_cov))
  estimates <- object$par
  z_val <- estimates / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
  
  coef_mat <- cbind(
    Estimate = estimates,
    `Std. Error` = se,
    `z value` = round(z_val, 3),
    `Pr(>|z|)` = p_val
  )
  
  branching_ratio <- estimates[2] / estimates[3]
  aic <- 2*length(estimates) - 2*object$loglik
  dev <- -2*object$loglik
  
  structure(
    list(
      coefficients = coef_mat,
      aic = aic,
      deviance = dev,
      n = object$n,
      branching_ratio = branching_ratio,
      stable = object$stable
    ),
    class = "summary.hawkes_fit"
  )
}

#' @export
print.summary.hawkes_fit <- function(x, ...) {
  cat("\n--- Hawkes Process MLE Summary ---\n")
  cat("Observations:          ", x$n, "\n")
  cat("AIC:                   ", round(x$aic, 2), "\n")
  cat("-2 Log-Lik:            ", round(x$deviance, 2), "\n")
  
  cat("Branching Ratio (est): ", round(x$branching_ratio, 4))
  if (!x$stable) {
    cat(" (Non-Stationary!)\n")
  } else {
    cat(" (Stationary)\n")
  }
  
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, 
                      P.values = TRUE, 
                      has.Pvalue = TRUE, 
                      digits = 4)
  cat("---\n")
}

#' @title Plot a Fitted Hawkes Process
#' @description S3 method for the 'hawkes_fit' class. Plots the estimated 
#'   conditional intensity function against the observed event times.
#' @param x An object of class 'hawkes_fit'.
#' @param n_grid Integer. The number of points used to calculate the 
#'   intensity curve. Default is 1000.
#' @param ... Additional arguments passed to the plot.
#' @export
plot.hawkes_fit <- function(x, n_grid = 1000, ...) {
  
  # when generating time grid we add the actual event times to the grid to show
  # the immdeiate jump in intensity
  t_grid <- sort(unique(c(
    seq(0, x$T_max, length.out = n_grid),
    x$H_t,
    x$H_t + 1e-9  # Tiny offset to show the vertical increase after a jump
  )))
  
  # calculate intensity on the grid (using our internal recursive helper)
  intensities <- .get_exp_intensity_grid(t_grid, x$H_t, x$par)
  
  plot_base_canvas(
    xlim  = c(0, x$T_max), 
    ylim  = c(0, max(intensities)),
    title = "Fitted Hawkes Intensity"
  )
  graphics::abline(h = x$par[1], col = "slategrey", lty = 3, lwd = 1.2)
  add_intensity(t_grid, intensities,)
  add_events(x$H_t)
  m_text <- sprintf("lambda: %.3f, alpha: %.3f, beta: %.3f", 
                    x$par[1], x$par[2], x$par[3])
  graphics::mtext(m_text, side = 3, line = 0.2, cex = 0.8)
  
  invisible(x)
}

#' @title Calculate Transformed Residuals
#' @description Calculates the Time-scale transformed residuals for a Hawkes fit
#'   using Thm's 9.1 and 9.2 from Laub, Taimre, Pollett (2021) by converting the
#'   Hawkes process into a unit rate Poisson process by taking the value of the 
#'   compensator at each point as the new transformed process.
#' @param object An object of class 'hawkes_fit'.
#' @param ... Additional arguments.
#' @importFrom stats residuals
#' @export
residuals.hawkes_fit <- function(object, ...) {
  H_t <- object$H_t
  theta <- object$par
  
  #### OLD SLOW (O(n^2)) LOGIC ####
  # For each event i, calculate the compensator value at time t_i
  # tau <- sapply(seq_along(H_t), function(i) {
  #   ti <- H_t[i]
  #   hi <- H_t[1:(i - 1)] # Only use history up to this point
  #   exp_hp_compensator(theta = theta, H_t = hi, t = ti)
  # })
  # tau
  
  #### NEW FAST LOGIC ####
  # on a test hp realisation with N_T = 40,000, went from 17s (old) to 0.5s (new)!
  .residuals_fast(theta = object$par, H_t = object$H_t)
  
}
#' @keywords internal
.residuals_fast <- function(theta, H_t) {
  n <- length(H_t)
  if (n == 0) return(numeric(0))
  
  # unpack parameters
  lambda <- theta[1]
  alpha <- theta[2]
  beta <- theta[3]
  
  # get intensities with O(n) method
  intensity <- .exp_hp_intensity_at_events(theta, H_t)
  
  # calculate Compensator at time t_i using markovian property 
  B <- numeric(n)
  B[1] <- lambda*H_t[1]
  if (n > 1) {
    for (i in 2:n) {
      dt <- H_t[i] - H_t[i-1]
      B[i] <- B[i-1] + lambda*dt + (intensity[i-1] - lambda + alpha) * 
        (1/beta) * (1 - exp(-beta*dt))
    }
  }
  return(B)
}


#' @title Q-Q Plot for Hawkes Residuals
#' @description Checks goodness-of-fit by comparing transformed inter-arrival 
#'   times to a theoretical Exponential(1) distribution.
#' @param object An object of class 'hawkes_fit'.
#' @export
plot_residuals <- function(object) {
  # uses Thm 9.1 (Random Time Change Theorem) from Laub, Taimre, Pollett (2021)
  tau <- residuals(object)

  # convert to inter-arrival times ()
  gaps <- diff(c(0, tau))
  n <- length(gaps)
  
  max_val <- max(stats::qexp(stats::ppoints(n)), gaps) * 1.05
  
  plot_base_canvas(
    xlim  = c(0, max_val), 
    ylim  = c(0, max_val), 
    title = "Residual Q-Q Plot",
    xlab  = "Theoretical Exp(1) Quantiles",
    ylab  = "Sample Residual Quantiles"
  )
  graphics::abline(0, 1, col = "slategrey", lty = 2, lwd = 2)
  theoretical_q <- stats::qexp(stats::ppoints(n))
  graphics::points(theoretical_q, sort(gaps), pch = 20, cex = 0.5)
  
  # subtitle for the Kolmogorov-Smirnov test
  ks_res <- stats::ks.test(gaps, "pexp", rate = 1)
  ks_text <- sprintf("K-S Test p-value: %.4f", ks_res$p.value)
  graphics::mtext(ks_text, side = 3, line = 0.2, cex = 0.8)
  
  invisible(object)
}
