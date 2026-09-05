################################################################################
# Univariate exponential Hawkes model
################################################################################

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

  lambda + alpha * sum(exp(-beta * (t - H_t)))
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
      A[i] <- exp(-beta * (H_t[i] - H_t[i - 1])) * (1 + A[i - 1])
    }
  }
  lambda + alpha * A
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

  if (t <= 0) {
    return(0)
  }
  lambda * t + (alpha / beta) * sum(1 - exp(-beta * (t - H_t)))
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
  if (any(theta <= 0)) {
    return(-Inf)
  }

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
.exp_hp_negloglik <- function(
  theta,
  H_t,
  T_max,
  lik_method = c("fast", "slow")
) {
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
exp_hp_fit <- function(
  object,
  T_max = NULL,
  lik_method = c("fast", "slow"),
  init = c(0.1, 0.5, 1.0)
) {
  # defaults to "fast"
  lik_method <- match.arg(lik_method)

  # handle dual-input ("point_process_sim" class vs. numeric vector)
  if (inherits(object, "point_process_sim")) {
    H_t <- object$events
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
    method = "L-BFGS-B",
    lower = c(1e-6, 1e-6, 1e-6), # Ensure parameters stay positive
    hessian = TRUE
  )

  # for the process to be stable, alpha / beta must be < 1
  is_stable <- fit_out$par[2] < fit_out$par[3]
  if (!is_stable) {
    warning(
      "Estimated process is non-stationary (alpha >= beta). Results may be unreliable."
    )
  }

  # return the 'hawkes_fit' Object
  structure(
    list(
      par = stats::setNames(fit_out$par, c("lambda", "alpha", "beta")),
      hessian = fit_out$hessian,
      loglik = -fit_out$value,
      H_t = H_t,
      T_max = T_max,
      convergence = fit_out$convergence,
      n = length(H_t),
      stable = is_stable
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
  aic <- 2 * length(estimates) - 2 * object$loglik
  dev <- -2 * object$loglik

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
  stats::printCoefmat(
    x$coefficients,
    P.values = TRUE,
    has.Pvalue = TRUE,
    digits = 4
  )
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
    x$H_t + 1e-9 # Tiny offset to show the vertical increase after a jump
  )))

  # calculate intensity on the grid (using our internal recursive helper)
  intensities <- .get_exp_intensity_grid(t_grid, x$H_t, x$par)

  plot_base_canvas(
    xlim = c(0, x$T_max),
    ylim = c(0, max(intensities)),
    title = "Fitted Hawkes Intensity"
  )
  graphics::abline(h = x$par[1], col = "slategrey", lty = 3, lwd = 1.2)
  add_intensity(t_grid, intensities, )
  add_events(x$H_t)
  m_text <- sprintf(
    "lambda: %.3f, alpha: %.3f, beta: %.3f",
    x$par[1],
    x$par[2],
    x$par[3]
  )
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
  if (n == 0) {
    return(numeric(0))
  }

  # unpack parameters
  lambda <- theta[1]
  alpha <- theta[2]
  beta <- theta[3]

  # get intensities with O(n) method
  intensity <- .exp_hp_intensity_at_events(theta, H_t)

  # calculate Compensator at time t_i using markovian property
  B <- numeric(n)
  B[1] <- lambda * H_t[1]
  if (n > 1) {
    for (i in 2:n) {
      dt <- H_t[i] - H_t[i - 1]
      B[i] <- B[i - 1] +
        lambda * dt +
        (intensity[i - 1] - lambda + alpha) *
          (1 / beta) *
          (1 - exp(-beta * dt))
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
    xlim = c(0, max_val),
    ylim = c(0, max_val),
    title = "Residual Q-Q Plot",
    xlab = "Theoretical Exp(1) Quantiles",
    ylab = "Sample Residual Quantiles"
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


################################################################################
# Marked exponential Hawkes model
#  - linear mark effect: A*m
################################################################################

#' @title Exponential Marked Hawkes Process Intensity
#' @description Calculates the conditional ground-process intensity at time
#'   \eqn{t} for a marked Hawkes process with excitation
#'   \eqn{A m \exp(-\beta (t-s))}.
#' @param t Numeric value. The time at which to evaluate the intensity.
#' @param H_t Numeric vector of event timestamps occurring before \eqn{t}.
#' @param marks Numeric vector of marks corresponding to `H_t`.
#' @param theta Numeric vector of parameters \eqn{c(\lambda_0, A, \beta)}.
#' @export
exp_mhp_intensity <- function(t, H_t, marks, theta) {
  lambda_0 <- theta[1]
  A <- theta[2]
  beta <- theta[3]

  if (length(H_t) == 0) {
    return(lambda_0)
  }

  lambda_0 + A * sum(marks * exp(-beta * (t - H_t)))
}

#' @keywords internal
.exp_mhp_intensity_at_events <- function(theta, H_t, marks) {
  lambda_0 <- theta[1]
  A <- theta[2]
  beta <- theta[3]

  n <- length(H_t)
  R <- numeric(n)

  if (n > 1) {
    for (i in 2:n) {
      dt <- H_t[i] - H_t[i - 1]
      R[i] <- exp(-beta * dt) * (marks[i - 1] + R[i - 1])
    }
  }

  lambda_0 + A * R
}

#' @title Exponential Marked Hawkes Process Compensator
#' @description Calculates the integrated ground-process intensity over
#'   \eqn{[0,t]}.
#' @param theta Numeric vector of parameters \eqn{c(\lambda_0, A, \beta)}.
#' @param H_t Numeric vector of event timestamps.
#' @param marks Numeric vector of marks corresponding to `H_t`.
#' @param t Numeric value. The end of the observation window.
#' @export
exp_mhp_compensator <- function(theta, H_t, marks, t) {
  lambda_0 <- theta[1]
  A <- theta[2]
  beta <- theta[3]

  if (t <= 0) {
    return(0)
  }

  lambda_0 * t + (A / beta) * sum(marks * (1 - exp(-beta * (t - H_t))))
}

#' @title Marked Hawkes Log-Likelihood (Exponential)
#' @description Calculates the ground-process log-likelihood for a marked
#'   Hawkes process with excitation \eqn{A m \exp(-\beta t)}.
#' @param theta Numeric vector of parameters \eqn{c(\lambda_0, A, \beta)}.
#' @param H_t Numeric vector of event timestamps.
#' @param marks Numeric vector of marks corresponding to `H_t`.
#' @param T_max Numeric value. The end of the observation window.
#' @param lik_method Character string, either `"fast"` or `"slow"`.
#' @export
exp_mhp_loglik <- function(
  theta,
  H_t,
  marks,
  T_max,
  lik_method = c("fast", "slow")
) {
  lik_method <- match.arg(lik_method)

  if (length(H_t) != length(marks)) {
    stop("H_t and marks must have the same length.")
  }

  if (any(theta <= 0)) {
    return(-Inf)
  }

  if (lik_method == "fast") {
    .exp_mhp_loglik_fast(theta, H_t, marks, T_max)
  } else {
    .exp_mhp_loglik_slow(theta, H_t, marks, T_max)
  }
}

#' @keywords internal
.exp_mhp_loglik_fast <- function(theta, H_t, marks, T_max) {
  intensities <- .exp_mhp_intensity_at_events(
    theta = theta,
    H_t = H_t,
    marks = marks
  )

  if (any(intensities <= 0) || any(!is.finite(intensities))) {
    return(-Inf)
  }

  sum(log(intensities)) -
    exp_mhp_compensator(
      theta = theta,
      H_t = H_t,
      marks = marks,
      t = T_max
    )
}

#' @keywords internal
.exp_mhp_loglik_slow <- function(theta, H_t, marks, T_max) {
  n <- length(H_t)

  if (n == 0) {
    return(-exp_mhp_compensator(theta, H_t, marks, T_max))
  }

  log_intensities <- vapply(
    seq_len(n),
    function(i) {
      if (i == 1) {
        lambda_0 <- theta[1]
        return(log(lambda_0))
      }

      past_events <- H_t[seq_len(i - 1)]
      past_marks <- marks[seq_len(i - 1)]

      intensity <- exp_mhp_intensity(
        t = H_t[i],
        H_t = past_events,
        marks = past_marks,
        theta = theta
      )

      if (intensity <= 0 || !is.finite(intensity)) {
        return(-Inf)
      }

      log(intensity)
    },
    numeric(1)
  )

  if (any(!is.finite(log_intensities))) {
    return(-Inf)
  }

  sum(log_intensities) -
    exp_mhp_compensator(theta, H_t, marks, T_max)
}

#' @keywords internal
.exp_mhp_negloglik <- function(
  theta,
  H_t,
  marks,
  T_max,
  lik_method = c("fast", "slow")
) {
  -exp_mhp_loglik(
    theta = theta,
    H_t = H_t,
    marks = marks,
    T_max = T_max,
    lik_method = match.arg(lik_method)
  )
}

#' @title Fit Marked Hawkes Process via MLE
#' @description Estimates \eqn{\lambda_0}, \eqn{A}, and \eqn{\beta} by maximum
#'   likelihood for a marked Hawkes process with exponential temporal decay
#'   and linear mark productivity.
#' @param object Either a numeric vector of event times or a
#'   `marked_pp_sim` object.
#' @param marks Numeric vector of event marks when `object` is a numeric vector.
#'   Not required when `object` is a `marked_pp_sim`.
#' @param T_max Numeric value. Required if `object` is a numeric vector.
#' @param lik_method Character string, either `"fast"` or `"slow"`.
#' @param init Numeric vector of initial values for
#'   \eqn{c(\lambda_0,A,\beta)}.
#' @return An object of class `marked_hawkes_fit`.
#' @export
exp_mhp_fit <- function(
  object,
  marks = NULL,
  T_max = NULL,
  lik_method = c("fast", "slow"),
  init = c(0.5, 1, 1)
) {
  lik_method <- match.arg(lik_method)

  if (inherits(object, "marked_pp_sim")) {
    H_t <- object$events
    marks <- object$marks
    T_max <- object$T_max
  } else {
    H_t <- object

    if (is.null(marks)) {
      stop("marks must be provided when object is a numeric vector.")
    }

    if (is.null(T_max)) {
      stop("T_max must be provided when object is a numeric vector.")
    }
  }

  if (length(H_t) != length(marks)) {
    stop("H_t and marks must have the same length.")
  }

  if (length(init) != 3L || any(!is.finite(init)) || any(init <= 0)) {
    stop("init must contain three positive finite values.")
  }

  fit_out <- stats::optim(
    par = init,
    fn = .exp_mhp_negloglik,
    H_t = H_t,
    marks = marks,
    T_max = T_max,
    lik_method = lik_method,
    method = "L-BFGS-B",
    lower = c(1e-8, 1e-8, 1e-8),
    hessian = TRUE
  )

  is_stable <- fit_out$par[2] * mean(marks) < fit_out$par[3]

  if (!is_stable) {
    warning(
      "Estimated process is non-stationary under the empirical mean mark."
    )
  }

  structure(
    list(
      par = stats::setNames(
        fit_out$par,
        c("lambda_0", "A", "beta")
      ),
      hessian = fit_out$hessian,
      loglik = -fit_out$value,
      H_t = H_t,
      marks = marks,
      T_max = T_max,
      convergence = fit_out$convergence,
      n = length(H_t),
      stable = is_stable
    ),
    class = "marked_hawkes_fit"
  )
}


#' @export
logLik.marked_hawkes_fit <- function(object, ...) {
  structure(
    object$loglik,
    df = length(object$par),
    nobs = object$n,
    class = "logLik"
  )
}


#' @export
summary.marked_hawkes_fit <- function(object, ...) {
  v_cov <- tryCatch(
    solve(object$hessian),
    error = function(e) NULL
  )

  if (is.null(v_cov) || any(!is.finite(v_cov))) {
    warning("Unable to invert Hessian; standard errors unavailable.")
    se <- rep(NA_real_, length(object$par))
  } else {
    se <- sqrt(pmax(diag(v_cov), 0))
  }

  estimates <- object$par
  z_val <- estimates / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  coef_mat <- cbind(
    Estimate = estimates,
    `Std. Error` = se,
    `z value` = round(z_val, 3),
    `Pr(>|z|)` = p_val
  )

  branching_ratio <- estimates["A"] * mean(object$marks) / estimates["beta"]
  aic <- 2 * length(estimates) - 2 * object$loglik
  dev <- -2 * object$loglik

  structure(
    list(
      coefficients = coef_mat,
      aic = aic,
      deviance = dev,
      n = object$n,
      branching_ratio = unname(branching_ratio),
      stable = object$stable
    ),
    class = "summary.marked_hawkes_fit"
  )
}


#' @export
print.summary.marked_hawkes_fit <- function(x, ...) {
  cat("\n--- Marked Hawkes Process MLE Summary ---\n")
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
  stats::printCoefmat(
    x$coefficients,
    P.values = TRUE,
    has.Pvalue = TRUE,
    digits = 4
  )
  cat("---\n")
}

#' @title Plot a Fitted Marked Hawkes Process
#' @description S3 method for the `marked_hawkes_fit` class. Plots the fitted
#'   ground-process conditional intensity against the observed event times.
#' @param x An object of class `marked_hawkes_fit`.
#' @param show_marks Boolean. If `TRUE` displays mark magnitude plot below
#'   fitted intensity
#' @param n_grid Integer. The number of points used to calculate the
#'   intensity curve. Default is 1000.
#' @param ... Additional arguments passed to the plot.
#' @export
plot.marked_hawkes_fit <- function(x, show_marks = TRUE, n_grid = 1000, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (show_marks) {
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  }

  t_grid <- sort(unique(c(
    seq(0, x$T_max, length.out = n_grid),
    x$H_t,
    x$H_t + 1e-9
  )))

  intensities <- .get_exp_mhp_intensity_grid(
    t_grid = t_grid,
    H_t = x$H_t,
    marks = x$marks,
    theta = x$par
  )

  plot_base_canvas(
    xlim = c(0, x$T_max),
    ylim = c(0, max(intensities)),
    title = "Fitted Marked Hawkes Intensity"
  )

  graphics::abline(h = x$par[1], col = "slategrey", lty = 3, lwd = 1.2)

  add_intensity(t_grid, intensities)
  add_events(x$H_t)

  m_text <- sprintf(
    "lambda_0: %.3f, A: %.3f, beta: %.3f",
    x$par[1],
    x$par[2],
    x$par[3]
  )

  if (show_marks) {
    plot(
      x$H_t,
      x$marks,
      type = "h",
      xlim = c(0, x$T_max),
      ylim = c(0, max(x$marks, na.rm = TRUE)),
      xlab = "Time (t)",
      ylab = "Mark (m)",
      main = "Mark Magnitudes",
    )

    points(
      x$H_t,
      x$marks,
      pch = 21,
      bg = "slategrey",
      col = "black",
      cex = 0.8
    )
  } else {
    graphics::mtext(m_text, side = 3, line = 0.2, cex = 0.8)
  }

  invisible(x)
}

#' @keywords internal
.residuals_marked_fast <- function(theta, H_t, marks) {
  lambda_0 <- theta[1]
  A <- theta[2]
  beta <- theta[3]

  n <- length(H_t)

  if (n == 0) {
    return(numeric(0))
  }

  tau <- numeric(n)

  # First event has no previous excitation
  tau[1] <- lambda_0 * H_t[1]

  R <- 0

  if (n > 1) {
    for (i in 2:n) {
      dt <- H_t[i] - H_t[i - 1]

      # Excitation immediately after T_{i-1}, excluding A
      R_after <- marks[i - 1] + R

      # Integrate intensity over (T_{i-1}, T_i]
      tau[i] <- tau[i - 1] +
        lambda_0 * dt +
        A * R_after * (1 / beta) * (1 - exp(-beta * dt))

      # Update state to just before T_i
      R <- exp(-beta * dt) * R_after
    }
  }

  tau
}

#' @title Calculate Transformed Residuals for a Marked Hawkes Fit
#' @description Calculates time-rescaled residuals using the compensator of the
#'   fitted ground process. For unpredictable marks, the mark density does not
#'   enter the transformed event-time residuals.
#' @param object An object of class `marked_hawkes_fit`.
#' @param ... Additional arguments.
#' @importFrom stats residuals
#' @export
residuals.marked_hawkes_fit <- function(object, ...) {
  .residuals_marked_fast(
    theta = object$par,
    H_t = object$H_t,
    marks = object$marks
  )
}
