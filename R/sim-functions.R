#' Internal Constructor for point_process_sim
#' @keywords internal
.new_point_process_sim <- function(events, T_max, params, type, 
                                   marks = NULL, full_history = NULL) {
  out <- list(
    events = sort(events),
    marks = marks, # New field
    T_max = T_max,
    params = params,
    n = length(events),
    type = type,
    full_history = full_history
  )
  
  # Add "marked_pp_sim" as a subclass if marks exist
  class_vector <- "point_process_sim"
  if (!is.null(marks)) class_vector <- c("marked_pp_sim", class_vector)
  
  structure(out, class = class_vector)
}


#' Simulation of a Poisson process
#' 
#' Simulation of a (possibly inhomogeneous) Poisson process using the thinning 
#' algorithm in  Laub, Taimre, and Pollett (2021) Chapter 4 (Algorithm 1). It 
#' works by generating a homogeneous Poisson process with rate `M` on 
#' \eqn{[0,T_{max}]}, and then, in the inhomogeneous case, thinning this process 
#' using the provided intensity function `intens`
#'
#' @param T_max A non-negative numeric value - the end of the interval 
#'   \eqn{[0,T_{max}]}.
#' @param intens A non-negative function (or numeric if simulating a homogeneous
#'   Poisson process) representing the intensity function, \eqn{\lambda (t)} of 
#'   the Poisson process.
#' @param M A non-negative numeric value - upper bound on `fun`. Ignored for 
#'   a homogeneous Poisson process.
#'
#' @return An object of class `point_process_sim`. A list containing:
#' \itemize{
#'   \item \code{events}: A sorted numeric vector of arrival times.
#'   \item \code{T_max}: The end of the simulation window.
#'   \item \code{params}: A list of the parameters used for simulation.
#'   \item \code{n}: The total number of events.
#'   \item \code{type}: The type of process simulated.
#' }
#' @export
#'
#' @examples
#' # simulate an inhomogeneous Poisson process with Exp(2) intensity
#' t <- sim_pp(T_max = 5, intens = function(t) {2*exp(-2*t)}, M = 2)
#' # Simulate a homogeneous Poisson process with λ = 2
#' x <- sim_pp(T_max = 5, intens = 2)
sim_pp <- function(T_max, intens, M = NULL) {
  
  # Homogeneous Case 
  if (is.numeric(intens)) {
    stopifnot(
      "'T_max' must be positive" = is.numeric(T_max) && length(T_max) == 1 && T_max > 0,
      "'intens' must be a single non-negative number" = 
        length(intens) == 1 && !is.na(intens) && intens >= 0
    )
    N_T <- stats::rpois(1, intens * T_max)
    arrival_times <- sort(stats::runif(N_T, 0, T_max))
  } else {
    # Inhomogeneous Case Guardrails
    # TODO: errors for missing 'M', could instead approximate across a grid
    stopifnot(
      "'T_max' must be positive" = is.numeric(T_max) && length(T_max) == 1 && T_max > 0,
      "'intens' must be a function" = is.function(intens),
      "Inhomogeneous simulation requires a single numeric 'M' >= 0" = 
        is.numeric(M) && length(M) == 1 && !is.na(M) && M >= 0
    )
    
    # simulate homogeneous poisson process with rate M on [0, T_max]
    N_T <- stats::rpois(1, M * T_max)
    if (N_T == 0) {
      arrival_times <- numeric(0) # Handle edge case of 0 arrivals
    } else {
      unthinned_arrival_times <- sort(stats::runif(N_T, 0, T_max))
      # keep arrivals with some probability
      # intens() must be able to handle the vector 'unthinned_arrival_times'
      u <- stats::runif(N_T, 0, M)
      keep <- u <= intens(unthinned_arrival_times)
      arrival_times <- unthinned_arrival_times[keep]
    }
  }
  
  # Wrap the results in a point_process_sim class
  .new_point_process_sim(
    events = arrival_times,
    T_max = T_max,
    params = list(intens = intens),
    type = "Poisson"
  )
}

#' Simulation of a Hawkes process by thinning
#' 
#' Simulation of a Hawkes process using Ogata's modified thinning algorithm, 
#' following Chapter 4 (Algorithm 2) of Laub, Taimre, and Pollett (2021). 
#' 
#' @details 
#' This method is similar to the thinning method for an inhomogeneous Poisson 
#' process. However, since the Hawkes intensity depends on the event history, a 
#' global upper bound for the conditional intensity is typically not available. 
#' Instead, we assume that the excitation function \eqn{\mu(\cdot)} is 
#' non-increasing. This ensures that, between arrivals, the conditional 
#' intensity is non-increasing in time. As a result, the current intensity can 
#' be used as a valid local upper bound for the next candidate arrival, which is
#' updated after each event.
#'
#' @param T_max A non-negative numeric value - the end of the interval 
#'   \eqn{[0,T_{max}]}.
#' @param lambda A non-negative numeric value - the background arrival 
#'   component \eqn{\lambda} of the Hawkes conditional intensity
#' @param mu A non-negative, non-increasing function - the excitation 
#'   function/kernel \eqn{\mu(\cdot)} of the Hawkes conditional intensity.
#'
#' @return An object of class `point_process_sim`. A list containing:
#' \itemize{
#'   \item \code{events}: A sorted numeric vector of arrival times.
#'   \item \code{T_max}: The end of the simulation window.
#'   \item \code{params}: A list of the parameters used for simulation.
#'   \item \code{n}: The total number of events.
#' }
#' @export
#'
#' @examples
#' # simulate a Hawkes process with background rate 1 and Exp(2) kernel
#' t <- sim_hp(T_max=5, lambda = 1, mu = function(t) {2*exp(-2*t)})
sim_hp <- function(T_max, lambda, mu) {
  
  # some guardrails
  stopifnot(
    "'T_max' must be a single numeric" = 
      is.numeric(T_max) && length(T_max) == 1 && !is.na(T_max),
    
    "'T_max' must be positive" 
    = T_max > 0,
    
    "'lambda' must be a single numeric" 
    = is.numeric(lambda) && length(lambda) == 1 && !is.na(lambda),
    
    "'lambda' must be non-negative" 
    = lambda >= 0,
    
    "'mu' must be a function" 
    = is.function(mu)
  )
  
  # Start with a sensible ~8KB buffer for vector (will double later if required)
  arrival_times <- numeric(1024)
  t <- 0
  count <- 0
  
  # initialise simulation objects
  t <- 0
  while (t < T_max) {
    
    # find new upper bound, using only the filled portion of arrival_times
    M <- lambda + sum(mu(t - arrival_times[seq_len(count)]))
    if (is.na(M)) {
      stop("Intensity calculation returned NA. Check your 'mu' function.")
    }
    
    # generate next candidate point
    t <- t + stats::rexp(1, M)
    if (t > T_max) break
    
    current_intensity <- lambda + sum(mu(t - arrival_times[seq_len(count)]))
    # keep it with some probability
    if (stats::runif(1) <= current_intensity/M) {
      count <- count + 1
      
      # double the capacity of arrival_times if out of space
      if (count > length(arrival_times)) {
        arrival_times <- c(arrival_times, numeric(length(arrival_times)))
      }
      arrival_times[count] <- t
    }
  }
  
  # trim trailing zeros
  arrival_times <- if (count == 0) numeric(0) else arrival_times[1:count]
  
  # Wrap the results in a point_process_sim class
  .new_point_process_sim(
    events = arrival_times,
    T_max = T_max,
    params = list(lambda = lambda, mu = mu),
    type = "Hawkes"
  )
}

#' Simulation of a Hawkes process by clusters
#' 
#' Simulation of a Hawkes process by clusters i.e. using the immigrant-birth 
#' representation of a Hawkes process. It follows Chapter 4 (Algorithm 3) of
#' Laub, Taimre, and Pollett (2021). For simplicity, we factor the usual 
#' Hawkes kernel as \eqn{\mu(t) = \eta g(t)} where 
#' \eqn{\eta = \int_0^{\infty} \mu(t) dt} is the branching ratio, and \eqn{g(t)}
#' is a probability density function on \eqn{\mathbb{R}^+}. 
#' 
#' Simulation via an immigrant-birth process proceeds in two stages. The first 
#' stage consists of randomly generating immigrant arrivals according to a 
#' Poisson process with intensity \eqn{\lambda} on the interval 
#' \eqn{[0, T_{max}]}.
#' The second stage then proceeds by randomly generating 
#' \eqn{\mathrm{Poisson}(\eta)} first-generation offspring for each immigrant 
#' generated, where conditional on an immigrant arriving at time \eqn{s}, the 
#' waiting times of its offspring are independent draws from the density 
#' \eqn{g(t)}. Thus, first-generation offspring arrival times are given by 
#' \deqn{s + E_j, \quad E_j \sim g(·)}
#' Each first-generation offspring in turn generates further descendants 
#' according to the same procedure, which produces the characteristic branching 
#' structure. The recursion continues until no further offspring are generated 
#' within \eqn{[0, T_{max}]}.
#'
#' @param T_max A non-negative numeric value - the end of the interval 
#'   \eqn{[0,T_{max}]}.
#' @param lambda A non-negative numeric value - the background arrival 
#'   component \eqn{\lambda} of the Hawkes conditional intensity
#' @param eta A non-negative numeric value - the branching ratio of the Hawkes
#'   process
#' @param family A character string "family" which corresponds to a distribution
#'   with a random generation function "rfamily". 
#' @param ... Additional arguments passed to the random generation function.
#'
#' @return An object of class `point_process_sim`. A list containing:
#' \itemize{
#'   \item \code{events}: A sorted numeric vector of arrival times.
#'   \item \code{T_max}: The end of the simulation window.
#'   \item \code{params}: A list of the parameters used for simulation.
#'   \item \code{n}: The total number of events.
#'   \item \code{history}: A dataframe containing the simulation with
#'     additional fields to rebuild the history of immigrant-offspring branches.
#' }
#' @export
#' 
sim_hp_clus <- function(T_max, lambda, eta, family, ...) {
  
  # some guardrails
  stopifnot(
    "'T_max' must be a single numeric" = 
      is.numeric(T_max) && length(T_max) == 1 && !is.na(T_max),
    
    "'T_max' must be positive" = 
      T_max > 0,
    
    "'lambda' must be a single numeric" = 
      is.numeric(lambda) && length(lambda) == 1 && !is.na(lambda),
    
    "'lambda' must be non-negative" = 
      lambda >= 0,
    
    "'eta' (branching ratio) must be a single numeric" = 
      is.numeric(eta) && length(eta) == 1 && !is.na(eta),
    
    "'eta' must be non-negative" = 
      eta >= 0,
    
    "'family' must be a single character string" = 
      is.character(family) && length(family) == 1
  )
  
  # warning for possible explosion
  if (eta >= 1) {
    warning("eta >= 1: The process is non-stationary and may not terminate.")
  }
  
  # Check for the existence of the distribution function
  dist_name <- paste0("r", family)
  if (!exists(dist_name, mode = "function")) {
    stop(sprintf("The distribution function '%s' does not exist.", dist_name))
  } else {
    random_fn <- get(paste0("r", family), mode = "function")
  }
  
  # Stage 1: Generate immigrants using sim_pp function
  # TODO: possibility for inghomogeneous background (not standard Hawkes formulation though)
  immigrant_obj <- sim_pp(T_max, lambda)
  immigrant_times <- immigrant_obj$events
  n_immigrants <- length(immigrant_times)
  
  # Store results: time, generation, id, parent_id
  res <- data.frame(
    time = immigrant_times,
    gen = 0,
    id = seq_len(n_immigrants),
    parent_id = 0
  )
  
  # 'recent_gen_df are the arrivals that can still "give birth"
  recent_gen_df <- res
  next_id <- n_immigrants + 1
  
  # Stage 2: Recursive branching (Generation by Generation)
  while (nrow(recent_gen_df) > 0) {
    # number of offspring for each event in the recent generation
    num_offspring <- stats::rpois(nrow(recent_gen_df), eta)
    
    # if no offspring in this entire generation, we are done
    if (sum(num_offspring) == 0) break
    
    # map children to their specific parents
    parent_indices <- rep(seq_len(nrow(recent_gen_df)), times = num_offspring)
    parents <- recent_gen_df[parent_indices, ]
    
    # generate child arrival times
    child_times <- parents$time + random_fn(nrow(parents), ...)
    
    # filter by T_max and build the next generation dataframe slice
    in_window <- child_times <= T_max
    if (!any(in_window)) break
    
    next_gen_df <- data.frame(
      time = child_times[in_window],
      gen = parents$gen[in_window] + 1,
      id = seq(next_id, length.out = sum(in_window)),
      parent_id = parents$id[in_window]
    )
    
    # update for next loop and storage
    res <- rbind(res, next_gen_df)
    recent_gen_df <- next_gen_df
    next_id <- next_id + sum(in_window)
  }
  
  # Handle case where no immigrants are generated
  if (n_immigrants == 0) {
    arrival_times <- numeric(0)
    full_history  <- data.frame()
  } else {
    # Sort the full dataframe by time
    res <- res[order(res$time), ]
    arrival_times <- res$time
    full_history  <- res
  }
  
  # Warp the results in a point_process_sim class
  .new_point_process_sim(
    events = arrival_times,
    T_max = T_max,
    params = list(lambda = lambda, eta = eta, family = family, dots = list(...)),
    type = "Hawkes (cluster)",
    full_history = res # This is the dataframe with gen/parent info
  )
}

#' Plot a Point Process Simulation
#' @description S3 method for the point_process_sim class.
#' @param x An object of class 'point_process_sim'.
#' @param ... Additional arguments (ignored for now).
#' @export
plot.point_process_sim <- function(x, ...) {
  # Setup Canvas
  # Height goes up to total number of events (n)
  plot_base_canvas(
    xlim  = c(0, x$T_max), 
    ylim  = c(0, x$n), 
    title = "Point Process Simulation: N(t)",
    ylab  = expression(N(t))
  )
  add_counting_process(x$events, x$T_max)
  add_events(x$events)
  invisible(x)
}

#' @export
print.point_process_sim <- function(x, ...) {
  cat("\n--- Point Process Simulation ---\n")
  cat("Type:   ", x$type, "\n")
  cat("Events: ", x$n, "\n")
  cat("Window: [0, ", x$T_max, "]\n", sep = "")
  if (x$n > 0) {
    cat(
      "Times:  ", 
      paste(round(utils::head(x$events, 3), 3), collapse = ", "), 
      "...\n"
    )
  }
  cat("---\n")
}

#' Simulation of a Marked Hawkes Process by thinning
#' 
#' Simulation of a marked Hawkes process using a generalized version of Ogata's 
#' thinning algorithm. This function allows for a flexible conditional intensity 
#' that can depend on the history of both event times and their associated marks.
#'
#' @details 
#' The simulation assumes that the `intensity_func` provided is non-increasing 
#' between events. This allows the intensity at the time of the most recent 
#' event (or the background rate) to serve as a local upper bound for the 
#' next candidate arrival. 
#' 
#' For each successful arrival, a mark is generated using `rmark_func`. The 
#' history of events (times and marks) is passed as separate vectors to the 
#' intensity function to ensure high performance.
#'
#' @param T_max A non-negative numeric value - the end of the interval 
#'   \eqn{[0,T_{max}]}.
#' @param intensity_func A function that calculates the conditional intensity. 
#'   It must accept arguments \code{t} (current time), \code{times} (vector of 
#'   previous event times), and \code{marks} (vector of previous marks).
#' @param rmark_func A function with no required arguments that returns a 
#'   single value representing a mark (e.g., \code{function() rexp(1, 1)}).
#' @param ... Additional arguments passed to \code{intensity_func}.
#'
#' @return An object of class \code{c("marked_pp_sim", "point_process_sim")}. 
#' A list containing:
#' \itemize{
#'   \item \code{events}: A sorted numeric vector of arrival times.
#'   \item \code{marks}: A numeric or character vector of marks associated with each event.
#'   \item \code{T_max}: The end of the simulation window.
#'   \item \code{params}: A list of the functions and parameters used for simulation.
#'   \item \code{n}: The total number of events.
#'   \item \code{type}: The string "Marked Hawkes".
#' }
#' @export
#'
#' @examples
#' # Simple marked process: intensity depends on the sum of previous marks
#' my_intens <- function(t, times, marks, lambda0) {
#'   if (length(times) == 0) return(lambda0)
#'   lambda0 + sum(marks * exp(-(t - times)))
#' }
#' 
#' sim_mhp(T_max = 5, 
#'         intensity_func = my_intens, 
#'         rmark_func = function() runif(1, 0, 1), 
#'         lambda0 = 1)
sim_mhp <- function(T_max, intensity_func, rmark_func, ...) {
  
  # Guardrails
  stopifnot(
    "'T_max' must be positive" = is.numeric(T_max) && length(T_max) == 1 && T_max > 0,
    "'intensity_func' must be a function" = is.function(intensity_func),
    "'rmark_func' must be a function" = is.function(rmark_func)
  )
  
  # Pre-allocate buffers (8KB start)
  times_buf <- numeric(1024)
  marks_buf <- numeric(1024)
  t <- 0
  count <- 0
  
  while (t < T_max) {
    # Extract current history
    h_times <- if (count == 0) numeric(0) else times_buf[1:count]
    h_marks <- if (count == 0) numeric(0) else marks_buf[1:count]
    
    # Find upper bound (intensity at current time t)
    M <- intensity_func(t = t, times = h_times, marks = h_marks, ...)
    
    if (is.na(M)) stop("Intensity calculation returned NA.")
    
    # Generate next candidate point
    t <- t + stats::rexp(1, M)
    if (t > T_max) break
    
    # Check actual intensity at candidate time
    current_int <- intensity_func(t = t, times = h_times, marks = h_marks, ...)
    
    if (stats::runif(1) <= current_int / M) {
      count <- count + 1
      # Expand buffers if needed
      if (count > length(times_buf)) {
        times_buf <- c(times_buf, numeric(length(times_buf)))
        marks_buf <- c(marks_buf, numeric(length(marks_buf)))
      }
      times_buf[count] <- t
      marks_buf[count] <- rmark_func()
    }
  }
  
  .new_point_process_sim(
    events = times_buf[seq_len(count)],
    marks  = marks_buf[seq_len(count)],
    T_max  = T_max,
    params = list(intensity = intensity_func, rmark = rmark_func, dots = list(...)),
    type   = "Marked Hawkes (iid marks)"
  )
}

#' @export
plot.marked_pp_sim <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  
  plot.point_process_sim(x)
  
  plot(x$events, x$marks, type = "h", 
       xlim = c(0, x$T_max), ylim = c(0, max(x$marks, na.rm = TRUE)),
       xlab = "Time (t)", ylab = "Mark (m)",
       main = "Mark Magnitudes over Time", 
       )
  
  points(x$events, x$marks, pch = 21, bg = "slategrey", col = "black", cex = 0.8)
  
  invisible(x)
}


#' @export
print.marked_pp_sim <- function(x, ...) {
  NextMethod("print")
  
  if (x$n > 0) {
    cat("Marks:  ")
    if (is.numeric(x$marks)) {
      m_range <- range(x$marks, na.rm = TRUE)
      cat("Numeric (Range: [", round(m_range[1], 2), ", ", round(m_range[2], 2), "])\n", sep = "")
      cat("        First few: ", paste(round(utils::head(x$marks, 3), 3), collapse = ", "), "...\n")
    } else {
      cat(class(x$marks)[1], "\n")
      cat("        First few: ", paste(utils::head(x$marks, 3), collapse = ", "), "...\n")
    }
  }
  cat("---\n")
  invisible(x)
}
