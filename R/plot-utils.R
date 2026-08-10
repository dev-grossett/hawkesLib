#' Create a Base R Canvas
#' @description Sets up a standard plot area using default base R styling.
#' @param xlim Numeric vector of length 2 (e.g., c(0, T)).
#' @param ylim Numeric vector of length 2 (e.g., c(0, max_y)).
#' @param title String for the plot title.
#' @param xlab String for the x-axis label. Default "Time (t)".
#' @param ylab String for the y-axis label. Default intensity symbol.
#' 
#' @export
plot_base_canvas <- function(
    xlim, 
    ylim, 
    title, 
    xlab = "Time (t)", 
    ylab = expression(lambda^"*"*(t))
) {
  plot(0, type = "n", xlim = xlim, ylim = ylim, 
       xlab = xlab, ylab = ylab, 
       main = title)
}

#' Add Intensity Line to Plot
#' @description Layers an intensity curve onto an existing plot.
#' @param t Numeric vector of time points (the grid).
#' @param intensities Numeric vector of intensity values.
#' @param col Line colour. Default is base R's default black.
#' 
#' @export
add_intensity <- function(t, intensities, col = "black") {
  graphics::lines(t, intensities, col = col)
}

#' Add Event Ticks to Plot
#' @description Adds rug ticks representing event times.
#' @param H_t Numeric vector of event timestamps.
#' @param col Tick colour. Default is base R's default red.
#' 
#' @export
add_events <- function(H_t, col = "red") {
  graphics::rug(H_t, col = col)
}

#' Add Counting Process (Step Plot)
#' @description Layers a step-function representing the counting process N(t).
#' @param H_t Numeric vector of event timestamps.
#' @param T_max The end of the observation window.
#' @param col Colour for the step line. Default is base R's default black.
#' 
#' @export
add_counting_process <- function(H_t, T_max, col = "black") {
  # We start at (0,0), then jump at each event time
  # x-coords: 0 followed by the event times
  # y-coords: 0 followed by 1, 2, ..., n
  x_vals <- c(0, H_t, T_max)
  y_vals <- c(0, seq_along(H_t), length(H_t))
  
  # type = "s" creates the step starting from the left
  graphics::lines(x_vals, y_vals, type = "s", col = col)
}

#' @keywords internal
.get_exp_intensity_grid <- function(t_grid, H_t, theta) {
  # Unpack
  lambda <- theta[1]
  alpha  <- theta[2]
  beta   <- theta[3]
  n      <- length(H_t)
  
  # exponential kernel speedup - same as in .exp_hp_intensity_at_events() 
  A <- numeric(n)
  if (n > 1) {
    for (i in 2:n) {
      A[i] <- exp(-beta * (H_t[i] - H_t[i-1])) * (1 + A[i-1])
    }
  }
  
  # map grid points to the most recent event index
  # idx[j] tells us which event in H_t happened just before t_grid[j]
  idx <- findInterval(t_grid, H_t)
  
  intensities <- numeric(length(t_grid))
  
  # points before the first event are just the baseline
  intensities[idx == 0] <- lambda
  
  # points after at least one event:
  # lambda(t) = lambda + alpha * (1 + A_prev) * exp(-beta * (t - t_prev))
  has_past <- idx > 0
  if (any(has_past)) {
    last_ev_idx <- idx[has_past]
    dt <- t_grid[has_past] - H_t[last_ev_idx]
    
    intensities[has_past] <- lambda + alpha * (1 + A[last_ev_idx]) * exp(-beta * dt)
  }
  return(intensities)
}
