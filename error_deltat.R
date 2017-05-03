# Creating a function to determine deltat's place in the error formula
# Conditions: Must have number of t, min and max as input
# -----------------------------------------------------------------------------
# Uses ideas from:
# - error_test2.R by RU
# - Trapezoidal Error Formula
# -----------------------------------------------------------------------------
# Set parameters
  nobs <- 12
  tmin <- 0
  tmax <- 24

# Create time schedule
  time.seq <- seq(tmin, tmax, by = (tmax-tmin)/nobs)

# don't want optim to optimise the boundaries, tmin and tmax
  init.par <- time.seq[-c(1, nobs + 1)]

# but you would include those back in, within the optim function
  time.seq %in% c(tmin, init.par, tmax)

# Determine delta t using the pieced together time sequence
  deltat <- c(init.par, tmax) - c(tmin, init.par)

# Load sumexp function
  secd.sumexp <- function(x, t, subtract.last) {
    n <- length(x)/2
    for (i in 1:n) {
      if (!exists("secd")) secd <- x[i]^2*exp(x[i]*t + x[n+i])
      else if (i != n | subtract.last == F) secd <- secd + x[i]^2*exp(x[i]*t + x[n+i])
      else if (subtract.last == T) secd <- secd - x[i]^2*exp(x[i]*t + x[n+i])
    }
    return(secd)
  }

# Function will take the original data and the function parameters for that data
# Data should be the sample points for the data. Function can be altered to
# accept a data.frame if necessary.
  optinterval <- function(time, fit.par, absorp) {
    x <- time[order(time)]
    init.par <- x[-c(1, length(x))]
    xmin <- min(x)
    xmax <- max(x)
    optim(
      init.par,
      trap.err,
      method = "L-BFGS-B",
      lower = xmin, upper = xmax,
      exp.par = fit.par, tmin = xmin, tmax = xmax, absorp = absorp
    )
  }

# Can either determine second derivative at parameter times, or choose the
# largest absolute second derivative between the two times
  trap.err <- function(par, exp.par, tmin, tmax, absorp) {
    times <- c(xmin, par, xmax)
    deltat <- diff(times)
    secd <- secd.sumexp(exp.par, times, absorp)
    err <- abs(deltat^3*secd[-1]/12)
    # there may be some issues with the collation of error...
    sum(err)
  }
