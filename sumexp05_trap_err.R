# Creating a function to determine deltat's place in the error formula
# Conditions: Must have number of t, min and max as input
# -----------------------------------------------------------------------------
# Uses ideas from:
# - error_test2.R by RU
# - Trapezoidal Error Formula
# -----------------------------------------------------------------------------
# Set parameters
  #nobs <- 12
  #tmin <- 0
  #tmax <- 24
# Create time schedule
  #time.seq <- seq(tmin, tmax, by = (tmax-tmin)/nobs)
# don't want optim to optimise the boundaries, tmin and tmax
  #init.par <- time.seq[-c(1, nobs + 1)]
# but you would include those back in, within the optim function
  #time.seq %in% c(tmin, init.par, tmax)
# Determine delta t using the pieced together time sequence
  #deltat <- c(init.par, tmax) - c(tmin, init.par)
# -----------------------------------------------------------------------------
# Edited pred.sumexp function from "sumexp03_..."
# designed to take into account which derivative you desire
# default is 0, should be compatible with all previous scripts
  pred.sumexp <- function(x, t, a, d = 0) {
    n <- length(x)/2
    for (i in 1:n) {
      if (!exists("y")) y <- x[i]^d*exp(x[i]*t + x[n+i])
      else if (i != n | a == F) y <- y + x[i]^d*exp(x[i]*t + x[n+i])
      else if (a == T) y <- y - x[i]^d*exp(x[i]*t + x[n+i])
    }
    return(y)
  }

# Function will take the original data and the function parameters for that data
# time should be the sample points for the data. Function can be altered to
# accept a data.frame if necessary.
  optinterval <- function(time, fit.par, absorp) {
    x <- time[order(time)]
    init.par <- x[-c(1, length(x))]  # remove first and last value (fixed values)
    xmin <- min(x)
    xmax <- max(x)
  # optim control argument used due to large number of iterations
    res <- optim(
      init.par,
      trap.err,
      method = "L-BFGS-B",
      lower = xmin, upper = xmax, control = c(maxit = 500),
      exp.par = fit.par, tmin = xmin, tmax = xmax, absorp = absorp
    )
  # Useful values to return are $par
  # $value is a sum of errors, errors would be more interesting when separated
    return(res)
  }

# Can either determine second derivative at parameter times, or choose the
# largest absolute second derivative between the two times
  trap.err <- function(par, exp.par, tmin, tmax, absorp) {
    times <- c(tmin, par, tmax)  # combine variable and fixed times
    deltat <- diff(times)
    secd <- pred.sumexp(exp.par, times, absorp, 2)  # d = 2 (second derivative)
    err <- abs(deltat^3*secd[-length(secd)]/12)  # trapezoidal error function
  # Might need to play around with different methods for error sum
    sum(err)
  }
# -----------------------------------------------------------------------------
# Set parameters
  nobs <- 12
  tmin <- 0
  tmax <- 24

# Set time sequence that is not far from the truth
# -2.5 is an arbitrary number, but does this relate in anyway to the exponentials?
  time.seq <- c(0, exp(seq(-2.5, 0, length.out = nobs))*tmax)

# Optimise
  optinterval(time.seq, c(-0.5, -0.05, 6, 5), absorp = F)
  optinterval(time.seq, c(-0.1, -0.2, 4, 4), absorp = T)
