# Updated function to determine deltat's place in the error formula
# Aim: Determine if an optim routine should be used to find the maximum secd
# -----------------------------------------------------------------------------
# Incorporates functions made in sumexp06 & sumexp07
# Verbose comments can be found in sumexp05
# Questions that remain:
#  - Does the error sum method matter?
#  - Does this apply to both "linear-up log-down" and "linear-up and -down?"
# -----------------------------------------------------------------------------
# pred.sumexp taken from "sumexp07..."
  pred.sumexp <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    for (i in 1:n) {
      if (i == 1) y <- x[i]^d*exp(x[i]*t + x[n+i])
      else if (i != n | a == 0) y <- y + x[i]^d*exp(x[i]*t + x[n+i])
      else if (a == 1) y <- y - x[i]^d*exp(x[i]*t)*sum(exp(x[(n+1):(2*n-1)]))
    }
    return(y)
  }

# Interval optimising function
  optinterval <- function(time, fit.par) {
    x <- time[order(time)]
    init.par <- x[-c(1, length(x))]
    xmin <- min(x)
    xmax <- max(x)
    res <- optim(
      init.par,
      trap.err,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = fit.par, tmin = xmin, tmax = xmax
    )
  # Useful values to return are $par
  # $value is a sum of errors, errors would be more interesting when separated
    return(res)
  }

# Error function for optimisation
# Can either determine second derivative at parameter times, or choose the
# largest absolute second derivative between the two times
  trap.err <- function(par, exp.par, tmin, tmax, absorp) {
    times <- c(tmin, par, tmax)
    deltat <- diff(times)
    theta <- c(NULL)
    for (i in 1:length(deltat)) {
      theta[i] <- optim(
        mean(c(times[i], times[i+1])),
        function(t, x) pred.sumexp(x, t, 2),
        method = "L-BFGS-B", control = c(maxit = 500),
        lower = times[i], upper = times[i+1], x = exp.par
      )$par
    }
    browser()
    secd <- pred.sumexp(exp.par, theta, 2)  # d = 2 (second derivative)
    err <- abs(deltat^3*secd/12)
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
  #optinterval(time.seq, c(-0.5, -0.05, 6, 5))  # 2 comp example
  #optinterval(time.seq, c(-0.2, -0.1, 4))      # absorp example

# -----------------------------------------------------------------------------
