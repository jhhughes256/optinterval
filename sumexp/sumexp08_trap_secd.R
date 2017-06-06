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

# error check
  #plot(1:24, pred.sumexp(c(-0.5, -0.05, 6, 5), 1:24, 0))
  #plot(1:24, exp(-0.5*(1:24) + 6) + exp(-0.05*(1:24) + 5))

  #plot(1:72, pred.sumexp(c(-0.5, -0.05, 6, 5), 1:72, 2))
  #plot(1:72, 0.25*exp(-0.5*(1:72) + 6) + 0.0025*exp(-0.05*(1:72) + 5))

  #plot(1:72, pred.sumexp(c(-0.1, -0.2, 4), 1:72, 2))
  #plot(1:72, 0.01*exp(-0.1*(1:72) + 4) - 0.04*exp(-0.2*(1:72) + 4))

# Interval optimising function
  optinterval <- function(times, fit.par) {
    x <- times[order(times)]
    init.par <- x[-c(1, length(x))]
    xmin <- min(x)
    xmax <- max(x)
    theta <- c(NULL)
    for (i in 1:(length(times)-1)) {
      theta[i] <- optim(
        mean(c(times[i], times[i+1])),
        function(t, x) -abs(pred.sumexp(x, t, 2)),
        method = "L-BFGS-B", control = c(maxit = 500),
        lower = times[i], upper = times[i+1], x = fit.par
      )$par
    }
    res <- optim(
      init.par,
      trap.err,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = fit.par, tmin = xmin, tmax = xmax, theta = theta
    )
  # Useful values to return are $par
  # $value is a sum of errors, errors would be more interesting when separated
    return(res)
  }

# Error function for optimisation
# Can either determine second derivative at parameter times, or choose the
# largest absolute second derivative between the two times
  trap.err <- function(par, exp.par, tmin, tmax, theta) {
    times <- c(tmin, par, tmax)
    deltat <- diff(times)
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
  #optinterval(time.seq, c(-0.1, -0.2, 4))      # absorp example

# -----------------------------------------------------------------------------
# Benchmarking:

# 2 compartment
## sumexp05 code
#  test                                                    replications elapsed
#  optinterval(time.seq, c(-0.5, -0.05, 6, 5), absorp = F)          100     4.6
## sumexp08 code
#  test                                                    replications elapsed
#  optinterval(time.seq, c(-0.5, -0.05, 6, 5))                      100    3.23

# absorption
## sumexp05 code
#  test                                                    replications elapsed
#  1 optinterval(time.seq, c(-0.1, -0.2, 4, 4), absorp = T)         100    7.05
## sumexp08 code
#  test                                                    replications elapsed
#  1 optinterval(time.seq, c(-0.1, -0.2, 4))                        100    2.48

# Some time shaved off with 2 compartment curves
# but alot of time removed from absorption curves
