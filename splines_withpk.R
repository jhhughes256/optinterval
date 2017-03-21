# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  library(ggplot2)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

  K <- c(5, 10)  # knots
  # NOT onecomp - > PK! WIP
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  ERR <- 0.3	# Standard deviation of error
  dose <- 50	# mg
  times <- seq(from = 0,to = 24,by = 0.25)	# Time sequence for simulating concentrations
  sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)	# Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))*(1+rnorm(n = length(times),mean = 0,sd = ERR))
  sample.conc <- conc[times %in% sample.times]
  onecomp <- data.frame(time = sample.times, conc = sample.conc)

# -----------------------------------------------------------------------------
# Interesting plot that shows where the best place to put a knot would be
  vk <- seq(0.05, 0.95, by = 0.05)
  SSR <- matrix(NA, length(vk))
  for(i in 1:(length(vk))){
    k <- vk[i]
    K <- min(onecomp$time) + k*(max(onecomp$time) - min(onecomp$time))
    reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
    SSR[i] <- sum(residuals(reg)^2)
  }
  dat <- data.frame(vk, c(SSR))
  plotobj <- ggplot(dat, aes(x = vk, y = SSR))
  plotobj <- plotobj + geom_point(shape = 21, size = 2, colour = "blue")
  plotobj <- plotobj + geom_line(linetype = 2, colour = "blue")
  plotobj

# ------------------------------------------------------------------------------
# Create a function to create a spline
  predict.spline.func <- function(x) {
    mod <- lm(
      x$Y ~ bs(x$X, knots = x$K, degree = x$n)
    )  # lm
    B <- data.frame(time = x$X)
    Y <- predict(mod, newdata = B)
  }

# -----------------------------------------------------------------------------
# Wrapper that calculates residuals to be minimised
  spline.oneknot.func <- function(par, Cobs, time, degree) {
    param <- list(
      X = time,
      Y = Cobs,
      K = par[1],
      n = degree
    )
    Chat <- predict.spline.func(param)
    err <- Cobs-Chat
    squ.err <- err^2
    sse <- sum(squ.err)  # Sum of squared errors to be minimised
    sse
  }

# Initial estimates
  tmax <- onecomp$time[which(onecomp$conc == max(onecomp$conc))]
	init.par <- c(tmax)
  #par <- init.par

# Minimise the sum of squared errors
# Using the optim function in R
	result <- optim(
		init.par,  # Initial parameter estimates
		spline.oneknot.func,  # Fitting function
    # method = "SANN",
    # method = "L-BFGS-B", lower = 0, upper = 24,
    method = "Brent", lower = 0, upper = 24,
		Cobs = onecomp$conc, time = onecomp$time, degree = 2  # Function arguments
	)
	result  # Print the result

# -----------------------------------------------------------------------------
# When using two knots, there are two ways to achieve this
# Have optim determine the two numbers
# Or determine the first knot, and then the distance of the second knot from
# that knot

# Wrapper using the distance from last knot method
  lk.spline.twoknot.func <- function(par, Cobs, time, degree) {
    param <- list(
      X = time,
      Y = Cobs,
      K = c(par[1], par[1] + par[2]),
      n = degree
    )
    Chat <- predict.spline.func(param)
    err <- Cobs-Chat
    squ.err <- err^2
    sse <- sum(squ.err)  # Sum of squared errors to be minimised
    sse
  }

# Wrapper using unique unrelated parameters
  uq.spline.twoknot.func <- function(par, Cobs, time, degree) {
    param <- list(
      X = time,
      Y = Cobs,
      K = c(par[1], par[2]),
      n = degree
    )
    Chat <- predict.spline.func(param)
    err <- Cobs-Chat
    squ.err <- err^2
    sse <- sum(squ.err)  # Sum of squared errors to be minimised
    sse
  }

# Initial estimates
  tmax <- onecomp$time[which(onecomp$conc == max(onecomp$conc))]
	init.par <- c(tmax, tmax)
  #par <- init.par

# Minimise the sum of squared errors
# Using the optim function in R
	result <- optim(
		init.par,  # Initial parameter estimates
		lk.spline.twoknot.func,  # Fitting function
    # method = "SANN",
    method = "L-BFGS-B", lower = 0, upper = 24,
		Cobs = onecomp$conc, time = onecomp$time, degree = 2  # Function arguments
	)
	result  # Print the result

  init.par <- c(tmax, tmax*2)

  result <- optim(
		init.par,  # Initial parameter estimates
		uq.spline.twoknot.func,  # Fitting function
    # method = "SANN",
    method = "L-BFGS-B", lower = 0, upper = 24,
		Cobs = onecomp$conc, time = onecomp$time, degree = 2  # Function arguments
	)
	result  # Print the result

# The differences between the two methods seems negligible, but this should be
# stress tested when all other functions are stress tested

# Further development will be done using the last knot method as the unique
# method is easier to reimplement later
# -----------------------------------------------------------------------------
# Deprecated:
#spline.threeknot.func <- function(par, Cobs, time, degree) {
#  param <- list(
#    X = time,
#    Y = Cobs,
#    K = c(par[1], par[1] + par[2], par[1] + par[2] + par[3]),
#    n = degree
#  )
#  Chat <- predict.spline.func(param)
#  err <- Cobs-Chat
#  squ.err <- err^2
#  sse <- sum(squ.err)  # Sum of squared errors to be minimised
#  sse
#}
# Three knots seems to be too many, but might be ok in more complex models?
# Doubtful... 2 knots is already alot of complexity for 7 internal points
# Is 3 knots ok if there are more internal points?
# Could this be calculated to form a rule for the max number of knots tested
# -----------------------------------------------------------------------------
# So we need a function that adds each value with the value before in case we
# ever need more than two knots.
# Also means that the oneknot and twoknot variants of the function can be
# combined into one uber-func

# Takes a vector x[1], x[2], ..., x[length(x)]
# Calculates sum(x[1]),  sum(x[1], x[2]), ...
  knot.func <- function(x) {
    y <- double(length(x))
    for (i in 1:length(x)) {
      y[i] <- sum(x[1:i])
    }
    return(y)
  }
# -----------------------------------------------------------------------------
# Using this knot.func (name pending) the functions above should be able to be
# combined into one function.
  error.spline.func <- function(par, Cobs, time, degree) {
    param <- list(
      X = time,
      Y = Cobs,
      K = knot.func(par),
      n = degree
    )
    Chat <- predict.spline.func(param)
    err <- Cobs-Chat
    squ.err <- err^2
    sse <- sum(squ.err)  # Sum of squared errors to be minimised
    sse
  }

# Initial estimates
  tmax <- onecomp$time[which(onecomp$conc == max(onecomp$conc))]
	init.par <- c(tmax, tmax)
  #par <- init.par

# Minimise the sum of squared errors
# Using the optim function in R
	result <- optim(
		init.par,  # Initial parameter estimates
		error.spline.func,  # Fitting function
    # method = "SANN",
    method = "L-BFGS-B", lower = 0, upper = 24,
		Cobs = onecomp$conc, time = onecomp$time, degree = 2  # Function arguments
	)
  print(K <- c(knot = knot.func(result$par)))
  result$value  # Print the result

# Success!
# -----------------------------------------------------------------------------
# Next we want to try to compare using varying amounts of knots...
# This would ideally be determined by the number of internal points
# A good roundabout figure would be 1 point for every 3.5 points?
# i.e > 3 points - 1 knot
#     > 6 points - 2 knots
#     > 10 points - 3 knots

# An issue here is there is currently no control for points being too close
# to one another, ideally knots would be at least.. 3 hours apart?

# Solution: First knot has standard boundaries c(0,24)
#           Successive knots have the previous knot as a lower boundary
#           and 24 - the previous knot for the upper boundary?

# This is not possible within the optim function, but a gradient could be
# designed in such a way that parameter values that break these rules cause a
# large malice on the objective function value

# This is something to address in the future, for now we compare varying knots
# -----------------------------------------------------------------------------
# Limitiation of current function:
#  - if xmax is more than 1/max.knots of max(x) then par will exceed the bounds
#  - Excessive warnings appear when multiple knots are between internal points
#  - this is due to a rank-deficient fit, I think this is because it's probably
#    not enough information to have two knots between two internal points
  which.knots <- function(x, y, max.knots, degree) {
    xmax <- x[which(y == max(y))]
    err <- NULL
    knots <- list(NULL)
    for (i in 1:max.knots) {
      par <- rep(xmax, i)
      opt.res <- optim(
        par,  # Initial parameter estimates
        error.spline.func,  # Fitting function
        # method = "SANN",
        method = "L-BFGS-B", lower = min(x), upper = max(x),
        Cobs = y, time = x, degree = 2  # Function arguments
      )
      err[i] <- opt.res$value
      knots[[i]] <- knot.func(opt.res$par)
    }
    knots[[which(err == min(err))]]
  }
  K <- which.knots(onecomp$time, onecomp$conc, 2, 2)

  reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
  u <- seq(0, 24, by = 0.1)  # test model at these times
  B <- data.frame(time = u)  # data.frame(times)
  Y <- predict(reg, newdata = B)  # predicted conc for the desired times
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = u, y = Y), size = 1, colour = "red")
  plotobj

# The values assigned to Y are what we will be carrying forwards to determine
# the best sampling times, but first, what about degrees?
# -----------------------------------------------------------------------------
