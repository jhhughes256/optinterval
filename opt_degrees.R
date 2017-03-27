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
  spline.degree.func <- function(par, Cobs, time, knots) {
    param <- list(
      X = time,
      Y = Cobs,
      K = knots,
      n = par
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
