# Initial Testing Script for Optimising Splines to Fit Data
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	set.seed(123)

# Load libraries
  library(dplyr)
  library(ggplot2)
  library(splines)

# ------------------------------------------------------------------------------
# Define some functions that will create data
	simulate.1comp.abs <- function(x) {
		k10 <- x$CL/x$V2
		A <- x$KA/(x$V2*(x$KA - k10))
		x$AMT*A*(exp(-k10*x$TIME) - exp(-x$KA*x$TIME))
	}

	simulate.2comp.abs <- function(x) {
		k10 <- x$CL/x$V2
		k12 <- x$Q/x$V2
		k21 <- x$Q/x$V3
		apb <- k10+k12+k21            # alpha + beta
		amb <- k10*k21                # alpha * beta
		alpha <- ((apb)+sqrt((apb)^2-4*amb))/2
		beta <- ((apb)-sqrt((apb)^2-4*amb))/2
		A <- x$KA*(k21-alpha)/(x$V2*(x$KA-alpha)*(beta-alpha))
		B <- x$KA*(k21-beta)/(x$V2*(x$KA-beta)*(alpha-beta))
		x$AMT*(A*exp(-alpha*x$TIME)+B*exp(-beta*x$TIME)-(A+B)*exp(-x$KA*x$TIME))
	}
# ------------------------------------------------------------------------------
# Create some "observed data"
  ERR <- 0.3
  times <- c(0:24)
  sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)  # Sampling times

  obs.param1 <- list(
    AMT = 50,  # mg
    CL = 10,  # Clearance, 10 L/h
    V2 = 50,  # Volume of distribution, 50 L
    KA = 0.5,  # Absorption rate constant, h^-1
    TIME = seq(from = 0,to = 24,by = 0.25)  # Time sequence for simulation
  )
  error1 <- 1 + rnorm(n = length(obs.param1$TIME),mean = 0,sd = ERR)
  conc1 <- simulate.1comp.abs(obs.param1) * error1
  sample.conc1 <- conc1[obs.param1$TIME %in% times]
  plot(sample.conc1 ~ times, col = "red")
  testdf1 <- data.frame(time = times, Cp = sample.conc1)

  obs.param2 <- list(
    AMT = 50,  # mg
    CL = 10,  # Clearance, 10 L/h
    V2 = 50,  # Volume of distribution, 50 L
    Q = 8,  #
    V3 = 100,  #
    KA = 0.5,  # Absorption rate constant, h^-1
    TIME = seq(from = 0,to = 24,by = 0.25)  # Time sequence for simulation
  )
  error2 <- 1 + rnorm(n = length(obs.param2$TIME),mean = 0,sd = ERR)
  conc2 <- simulate.2comp.abs(obs.param2) * error2
  sample.conc2 <- conc2[obs.param2$TIME %in% times]
  plot(sample.conc2 ~ times, col = "red")
  testdf2 <- data.frame(time = times, Cp = sample.conc2)
# ------------------------------------------------------------------------------
# Fit a spline to the data using stats::spline
  spline.df1 <- testdf1 %>%
    spline(n = 100) %>%
    data.frame() %>%
    rename(time = x, Cp = y)
  plotobj0 <- NULL
  plotobj0 <- ggplot(testdf1, aes(x = time, y = Cp))
  plotobj0 <- plotobj0 + geom_point()
  plotobj0 <- plotobj0 + geom_line(data = spline.df1)
  plotobj0

  spline.fn1 <- testdf1 %>%
    splinefun(method = "fmm", ties = mean)

  test <- data.frame(times = sample.times) %>%
    transform(yhat = spline.fn1(sample.times))

  plotobj0 + geom_point(data = test, aes(x = times, y = yhat), colour = "red")

# The spline here fits too closely, a spline with less knots is required
# This is due to the fact that a single outlier would throw the entire thing off
# As spline() does not have the option to define knots it won't be used
# ------------------------------------------------------------------------------
# Fit a spline to the data using splines::bs
# Works a little differently as instead of producing a spline it produces
#  a basis function instead. Not as easy to understand in practice...
  bs(testdf1$time, df = 4)

# Here is some code that I found
  summary(fm1 <- lm(Cp ~ bs(times, df = 4), data = testdf1))
	predict(fm1, time = times)

# However this does not appear to work as it seems to
  predict(fm1, time = sample.times)
  predict(fm1)

# Further searching has brought light to the situation
  set.seed(1)
  n <- 400
  x <- 0:(n-1)/(n-1)
  f <- 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
  y <- f + rnorm(n, 0, sd = 2)
  mod <- lm(y ~ bs(x, knots = seq(0.1, 0.9, by = 0.1)))
  pdat <- data.frame(x = seq(min(x), max(x), length = 100))
  ## predict for new `x`
  pdat <- transform(pdat, yhat = predict(mod, newdata = pdat))

  ## now plot
  ylim <- range(pdat$y, y) ## not needed, but may be if plotting CIs too
  plot(y ~ x)
  lines(yhat ~ x, data = pdat, lwd = 2, col = "red")

  rm(n, x, f, y, mod, pdat, ylim)

# This coder in this situation has multiple original observations
# He then interpolates within these original observations, grabbing the closest
#  concentration and saying it correlates with the value you gave...
# This doesn't seem ideal
# ------------------------------------------------------------------------------
# I may be approaching this incorrectly. What am I using the splines to achieve?
# It's primary use is to determine the second derivative at any time
# Uneasy about using stats::spline as the second derivative fluctuates a lot
# ------------------------------------------------------------------------------
# Fit a spline using splines::interpSpline
  interpSpline(testdf1$time, testdf1$Cp, bSpline = T)

# Hmmm that's no good either...
