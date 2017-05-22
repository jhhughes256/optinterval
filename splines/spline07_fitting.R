# Working out how to fit a linear model without using lm()
# Inspired by "Fitting a Model by Maximum Likelihood"
# by Andrew of Exegetic Analytics
# https://www.r-bloggers.com/fitting-a-model-by-maximum-likelihood/

# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

# -----------------------------------------------------------------------------
#  So, consider the following dataset, with the following spline regression,
  library(splines)
  K <- c(14, 20)  # knots
  plot(cars)
  lm.mod <- lm(dist ~ bs(speed, knots = c(K), degree = 2), data = cars)
  u <- seq(4, 25, by = 0.1)  # test model at these speeds
  B <- data.frame(speed = u)  # data.frame(speeds)
  lm.Y <- predict(lm.mod, newdata = B)  # predicted dist for the desired speeds
  lines(u, lm.Y, lwd = 2, col = "red")
# We have the following (nice) picture
# But this is done using the lm() function which uses QR decomposition
# While this method is efficient, it uses least squares

# The aim is to determine how to fit a model using gradient descent so that
# maximum likelihood can be used instead
# So first lets max a basic example to work on
  test.mod <- lm(dist ~ speed, data = cars)
  test.Y <- predict(test.mod, newdata = B)
  plot(cars)
  lines(u, test.Y, lwd = 2, col = "red")
  summary(test.mod)

# This model is a straight line so its represented by the function:
# y = beta1*x + beta0
# This is effectively a 1 degree spline with no knots
# lm() looks to find the model with the best residuals, using least squares
# e.g. sum((C[obs] - C[hat])^2)
# The same can be done with optim()

# function(Cobs, Chat) {
#   err <- Cobs - Chat
#   squ.err <- err^2
#   sse <- sum(squ.err)  # Sum of squared errors to be minimised
#   sse
# }

# This final values is what we want to optimise
# To reach the lowest possible squared error, we request that optim changes
# beta1 and beta0 of our function

  ls.func <- function(par, x, y) {
    b0 <- par[1]
    b1 <- par[2]
    err <- y - b1*x - b0
    sse <- sum(err^2)
  }

  ls.func.result <- optim(
    c(-17, 4),  # Initial parameter estimates
    ls.func,  # Fitting function
    method = "BFGS",
    x = cars$speed, y = cars$dist  # Function arguments
  )

# Here we appear to have replicated lm() to within four significant figures
# To make predictions using this model we simply replicate the function
  plot(cars)
  ls.func.Y <- ls.func.result$par[2]*u + ls.func.result$par[1]
  lines(u, ls.func.Y, lwd = 2, col = "red")

# The trick is to apply this to splines
# -----------------------------------------------------------------------------
# Splines are represented by an additive linear model, where each component of
# the spline does not interact with one another.

# To determine the number of parameters for optimisation, the degree of the
# spline and the number of knots must be considered. These are added together to
# determine the number of components there are to the linear model
# 1 degree + 0 knots = 1 component = 2 parameters
# 2 degree + 2 knots = 4 components = 5 parameters

# The basis matrix for a spline is used when fitting a model, therefore we use
# each part of the matrix as an independent variable to determine our dependent
# variable.

# Function for a 2 degree spline with 2 knots
# y = beta1*x[,1] + beta2*x[,2] + beta3*x[,3] + beta4*x[,4] + beta0
  ls.spline2d2k.func <- function(par, x, y) {
    b0 <- par[1]
    b1 <- par[2]
    b2 <- par[3]
    b3 <- par[4]
    b4 <- par[5]
    err <- y - b1*x[,1] - b2*x[,2] - b3*x[,3] - b4*x[,4] - b0
    sse <- sum(err^2)
  }

  ls.spline2d2k.res <- optim(
    c(6, 8, 45, 49, 95),  # Initial parameter estimates
    ls.spline2d2k.func,  # Fitting function
    method = "BFGS",
    x = bs(cars$speed, knots = c(K), degree = 2), y = cars$dist # Function arg.
  )

# To predict a model with a spline in it the new x value must be transformed
# into a basis spline matrix, as dictated by the spline we were fitting
  plot(cars)
  ls.spline.u <- bs(u, knots = c(K), degree = 2)
  ls.spline.Y <- ls.spline2d2k.res$par[2]*ls.spline.u[,1] +
    ls.spline2d2k.res$par[3]*ls.spline.u[,2] +
    ls.spline2d2k.res$par[4]*ls.spline.u[,3] +
    ls.spline2d2k.res$par[5]*ls.spline.u[,4] + ls.spline2d2k.res$par[1]
  lines(u, ls.spline.Y, lwd = 2, col = "red")

# -----------------------------------------------------------------------------
# So how do we make a flexible function to complete this process?
# Well first we should move around our parameters as determined by optim()
# This should be done so that for:
# par[1] == b[1], par[2] == b[2], ..., par[n-1] == b[last], par[n] == intercept
# Where: n = number of parameters
# That way the basis matrix will match up with the numbering of the parameters
  flex.u <- bs(u, knots = c(K), degree = 2)
  flex.optres <- list(
    par = c(
      ls.spline2d2k.res$par[2],
      ls.spline2d2k.res$par[3],
      ls.spline2d2k.res$par[4],
      ls.spline2d2k.res$par[5],
      ls.spline2d2k.res$par[1]
    )
  )
  flex.optres.l <- length(flex.optres$par)
  flex.x <- bs(cars$speed, knots = c(K), degree = 2)

# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# For fitting
# x = bs(cars$speed, knots = c(K), degree = 2), y = cars$dist
  flex.err <- cars$dist
  for (j in 1:flex.optres.l) {
    if (j != flex.optres.l) {
      flex.err <- flex.err - flex.optres$par[i]*flex.x[,i]
    } else {
      flex.err <- flex.err - flex.optres$par[i]
    }
  }
  flex.sse <- sum(flex.err^2)
  flex.sse != ls.spline2d2k.res$value

  lm.spline <- function(par, x, y) {
    err <- y
    for (i in 1:length(par)) {
      if (i != length(par)) {
        err <- err - par[i]*x[,i]
      } else if (i == length(par)) {
        err <- err - par[i]
      }
    }
    return(sum(err^2))
  }

  test <- optim(
    c(8, 45, 49, 95, 6),  # Initial parameter estimates
    lm.spline,  # Fitting function
    method = "BFGS",
    x = bs(cars$speed, knots = c(K), degree = 2), y = cars$dist # Function arg.
  )
  test$value != ls.spline2d2k.res$value

# -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
# For predicting
  flex.optres.l <- length(flex.optres$par)
  flex.Y <- double(length(u))
  for (j in 1:flex.optres.l) {
    if (j != flex.optres.l) {
      flex.Y <- flex.Y + flex.optres$par[i]*flex.u[,i]
    } else {
      flex.Y <- flex.Y + flex.optres$par[i]
    }
  }
  any(flex.Y != ls.spline.Y)  # Test that this has worked

# For the following function:
# x = the new x values for which you want to predict y
# par = the coefficients of the spline you have fitted
# K = the knots used to fit the spline
# n = the degree used to fit the spline
  predict.spline <- function(x, par, K, n) {
    bs.mat <- bs(x, knots = c(K), degree = n)
    out <- double(length(x))
    for (i in 1:length(par)) {
      if (i != length(par)) {
        out <- out + par[i]*bs.mat[,i]
      } else if (i == length(par)) {
        out <- out + par[i]
      }
    }
    return(out)
  }
# -----------------------------------------------------------------------------
# Lastly lets adjust these functions so that they can be used with maximum
# likelihood estimation and see how much of a difference it makes for the spline
# The function that is given to optim() needs to also optimise sigma so we add
# another parameter
  mle.spline <- function(par, x, y) {
    err <- y
    for (i in 1:length(par)) {
      if (i < length(par) - 1) {
        err <- err - par[i]*x[,i]
      } else if (i == length(par) - 1) {
        err <- err - par[i]
      }
    }
    loglik <- dnorm(err, 0, tail(par, 1), log = T)
    return(-1*sum(loglik))
  }

  mle.test.res <- optim(
    c(50, 2, 8, 60, 35, 1),  # Initial parameter estimates
    mle.spline,  # Fitting function
    method = "BFGS",
    x = bs(cars$speed, knots = c(K), degree = 2), y = cars$dist # Function arg.
  )

  mle.test.Y <- predict.spline(u, mle.test.res$par[-length(mle.test.res$par)], K, 2)
  plot(cars)
  lines(u, mle.test.Y, lwd = 2, col = "red")

# This isn't perfect though if you don't have good initial estimates

  mle.test.res2 <- optim(
    c(50, 2, 8, 60, 35, 1),  # Initial parameter estimates
    mle.spline,  # Fitting function
    method = "BFGS",
    x = bs(cars$speed, knots = c(K), degree = 2), y = cars$dist # Function arg.
  )

  mle.test.Y <- predict.spline(u, mle.test.res2$par[-length(mle.test.res2$par)], K, 2)
  plot(cars)
  lines(u, mle.test.Y, lwd = 2, col = "red")

# Splines appear to be highly dependent on initial estimates!
# But if we have to run lm() to run mle.spline() then why run mle.spline()?
# For this question we need to start looking at PK examples
# -----------------------------------------------------------------------------
  times <- seq(from = 0, to = 24, by = 0.25)	# Time sequence for simulating concentrations

# Simulate dataset
  K <- c(5, 10)  # knots
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  ERR <- 1 + rnorm(n = length(times), mean = 0, sd = 0.3)
  dose <- 50	# mg
  sample.times <- c(0, 0.25, 0.5, 1, 2, 4, 8, 12, 24)  # Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))
  plot(times, log(conc))
  sample.conc <- (conc*ERR)[times %in% sample.times]
  onecomp <- data.frame(time = sample.times, conc = sample.conc)

# Find our initial estimates using QR decomposition via lm()
  pk.modlm <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
  init.par <- unname(pk.modlm$coefficients[c(2:5, 1)])

# Supress warnings due to NaN being produced by estimating sigma of < 0
  suppressWarnings(
    pk.optres <- optim(
      c(init.par, 1),  # Initial parameter estimates
      mle.spline,  # Fitting function
      method = "BFGS",
      x = bs(onecomp$time, knots = c(K), degree = 2), y = onecomp$conc # Function arg.
    )
  )

  u <- times
  B <- data.frame(time = u)
  pk.lm.Y <- predict(pk.modlm, newdata = B)
  pk.mle.Y <- predict.spline(u, pk.optres$par[-length(pk.optres$par)], K, 2)
  plot(onecomp)
  lines(u, pk.lm.Y, lwd = 2, col = "red")
  lines(u, pk.mle.Y, lwd = 2, col = "blue")

# The two lines aren't visibly different when plotted
# When checking with least squares lm is better when comparing to original data
  list(
    lm = sum((conc - pk.lm.Y)^2),
    mle = sum((conc - pk.mle.Y)^2)
  )

# And mle is slightly better when looking at mle
  list(
    lm = -1*sum(dnorm(conc, pk.lm.Y, sd(conc - pk.lm.Y), log = T)),
    mle = -1*sum(dnorm(conc, pk.mle.Y, sd(conc - pk.lm.Y), log = T))
  )

# This phenomemon is the same when observing drugs with two-compartment kinetics
# Will use mle going forwards as it has a big difference when optimising knots
# Best to keep with one method rather than using two
