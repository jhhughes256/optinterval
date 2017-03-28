# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  library(ggplot2)
  library(plyr)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# -----------------------------------------------------------------------------
# Create a function to create a spline with limited predictions
  predict.spline.func <- function(x) {
    mod <- lm(
      x$Y ~ bs(x$X, knots = x$K, degree = x$n)
    )  # lm
    B <- data.frame(time = x$X)
    Y <- predict(mod, newdata = B)
  }

  uber.predict.spline.func <- function(x) {
    browser()
    mod <- lm(
      x$Y ~ bs(x$X, knots = x$K, degree = x$n)
    )  # lm
    B <- data.frame(time = seq(0, 24, by = 0.25))
    Y <- predict(mod, newdata = B)
  }

# -----------------------------------------------------------------------------
# Utility function to be used with knot optimisation code
  knot.func <- function(x) {
    y <- double(length(x))
    for (i in 1:length(x)) {
      y[i] <- sum(x[1:i])
    }
    return(y)
  }

# -----------------------------------------------------------------------------
# Wrapper that calculates residuals to be minimised
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

# -----------------------------------------------------------------------------
# Set up the knot optimising function
  which.knots <- function(x, y, max.knots, degree) {
    xmax <- x[which(y == max(y))]
    err <- NULL
    knots <- list(NULL)
    warn <- list(NULL)
    for (i in 1:max.knots) {
      par <- rep(xmax, i)
      opt.res <- optim(
        par,  # Initial parameter estimates
        error.spline.func,  # Fitting function
        # method = "SANN",
        method = "L-BFGS-B", lower = min(x), upper = max(x),
        Cobs = y, time = x, degree = degree  # Function arguments
      )
      err[i] <- opt.res$value
      knots[[i]] <- knot.func(opt.res$par)
    }
    list(
      error = err,
      knots = knots
    )
  }

# -----------------------------------------------------------------------------
# Simulate dataset
  K <- c(5, 10)  # knots
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

# Determine if there is a difference between degrees
  spline.data <- llply(seq_len(6), function(x) {
    which.knots(onecomp$time, onecomp$conc, 4, x)
  })

  llply(seq_len(length(spline.data)), function(x) {
    degree <- x
    knots <- spline.data[[x]]$knots
    llply(knots, function(x) {
      param <- list(
        X = onecomp$time,
        Y = onecomp$conc,
        K = x,
        n = degree
      )
      x.vec <- c(x, rep(0, length(knots) - length(x)))
      browser()
      out <- data.frame(
        time = seq(from = 0, to = 24, by = 0.25),
        conc = uber.predict.spline.func(param),
        degree = degree,
        nknots = length(x)
      )
      for (i in length(knots)) {
        out[[paste0("knot", i)]] <- x.vec[i]
      }
      out
    })
  })

# -----------------------------------------------------------------------------
# Simulate second dataset with more rich sampling
  sample.times <- c(0,0.25,0.5,1,2,3,4,6,8,10,12,18,24)	# Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))*(1+rnorm(n = length(times),mean = 0,sd = ERR))
  sample.conc <- conc[times %in% sample.times]
  onecomp <- data.frame(time = sample.times, conc = sample.conc)

# Determine if there is a difference between degrees
  splines <- llply(seq_len(6), function(x) {
    which.knots(onecomp$time, onecomp$conc, 4, x)
  })




# Naturally a 6 degree spline fits the data the best but that should be taken
# with a pinch of salt
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = u, y = Y), size = 1, colour = "red")
  plotobj


# -----------------------------------------------------------------------------
