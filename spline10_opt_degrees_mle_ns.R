# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  library(ggplot2)
  library(plyr)
  library(stringr)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# -----------------------------------------------------------------------------
# Create a function to create a spline with limited predictions
  predict.spline.func <- function(x) {
    time <- x$X
    conc <- x$Y
    mod <- lm(
      conc ~ ns(time, knots = x$K)
    )  # lm
    B <- data.frame(time = x$X)
    Y <- predict(mod, newdata = B)
  }

  uber.predict.spline.func <- function(x) {
    time <- x$X
    conc <- x$Y
    mod <- lm(
      conc ~ ns(time, knots = x$K)
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
  error.spline.func <- function(par, Cobs, time) {
    param <- list(
      X = time,
      Y = Cobs,
      K = knot.func(par)
    )
    Chat <- predict.spline.func(param)
    err <- Cobs-Chat
    sigma <- sd(err)
    loglik <- dnorm(Cobs, Chat, sigma, log = T)
    OFVmle <- -1*sum(loglik)
    OFVmle
  }

# -----------------------------------------------------------------------------
# Set up the knot optimising function
  which.knots <- function(x, y, max.knots) {
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
        method = "L-BFGS-B", lower = min(x) + 0.0001, upper = max(x),
        Cobs = y, time = x  # Function arguments
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
# Determine times for all datasets
  times <- seq(from = 0, to = 24, by = 0.25)	# Time sequence for simulating concentrations

# Simulate dataset
  K <- c(5, 10)  # knots
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  ERR <- 1 + rnorm(n = length(times), mean = 0, sd = 0.3)
  dose <- 50	# mg

  sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)	# Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))
  plot(times, log(conc))
  sample.conc <- (conc*ERR)[times %in% sample.times]
  onecomp <- data.frame(time = sample.times, conc = sample.conc)

# Determine if there is a difference between degrees
  which.knots(onecomp$time, onecomp$conc, 3)


  spline.to.plot <- function(data, spline) {
    plot.df <- ldply(seq_len(length(spline)), function(n) {
      knots <- spline[[n]]$knots
      mle.vec <- spline[[n]]$error
      ldply(knots, function(K) {
        param <- list(
          X = data$time,
          Y = data$conc,
          K = K
        )
        K.len <- length(K)
        K.vec <- c(K, rep(NULL, length(knots) - K.len))
        mle <- signif(mle.vec[[K.len]], 3)
        out <- data.frame(
          time = seq(from = 0, to = 24, by = 0.25),
          conc = uber.predict.spline.func(param)
          nknots = length(K),
          spline.err = mle
        )
        err <- conc - out$conc
        sigma <- sd(err)
        loglik <- dnorm(conc, out$conc, sigma, log = T)
        out$error <- signif(-1*sum(loglik), 3)
        for (i in 1:length(knots)) {
          out[[paste0("knot", i)]] <- K.vec[i]
        }
        out
      })
    })  # ldply
    plot.df$knots <- unlist(
      dlply(plot.df, .(nknots), function(x) {
        K <- unique(x$nknots)
        df.knots <- NULL
        for (i in 1:K) {
          df.knots <- paste(df.knots, round(x[[paste0("knot", i)]], 1))
        }
        df.knots
      })
    )
    plot.df$facet <- factor(
      paste0("Knots: ", plot.df$nknots, "\n",
        plot.df$knots
      )
    )
    plot.df
  }  # spline.to.plot

# -----------------------------------------------------------------------------
# Naturally a 6 degree spline fits the data the best but we aren't interested in
# the sse for spline against the observed, but infact the spline against the
# truth!
  plotdata <- spline.to.plot(onecomp, spline.data)
  truedata <- data.frame(
    times = rep(times, times = length(levels(plotdata$facet))),
    conc = rep(conc, times = length(levels(plotdata$facet))),
    facet = factor(rep(levels(plotdata$facet), each = length(times)))
  )

  plotobj <- NULL
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = time, y = conc), data = plotdata, size = 1, colour = "red")
  plotobj <- plotobj + geom_line(aes(x = times, y = conc), data = truedata, size = 0.5, linetype = 2)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.6, label = error), data = plotdata)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.3, label = spline.err), data = plotdata, colour = "red")
  plotobj <- plotobj + scale_y_log10("Concentration (mg/mL) \n", lim = c(0.01, 0.8))
  plotobj <- plotobj + scale_x_continuous("\n Time (hours)", lim = c(0, 24))
  plotobj <- plotobj + facet_wrap(~facet, ncol = 4)
  plotobj

# -----------------------------------------------------------------------------
# Simulate second dataset with more rich sampling
  sample.times <- c(0,0.25,0.5,1,2,3,4,6,8,10,12,18,24)	# Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))
  sample.conc <- (conc*ERR)[times %in% sample.times]
  onecomp <- data.frame(time = sample.times, conc = sample.conc)

# Determine if there is a difference between degrees
  spline.data <- llply(seq_len(4), function(x) {
    which.knots(onecomp$time, onecomp$conc, 4, x)
  })

  plotdata <- spline.to.plot(onecomp, spline.data)
  truedata <- data.frame(
    times = rep(times, times = length(levels(plotdata$facet))),
    conc = rep(conc, times = length(levels(plotdata$facet))),
    facet = factor(rep(levels(plotdata$facet), each = length(times)))
  )

  plotobj <- NULL
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = time, y = conc), data = plotdata, size = 1, colour = "red")
  plotobj <- plotobj + geom_line(aes(x = times, y = conc), data = truedata, size = 0.5, linetype = 2)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.6, label = error), data = plotdata)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.3, label = spline.err), data = plotdata, colour = "red")
  plotobj <- plotobj + scale_y_log10("Concentration (mg/mL) \n", lim = c(0.01, 0.8))
  plotobj <- plotobj + scale_x_continuous("\n Time (hours)", lim = c(0, 24))
  plotobj <- plotobj + facet_wrap(~facet, ncol = 4)
  plotobj

# -----------------------------------------------------------------------------
#times <- seq(from = 0, to = 24, by = 0.25)	# Time sequence for simulation

# Simulate two-compartmental kinetic drug
  CL <- 10
  V2 <- 30
  Q <- 15
  V3 <- 100
  KA <- 0.5
  ERR <- 1 + rnorm(n = length(times), mean = 0, sd = 0.3)
  dose <- 50	# mg
  k10 <- CL/V2
  k12 <- Q/V2
  k21 <- Q/V3
  apb <- k10 + k12 + k21            # alpha + beta
  amb <- k10*k21                # alpha * beta
  alpha <- (apb + sqrt(apb^2 - 4*amb))/2
  beta <- (apb - sqrt(apb^2 - 4*amb))/2
  A <- KA*(k21 - alpha)/(V2*(KA - alpha)*(beta - alpha))
  B <- KA*(k21 - beta)/(V2*(KA - beta)*(alpha - beta))

  sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)	# Sampling times
  conc <- dose*(A*exp(-alpha*times) + B*exp(-beta*times) - (A+B)*exp(-KA*times))
  plot(times, log(conc))
  sample.conc <- (conc*ERR)[times %in% sample.times]
  twocomp <- data.frame(time = sample.times, conc = sample.conc)

# Determine if there is a difference between degrees
  spline.data <- llply(seq_len(4), function(x) {
    which.knots(twocomp$time, twocomp$conc, 4, x)
  })

  plotdata <- spline.to.plot(twocomp, spline.data)
  truedata <- data.frame(
    times = rep(times, times = length(levels(plotdata$facet))),
    conc = rep(conc, times = length(levels(plotdata$facet))),
    facet = factor(rep(levels(plotdata$facet), each = length(times)))
  )

  plotobj <- NULL
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = time, y = conc), data = plotdata, size = 1, colour = "red")
  plotobj <- plotobj + geom_line(aes(x = times, y = conc), data = truedata, size = 0.5, linetype = 2)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.6, label = error), data = plotdata)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.3, label = spline.err), data = plotdata, colour = "red")
  plotobj <- plotobj + scale_y_log10("Concentration (mg/mL) \n", lim = c(0.01, 0.8))
  plotobj <- plotobj + scale_x_continuous("\n Time (hours)", lim = c(0, 24))
  plotobj <- plotobj + facet_wrap(~facet, ncol = 4)
  plotobj

# -----------------------------------------------------------------------------
# Simulate second dataset with more rich sampling
  sample.times <- c(0,0.25,0.5,1,2,3,4,6,8,10,12,18,24)	# Sampling times
  conc <- dose*(A*exp(-alpha*times) + B*exp(-beta*times) - (A+B)*exp(-KA*times))
  sample.conc <- (conc*ERR)[times %in% sample.times]
  twocomp <- data.frame(time = sample.times, conc = sample.conc)

# Determine if there is a difference between degrees
  spline.data <- llply(seq_len(4), function(x) {
    which.knots(twocomp$time, twocomp$conc, 4, x)
  })

  plotdata <- spline.to.plot(twocomp, spline.data)
  truedata <- data.frame(
    times = rep(times, times = length(levels(plotdata$facet))),
    conc = rep(conc, times = length(levels(plotdata$facet))),
    facet = factor(rep(levels(plotdata$facet), each = length(times)))
  )

  plotobj <- NULL
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = time, y = conc), data = plotdata, size = 1, colour = "red")
  plotobj <- plotobj + geom_line(aes(x = times, y = conc), data = truedata, size = 0.5, linetype = 2)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.6, label = error), data = plotdata)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.3, label = spline.err), data = plotdata, colour = "red")
  plotobj <- plotobj + scale_y_log10("Concentration (mg/mL) \n", lim = c(0.01, 0.8))
  plotobj <- plotobj + scale_x_continuous("\n Time (hours)", lim = c(0, 24))
  plotobj <- plotobj + facet_wrap(~facet, ncol = 4)
  plotobj
