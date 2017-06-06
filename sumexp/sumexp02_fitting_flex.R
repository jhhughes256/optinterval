# Working out how to fit to PK data using sum of exponentials
# -----------------------------------------------------------------------------
# Email from Richard
# included excel spreadsheet using formulas to fit coefficients of linear
# equations

# Absorption Curve
  time.samp <- seq(0, 48, by = 2)
  absdata <- data.frame(
    time = time.samp,
    line1 = -0.2*time.samp + 4,
    line2 = -0.1*time.samp + 4
  )
  absdata$sumexp <- exp(absdata$line2) - exp(absdata$line1)
  with(absdata, plot(time, sumexp))

# 2 Compartment Curve
  twodata <- data.frame(
    time = time.samp,
    line1 = -0.5*time.samp + 6,
    line2 = -0.05*time.samp + 5
  )
  twodata$sumexp <- exp(twodata$line1) + exp(twodata$line2)
  with(twodata, plot(time, log(sumexp)))
# -----------------------------------------------------------------------------
# Setup OFV functions
# Flexible least squares function
  ls.sumexp <- function(par, x, y, abs = F) {
    n <- length(par)/2
    for (i in 1:n) {
      if (!exists("yhat")) yhat <- exp(par[i]*x + par[n+i])
      else if (i != n | abs == F) yhat <- yhat + exp(par[i]*x + par[n+i])
      else if (abs == T) yhat <- yhat - exp(par[i]*x + par[n+i])
    }
    err <- y - yhat
    return(sum(err^2))
  }

  # Flexible maximum likelihood estimmation
  mle.sumexp <- function(par, x, y, abs = F) {
    n <- length(par)/2
    for (i in 1:n) {
      if (!exists("yhat")) yhat <- exp(par[i]*x + par[n+i])
      else if (i != n | abs == F) yhat <- yhat + exp(par[i]*x + par[n+i])
      else if (abs == T) yhat <- yhat - exp(par[i]*x + par[n+i])
    }
    err <- y - yhat
    sigma <- sd(err)
    loglik <- dnorm(y, yhat, sigma, log = T)
    return(-1*sum(loglik))
  }

# -----------------------------------------------------------------------------
# Will require fitting to real data, parameters are slope and intercept
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data1 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp * err
  )
  with(data1, plot(time, log(conc)))

  # will need to somehow guess initial parameters, but first get optim working
  init.par <- c(-0.5, -0.05, 6, 5)

  ls.result <- optim(
    init.par,  # Initial parameter estimates
    ls.sumexp,  # Fitting function
    method = "BFGS",
    x = data1$time, y = data1$conc # Function arguments
  )
  mle.result <- optim(
    init.par,  # Initial parameter estimates
    mle.sumexp,  # Fitting function
    method = "BFGS",
    x = data1$time, y = data1$conc # Function arguments
  )

  result <- mle.result
  data1$line1 <- data1$time*result$par[1] + result$par[3]
  data1$line2 <- data1$time*result$par[2] + result$par[4]
  data1$pred <- with(data1, exp(line1) + exp(line2))
  with(data1, plot(time, log(conc)))
  with(data1, lines(time, line1, lwd = 2, col = "blue"))
  with(data1, lines(time, line2, lwd = 2, col = "blue"))
  with(data1, lines(time, log(pred), lwd = 2, col = "red"))

  with(data1, plot(time, conc))
  with(data1, lines(time, exp(line1), lwd = 2, col = "blue"))
  with(data1, lines(time, exp(line2), lwd = 2, col = "blue"))
  with(data1, lines(time, pred, lwd = 2, col = "red"))

# -----------------------------------------------------------------------------
# Will require fitting to real data, parameters are slope and intercept
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data2 <- data.frame(
    time = time.samp,
    conc = absdata$sumexp*err
  )
  with(data2, plot(time, log(conc)))

  init.par <- c(-0.1, -0.2, 4, 4)

  ls.result <- optim(
    init.par,  # Initial parameter estimates
    ls.sumexp,  # Fitting function
    method = "BFGS",
    x = data2$time, y = data2$conc, abs = T
  )
  mle.result <- optim(
    init.par,  # Initial parameter estimates
    mle.sumexp,  # Fitting function
    method = "BFGS",
    x = data2$time, y = data2$conc, abs = T
  )

  result <- mle.result
  data2$line1 <- data2$time*result$par[1] + result$par[3]
  data2$line2 <- data2$time*result$par[2] + result$par[4]
  data2$pred <- with(data2, exp(line1) + -1*exp(line2))
  with(data2, plot(time, log(conc)))
  with(data2, lines(time, line1, lwd = 2, col = "blue"))
  with(data2, lines(time, line2, lwd = 2, col = "blue"))
  with(data2, lines(time, log(pred), lwd = 2, col = "red"))

  with(data2, plot(time, conc))
  with(data2, lines(time, exp(line1), lwd = 2, col = "blue"))
  with(data2, lines(time, -exp(line2), lwd = 2, col = "blue"))
  with(data2, lines(time, pred, lwd = 2, col = "red"))
