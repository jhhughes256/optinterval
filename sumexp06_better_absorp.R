# Script with improved method of optimising for oral data
# -----------------------------------------------------------------------------
# It was possible that the sum of all intercepts is not equal to zero
# For absorption curves this is wrong, as they should always pass through zero
#
# An issue that needs to be fixed to allow this (which should also be fixed as
# it is wrong) the boundaries on intercepts used don't apply to log-transformed
# intercepts
# This means that intercepts should go into optim without a transform applied
# and then be transformed before sum of exponentials occurs
# -----------------------------------------------------------------------------
# Setup workspace
  library(plyr)
  library(ggplot2)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# -----------------------------------------------------------------------------
# Setup OFV functions
# Flexible sum of coefficients function (updated!)
# New additions from "sumexp05" added + changed for new parameter layout
  pred.sumexp <- function(x, t, d = 0) {
  # Check length of x to see how many exponentials will be used
  # if l is odd, then an absorption curve is being made
    l <- length(x)/2
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l)
    for (i in 1:n) {
      if (i == 1) y <- x[i]^d*exp(x[i]*t + log(x[n+i]))
      else if (i != n | a == 0) y <- y + x[i]^d*exp(x[i]*t + log(x[n+i]))
      else if (a == 1) y <- y - x[i]^d*exp(x[i]*t + log(sum(x[(n+1):(2*n-1)])))
    }
    return(y)
  }

# Flexible maximum likelihood estimmation
  mle.sumexp <- function(par, x, y) {
    yhat <- pred.sumexp(par, x)
    err <- y - yhat
    sigma <- sd(err)
    loglik <- dnorm(y, yhat, sigma, log = T)
    return(-1*sum(loglik))
  }

  multi.mle.sumexp <- function(data, oral = F, nexp = 3) {
    x <- data[which(data[, 2] != 0), 1]
    y <- data[which(data[, 2] != 0), 2]
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    sub.y <- which(y == max(y)):length(y)
    lm.par <- unname(lm(log(y[sub.y]) ~ x[sub.y])$coefficients)
    for (i in 1:nexp) {
      if (i == 1) {
        init.par <- c(lm.par[2], exp(lm.par[1]))
      } else if (i == 2) {
        browser()
        if (oral) {
          init.par <- c(lm.par[2]*c(0.8, 1.2), exp(lm.par[1]*c(1, 1.2))[-i])
        } else {
          init.par <- c(lm.par[2]*c(0.8, 1.2), exp(lm.par[1]*c(1, 1.2)))
        }
      } else {
        seq(1-(i-(oral+1))*0.05, 1+(i-(oral+1))*0.05, length.out = i-oral-1)
        if (oral) {
          init.par <- c(
            mean(optres$par[1:(i-1)]), optres$par[1:(i-1)],
            sum(optres$par[i:(2*i-3)])/(i-2)
          )
        } else {
          init.par <- c(
            mean(optres$par[1:(i-1)]), optres$par[1:(i-1)],
            mean(optres$par[i:(2*i-2)]), optres$par[i:(2*i-2)]
          )
        }
      }
      optres <- optim(
        init.par,
        mle.sumexp,  # maximum likelihood fitting function
        method = "L-BFGS-B",
        lower = c(rep(-Inf, i), rep(1/10^10, i)),
        upper = c(rep(-1/10^10, i), rep(Inf, i)),
        x = x, y = y
      )
      opt.par[[i]] <- optres$par
      opt.val[[i]] <- optres$value
      opt.gra[[i]] <- optres$counts
      opt.con[[i]] <- optres$convergence
      opt.mes[[i]] <- ifelse(is.null(optres$message),
        "NULL", optres$message)
    }
    res <- list(par = opt.par, value = opt.val, counts = opt.gra,
      convergence = opt.con, message = opt.mes)
    res
  }

  plot.sumexp <- function(res, data) {
    plotdata <- ldply(res$par, function(x) {
      data.frame(
        nexp = as.factor(length(x)/2),
        time = data$time,
        cobs = data$conc,
        pred = pred.sumexp(x, data$time, 0),
        objv = res$value[[length(x)/2]]
      )
    })
    levels(plotdata$nexp) <- paste(levels(plotdata$nexp), "exponent(s)")
    ylim <- c(0, 1.1*max(plotdata$cobs))
    xlim <- c(0, max(plotdata$time))

    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Comparison of Number of Exponents Used")
    plotobj <- plotobj + geom_point(aes(x = time, y = cobs))
    plotobj <- plotobj + geom_line(aes(x = time, y = pred), colour = "red")
    plotobj <- plotobj + geom_text(aes(x = median(time), y = max(cobs), label = objv))
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n", lim = ylim)
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
    plotobj <- plotobj + facet_wrap(~nexp)
    plotobj
  }
# -----------------------------------------------------------------------------
# Simulate data
# Absorption Curve
  time.samp <- seq(0, 48, by = 2)
  absdata <- data.frame(
    time = time.samp,
    line1 = -0.2*time.samp + 4,
    line2 = -0.1*time.samp + 4
  )
  absdata$sumexp <- exp(absdata$line2) - exp(absdata$line1)
  #with(absdata, plot(time, sumexp))

# 2 Compartment Curve
  twodata <- data.frame(
    time = time.samp,
    line1 = -0.5*time.samp + 6,
    line2 = -0.05*time.samp + 5
  )
  twodata$sumexp <- exp(twodata$line1) + exp(twodata$line2)
  #with(twodata, plot(time, log(sumexp)))
# -----------------------------------------------------------------------------
# Add random error and optimise using rich data
# Absorption
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data1 <- data.frame(
    time = time.samp,
    conc = absdata$sumexp*err
  )
  with(data1, plot(time, log(conc)))

  all.res <- multi.mle.sumexp(data1, nexp = 4, oral = T)

  plot.sumexp(all.res, data1, oral = T)

# Two compartment
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data2 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp*err
  )
  with(data2, plot(time, log(conc)))

  all.res <- multi.mle.sumexp(data2, nexp = 4)
  #with(all.res, which(unlist(value) == min(unlist(value))))  # best fit
  plot.sumexp(all.res, data2)
