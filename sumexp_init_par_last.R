# Further work on determining initial parameters
# -----------------------------------------------------------------------------
# Setup workspace
  library(plyr)
  library(ggplot2)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# -----------------------------------------------------------------------------
# Setup OFV functions
# Flexible sum of coefficients function
  pred.sumexp <- function(x, t, subtract.last) {
    n <- length(x)/2
    for (i in 1:n) {
      if (!exists("yhat")) yhat <- exp(x[i]*t + x[n+i])
      else if (i != n | subtract.last == F) yhat <- yhat + exp(x[i]*t + x[n+i])
      else if (subtract.last == T) yhat <- yhat - exp(x[i]*t + x[n+i])
    }
    return(yhat)
  }

  # Flexible maximum likelihood estimmation
  mle.sumexp <- function(par, x, y, abs = F) {
    yhat <- pred.sumexp(par, x, abs)
    err <- y - yhat
    sigma <- sd(err)
    loglik <- dnorm(y, yhat, sigma, log = T)
    return(-1*sum(loglik))
  }

  multi.mle.sumexp <- function(data, absorp = F, nexp = 3) {
    x <- data[which(data[, 2] != 0), 1]
    y <- data[which(data[, 2] != 0), 2]
    sub.y <- which(y == max(y)):length(y)
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    lm.par <- unname(lm(log(y[sub.y]) ~ x[sub.y])$coefficients)[c(2, 1)]
    for (i in 1:nexp) {
      if (i == 1) {
        init.par <- lm.par
      } else if (i == 2) {
        if (absorp) {
          init.par <- rep(lm.par, each = 2)*c(0.8, 1.2)
        } else {
          init.par <- rep(lm.par, each = 2)*c(0.9, 1.1)
        }
      } else {
        init.par <- c(
          mean(optres$par[1:(i-1)]), optres$par[1:(i-1)],
          mean(optres$par[i:(2*i-2)]), optres$par[i:(2*i-2)]
        )
      }
      optres <- optim(
        init.par,
        mle.sumexp,  # maximum likelihood fitting function
        method = "BFGS",
        x = x, y = y, abs = absorp
      )
      opt.par[[i]] <- optres$par
      opt.val[[i]] <- optres$value
      opt.gra[[i]] <- optres$counts
      opt.con[[i]] <- optres$convergence
    }
    res <- list(par = opt.par, value = opt.val, counts = opt.gra,
      convergence = opt.con)
    res
  }

  plot.sumexp <- function(res, data, absorp = F) {
    plotdata <- ldply(res$par, function(x) {
      data.frame(
        nexp = as.factor(length(x)/2),
        time = data$time,
        cobs = data$conc,
        pred = pred.sumexp(x, data$time, absorp),
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

  all.res <- multi.mle.sumexp(data1, nexp = 4, absorp = T)

  plot.sumexp(all.res, data1, absorp = T)

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

# -----------------------------------------------------------------------------
# Add random error and optimise using sparse
  sub.time <- c(1:5, (1:4+3)*2-1, (1:3+3)*4-1)
# Absorption sub-data
  err <- 1 + rnorm(n = length(sub.time), mean = 0, sd = 0.3)
  data3 <- data.frame(
    time = time.samp[sub.time],
    conc = absdata$sumexp[sub.time]*err
  )
  with(data3, plot(time, log(conc)))

  all.res <- multi.mle.sumexp(data3, nexp = 4, absorp = T)

  plot.sumexp(all.res, data3, absorp = T)

# Two compartment sub-data
  err <- 1 + rnorm(n = length(sub.time), mean = 0, sd = 0.3)
  data4 <- data.frame(
    time = time.samp[sub.time],
    conc = twodata$sumexp[sub.time]*err
  )
  with(data4, plot(time, log(conc)))

  all.res <- multi.mle.sumexp(data4, nexp = 4)

  plot.sumexp(all.res, data4)
