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

# Flexible testing of number of exponents
multi.mle.sumexp <- function(data, absorp = F, nexp = 3) {
  x <- data[which(data[, 2] != 0), 1]
  y <- data[which(data[, 2] != 0), 2]
  sub.y <- which(y == max(y)):length(y)
  opt.par <- list(NULL)
  opt.val <- list(NULL)
  opt.gra <- list(NULL)
  opt.con <- list(NULL)
  opt.mes <- list(NULL)
  lm.par <- unname(lm(log(y[sub.y]) ~ x[sub.y])$coefficients)[c(2, 1)]
  for (i in 1:nexp) {
    if (i == 1) {
      init.par <- lm.par
    } else if (i == 2) {
      if (absorp) {
        init.par <- c(lm.par[1]*c(0.8, 1.2), rep(lm.par[2], 2))
      } else {
        init.par <- c(lm.par[1]*c(0.8, 1.2), lm.par[2]*c(0.9, 1.1))
      }
    } else {
      if (absorp) {
        init.par <- c(
          mean(optres$par[1:(i-1)]), optres$par[1:(i-1)],
          mean(optres$par[i:(2*i-2)])/2, optres$par[i:(2*i-2)]
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
      lower = c(rep(-Inf, i), rep(-1/10^10, i)),
      upper = c(rep(1/10^10, i), rep(Inf, i)),
      x = x, y = y, abs = absorp
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

# Plotting for function multi.mle.sumexp
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
