# Sum of exponential functions for sourcing
# -----------------------------------------------------------------------------
# The functions in order are:
#   pred.sumexp - gives sum of exponentials for given set of parameters
#   mle.sumexp - maximum likelihood estimation function to be minimised
#   optim.sumexp - determine parameters for differing numbers of exponentials
#   err.interv - trapezoidal error function to be minimised
#   optim.interv - determine intervals that give the smallest error
# -----------------------------------------------------------------------------
# Sum of exponentials predicted concentrations function
  pred.sumexp <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    m <- x[1:n]  # new
    ord <- order(m, decreasing = T)  # new
    p <- c(m[ord], x[(n+1):l])  # new
    for (i in 1:n) {
      if (i == 1) y <- p[i]^d*exp(p[i]*t + p[n+i])
      else if (i != n | a == 0) y <- y + p[i]^d*exp(p[i]*t + p[n+i])
      else if (a == 1) y <- y - p[i]^d*exp(p[i]*t)*sum(exp(p[(n+1):(2*n-1)]))
    }
    return(y)
  }

# Old function
# pred.sumexp <- function(x, t, d = 0) {
#   l <- length(x)
#   a <- ifelse(l %% 2 == 0, 0, 1)
#   n <- ceiling(l/2)
#   for (i in 1:n) {
#     if (i == 1) y <- x[i]^d*exp(x[i]*t + x[n+i])
#     else if (i != n | a == 0) y <- y + x[i]^d*exp(x[i]*t + x[n+i])
#     else if (a == 1) y <- y - x[i]^d*exp(x[i]*t)*sum(exp(x[(n+1):(2*n-1)]))
#   }
#   return(y)
# }

# Maximum likelihood estimation function for parameter optimisation
  mle.sumexp <- function(par, x, y, sigma, ga = F) {
    z <- ifelse(ga, 1, -1)
    yhat <- pred.sumexp(par, x)
    loglik <- dnorm(y, yhat, abs(yhat)*sigma, log = T)
    return(z*sum(loglik))
  }

# Fit sum of exponentials to curve for different numbers of exponentials
  optim.sumexp <- function(data, oral = F, nexp = 3) {
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    lmres <- unname(lm(log(y[lm.sub]) ~ x[lm.sub])$coefficients)
    if (lmres[2] > 0) lmres[2] <- -lmres[2]  # NOT FOR REAL FUNCTION
    # if (lmres[2] > 0) stop("Data Error: No terminal phase found")
    for (i in 1:nexp) {
      if (i == 1 & !oral) {
        optres <- list(
          par = c(lmres[2], lmres[1]),
          value = mle.sumexp(unname(c(lmres[2], lmres[1])), x, y, 0.01),
          counts = NULL, convergence = 0, message = NULL
        )
      } else {
        gares <- ga("real-valued",
          mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 50,
          popSize = 1250,
          monitor = F
        )
        optres <- optim(
          gares@solution[1, ],
          mle.sumexp,
          method = "BFGS",
          x = x, y = y, sigma = 0.01
        )
      }
      slope.par <- optres$par[1:(i+oral)]
      slope.ord <- order(slope.par, decreasing = T)
      par.ord <- unname(c(
        slope.par[slope.ord],
        optres$par[(i+oral+1):length(optres$par)]
      ))
      opt.par[[i]] <- par.ord
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

# Chi-squared difference test
# Takes a list of optim results and gives the best optim result
  chisq.sumexp <- function(opt) {
    i <- 1
    for (j in 2:length(opt$par)) {
      degf <- length(opt$par[[j]]) - length(opt$par[[i]])
      x <- opt$value[[i]] - opt$value[[j]]
      p <- pchisq(x, degf, lower.tail = F)
      if (p < 0.05) {
        i <- i + 1
      }
    }
    return(sapply(opt, function(x) x[i]))
  }

# Trapezoidal error function for interval optimisation
  err.interv <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F) {
    times <- c(tfirst, par, tlast, tmax)
    if (!is.null(tmax)) times <- times[order(times)]
    deltat <- diff(times)
    if (a) {
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
      }
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    sum(err^2)
  }

# Old function
# err.interv <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F) {
#   times <- c(tfirst, par, tlast, tmax)
#   deltat <- diff(times)
#   if (a) {
#     all.secd <- abs(pred.sumexp(exp.par, times, 2))
#     secd <- c(NULL)
#     for (i in 1:(length(times)-1)) {
#       secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
#     }
#   } else {
#     secd <- pred.sumexp(exp.par, times[-length(times)], 2)
#   }
#   err <- deltat^3*secd/12
#   sum(err^2)
# }

# Interval optimising function
  optim.interv.old <- function(times, fit.par, tmax = NULL) {
    x <- times[order(times)]
    init.par <- x[-c(1, length(x))]
    if (!is.null(tmax)) init.par <- init.par[-length(init.par)]
    xmin <- min(x)
    xmax <- max(x)
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    res <- optim(
      init.par,
      err.interv,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = fit.par, tfirst = xmin + 0.01, tlast = xmax - 0.01, tmax = tmax,
      a = absorp
    )
    return(res)
  }

  optim.interv.rep <- function(fit.par, jmax, tmax = NULL) {
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    best.res <- list(NULL)
    j <- 0
    repeat {
      x <- c(0, exp(runif(7, 0, log(24))), 24)
      x <- x[order(x)]
      init.par <- x[-c(1, length(x))]
      if (!is.null(tmax)) init.par <- init.par[-length(init.par)]
      res <- optim(
        init.par,
        err.interv,
        method = "L-BFGS-B", control = c(maxit = 500),
        lower = min(x), upper = max(x),
        exp.par = fit.par, tfirst = min(x) + 0.01, tlast = max(x) - 0.01,
        tmax = tmax, a = absorp
      )
      if (res$convergence == 0) {
        if (j == 0) {
          best.res <- res
          j <- 1
        } else if ((res$value + 0.1) < best.res$value) {
          best.res <- res
          j <- 1
        } else {
          j <- j + 1
        }
      }
      if (j == jmax) break
    }
    return(best.res)
  }

# -----------------------------------------------------------------------------
# Determine tmax given a set of sumexp parameters
  tmax.sumexp <- function(fit.par, tlast = 24, res = 0.1) {
    times <- seq(0, tlast, by = res)
    yhat <- pred.sumexp(fit.par, times)
    return(times[which(yhat == max(yhat))])
  }

# Determine auc given a set of intervals
  auc.interv <- function(times, fit.par, fn, log = F) {
    C <- do.call(fn, list(x = times, p = fit.par))
    auc <- c(NULL)
    for (i in 2:length(C)) {
      h <- times[i] - times[i-1]
      dC <- C[i-1] - C[i]
      if (log & dC > 0) auc[i-1] <- dC*h/log(C[i-1]/C[i])
      else auc[i-1] <- (C[i-1] + C[i])*h/2
    }
    return(sum(auc))
  }

# Without for loop
#  auc.interv <- function(times, fit.par, fn, log = F) {
#    C <- do.call(fn, list(x = times, p = fit.par))
#    h <- diff(times)
#    EC <- psum(C[-1], C[-length(C)])
#    dC <- diff(-C)
#    if (!log) auc <- EC*h/2
#    else auc[i-1] <- dC*h/log(C[i-1]/C[i])
#    return(sum(auc))
#  }

# Plot random data
  plot.rdata <- function(data, t, n, interv) {
    plotdata <- data.frame(
      id = rep(1:n, each = length(t)),
      time = rep(t, times = n),
      dv = as.vector(data)
    )
    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
    plotobj <- plotobj + geom_line(aes(x = time, y = dv), colour = "red")
    plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)")
    plotobj <- plotobj + facet_wrap(~id, ncol = round(sqrt(n)), scales = "free")
    return(plotobj)
  }
