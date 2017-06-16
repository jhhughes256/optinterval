# Fixing optiminterv as per Issue #9
# -----------------------------------------------------------------------------
# Set up old functions
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

# Trapezoidal error function for interval optimisation
  err.interv <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F) {
    times <- c(tfirst, par, tlast, tmax)
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

# Interval optimising function
  optim.interv <- function(times, fit.par, tmax = NULL) {
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

  tmax.sumexp <- function(fit.par, tlast = 24, res = 0.1) {
    times <- seq(0, tlast, by = res)
    yhat <- pred.sumexp(fit.par, times)
    return(times[which(yhat == max(yhat))])
  }

# -----------------------------------------------------------------------------
# New functions
# Trapezoidal error function for interval optimisation
  err.interv.new <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F) {
    if (a) {
      browser()
      times <- c(tfirst, par, tlast, tmax)
      times <- times[order(times)]
      deltat <- diff(times)
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))][1]
      }
    } else {
      times <- c(tfirst, par, tlast)
      deltat <- diff(times)
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    sum(err^2)
  }

# Interval optimising function
  optim.interv.new <- function(times, fit.par, tmax = NULL) {
    x <- times[order(times)]
    l <- length(x)
    init.par <- x[-c(1, (l-1), l)]
    if (!is.null(tmax)) init.par <- init.par[-length(init.par)]
    xmin <- min(x)
    xmax <- max(x)
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    res <- optim(
      init.par,
      err.interv.new,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = fit.par, tfirst = xmin + 0.01, tlast = xmax - 0.01, tmax = tmax,
      a = absorp
    )
    return(res)
  }

# -----------------------------------------------------------------------------
# Reproducible example from issue #9
  t1 <- seq(0, 24, length.out = 9)
  p <- c(-0.8243497, -0.8205869, -1.2830615, 3.6208785, 4.5177880)
  optim.interv(t1, p)
  optim.interv.new(t1, p)

  optim.interv(t1, p, tmax = tmax.sumexp(p))
  optim.interv.new(t1, p, tmax = tmax.sumexp(p))
