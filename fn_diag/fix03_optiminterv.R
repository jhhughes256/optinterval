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

  sumfuncPercentile <- function(x) {
    stat1 <-  quantile(x, probs = 0.05, na.rm = T, names = F)
    stat2 <-  quantile(x, probs = 0.10, na.rm = T, names = F)
    stat3 <-  quantile(x, probs = 0.25, na.rm = T, names = F)
    stat4 <-  quantile(x, probs = 0.50, na.rm = T, names = F)
    stat5 <-  quantile(x, probs = 0.75, na.rm = T, names = F)
    stat6 <-  quantile(x, probs = 0.90, na.rm = T, names = F)
    stat7 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    result <- c("05perct" = stat1, "10perct" = stat2,
      "25perct" = stat3, "50perct" = stat4, "75perct" = stat5,
      "90perct" = stat6, "95perct" = stat7)
    result
  }

# -----------------------------------------------------------------------------
# New functions
# Trapezoidal error function for interval optimisation
  err.interv.new <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F) {
    times <- c(tfirst, par, tlast, tmax)
    if (!is.null(tmax)) times <- times[order(times)]
    # browser()
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
  optim.interv.alt <- function(times, fit.par, tmax = NULL) {
    x <- times[order(times)]
    init.par <- x[-c(1, length(x))]
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

  optim.interv.new <- function(fit.par, jmax, tmax = NULL) {
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
        err.interv.new,
        method = "L-BFGS-B", control = c(maxit = 500),
        lower = min(x), upper = max(x),
        exp.par = fit.par, tfirst = min(x) + 0.01, tlast = max(x) - 0.01, tmax = tmax,
        a = absorp
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
    return(res)
  }

# -----------------------------------------------------------------------------
# Reproducible example from issue #9, no tmax
  t1 <- seq(0, 24, length.out = 9)
  t2 <- c(0, 1, 2, 3, 4, 5, 6, 7, 24)
  p <- c(-0.8243497, -0.8205869, -1.2830615, 3.6208785, 4.5177880)
  # optim.interv(t1, p)
  # optim.interv(t2, p)
  # optim.interv.new(t1, p)

  res.val.lin <- c(NULL)
  res.con.lin <- c(NULL)
  res.val.log <- c(NULL)
  res.con.log <- c(NULL)
  for (i in 1:1000) {
    t3 <- c(0, exp(runif(7, 0, log(24))), 24)
    t3 <- t3[order(t3)]
    opt <- optim.interv.alt(t3, p)
    res.val.log[i] <- opt$value
    res.con.log[i] <- opt$convergence
    t3 <- c(0, runif(7, 0, 24), 24)
    t3 <- t3[order(t3)]
    opt <- optim.interv.alt(t3, p)
    res.val.lin[i] <- opt$value
    res.con.lin[i] <- opt$convergence
  }
  sumfuncPercentile(res.val.lin)
  sumfuncPercentile(res.val.log)
  table(res.con.lin)
  table(res.con.log)

  res.val.log <- c(NULL)
  res.con.log <- c(NULL)
  t.m4 <- c(0, exp(seq(log(4), log(24), length.out = 8)))
  t.m2 <- c(0, exp(seq(log(2), log(24), length.out = 8)))
  t.m1 <- c(0, exp(seq(log(1), log(24), length.out = 8)))
  t.d2 <- c(0, exp(seq(log(1/2), log(24), length.out = 8)))
  t.d4 <- c(0, exp(seq(log(1/4), log(24), length.out = 8)))
  t.d8 <- c(0, exp(seq(log(1/8), log(24), length.out = 8)))
  t.d16 <- c(0, exp(seq(log(1/24), log(24), length.out = 8)))
  optim.interv.alt(t.m4, p)
  optim.interv.alt(t.m2, p)
  optim.interv.alt(t.m1, p)
  optim.interv.alt(t.d2, p)
  optim.interv.alt(t.d4, p)
  optim.interv.alt(t.d8, p)
  optim.interv.alt(t.d16, p)
  for (i in 1:1000) {
    t3 <- c(0, exp(runif(7, 0, log(24))), 24)
    t3 <- t3[order(t3)]
    opt <- optim.interv.alt(t3, p)
    res.val.log[i] <- opt$value
    res.con.log[i] <- opt$convergence
  }
  sumfuncPercentile(res.val.nran)
  sumfuncPercentile(res.val.log)
  table(res.con.nran)
  table(res.con.log)

# -----------------------------------------------------------------------------
# Testing tmax
  # optim.interv(t1, p, tmax = tmax.sumexp(p))
  # optim.interv.alt(t1, p, tmax = tmax.sumexp(p))

  res.val.lin <- c(NULL)
  res.con.lin <- c(NULL)
  res.val.log <- c(NULL)
  res.con.log <- c(NULL)
  for (i in 1:1000) {
    t3 <- c(0, exp(runif(7, 0, log(24))), 24)
    t3 <- t3[order(t3)]
    opt <- optim.interv.alt(t3, p, tmax = 1)
    res.val.log[i] <- opt$value
    res.con.log[i] <- opt$convergence
    t3 <- c(0, runif(7, 0, 24), 24)
    t3 <- t3[order(t3)]
    opt <- optim.interv.alt(t3, p, tmax = 1)
    res.val.lin[i] <- opt$value
    res.con.lin[i] <- opt$convergence
  }
  sumfuncPercentile(res.val.lin)
  sumfuncPercentile(res.val.log)
  table(res.con.lin)
  table(res.con.log)

# -----------------------------------------------------------------------------
# All this testing is useless, function returns the last result, not the best!
# See fix04... for proper analysis
  optim.interv.new(p, 15)
  res.val <- c(NULL)
  res.con <- c(NULL)
  for (i in 1:1000) {
    opt <- optim.interv.new(p, 3)
    res.val[i] <- opt$value
    res.con[i] <- opt$convergence
  }
  sumfuncPercentile(res.val)
  table(res.con)
  # jmax values
  # - 5
  # 05perct   10perct   25perct   50perct   75perct   90perct   95perct
  # 6.775586  6.775586  6.775586  6.775586  6.775586 38.315197 41.629376
  # - 8
  # 05perct   10perct   25perct   50perct   75perct   90perct   95perct
  # 6.775586  6.775586  6.775586  6.775586  6.775586 38.315139 41.217015
  # - 10
  # 05perct   10perct   25perct   50perct   75perct   90perct   95perct
  # 6.775586  6.775586  6.775586  6.775586  6.775586 38.320782 41.970952
  # - 15  # note: only 722 loops
  # 05perct   10perct   25perct   50perct   75perct   90perct   95perct
  # 2.931696  6.775586  6.775586  6.775586  6.775586 38.353301 41.510440

  res.val <- c(NULL)
  res.con <- c(NULL)
  for (i in 1:1000) {
    opt <- optim.interv.new(p, 3, tmax = 1)
    res.val[i] <- opt$value
    res.con[i] <- opt$convergence
  }
  sumfuncPercentile(res.val)
  table(res.con)
  # jmax values
  # - 5
  # 05perct    10perct    25perct    50perct    75perct    90perct    95perct
  # 3.388907  98.833650 100.839873 100.839944 100.840255 100.842071 101.479957
  # - 8
  # 05perct    10perct    25perct    50perct    75perct    90perct    95perct
  # 3.315665  98.830237 100.839850 100.839933 100.840250 100.843184 101.480487
  # - 10
  # 05perct    10perct    25perct    50perct    75perct    90perct    95perct
  # 2.164885  98.830761 100.839855 100.839935 100.840272 100.844184 101.481981
  # - 15  # note: only 838 loops
  # 05perct    10perct    25perct    50perct    75perct    90perct    95perct
  # 1.462877   6.795721 100.839844 100.839927 100.840197 100.842331 101.479591
