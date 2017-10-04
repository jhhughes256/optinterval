# Fixing optiminterv as per Issue #9
# -----------------------------------------------------------------------------
# Load libraries
  library(GA)
  library(ggplot2)

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
  err.interv.ga <- function(par, exp.par, tfirst, tlast, tmax = NULL, a = F, ga = F) {
    z <- ifelse(ga, -1, 1)
    times <- c(tfirst, par, tlast, tmax)
    times <- times[order(times)]
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
    return(z*sum(err^2))
  }

# Interval optimising function
  optim.interv.ga <- function(fit.par, times, tmax = NULL) {
    tfirst <- min(times)
    tlast <- max(times)
    npar <- length(times)
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    gares <- ga("real-valued",
      err.interv.ga, exp.par = fit.par, a = absorp,
      tfirst = tfirst, tlast = tlast, tmax = tmax, ga = T,
      min = rep(tfirst + 0.01, npar), max = rep(tlast - 0.01, npar),
      selection = gareal_lrSelection,
      crossover = gareal_spCrossover,
      mutation = gareal_raMutation,
      maxiter = 50,
      popSize = 250,
      monitor = F
    )
    res <- optim(
      gares@solution[order(gares@solution)],
      err.interv.ga,
      method = "BFGS", hessian = T,
      exp.par = fit.par, tfirst = tfirst, tlast = tlast, tmax = tmax,
      a = absorp
    )
    return(res)
  }

# -----------------------------------------------------------------------------
# Reproducible example from issue #9, no tmax
  nsamp <- 9
  xmin <- 0
  xmax <- 24
  p <- c(-0.8243497, -0.8205869, -1.2830615, 3.6208785, 4.5177880)
  t1 <- c(rep(xmin, nsamp - 1), xmax)
  res.par <- list(NULL)
  res.val <- c(NULL)
  res.cou <- list(NULL)
  res.con <- c(NULL)
  res.mes <- list(NULL)
  res.hes <- list(NULL)
  for (i in 1:100) {
    opt <- optim.interv.ga(p, t1)
    res.par[[i]] <- opt$par
    res.val[i] <- opt$value
    res.cou[[i]] <- opt$counts
    res.con[i] <- opt$convergence
    res.mes[[i]] <- ifelse(is.null(opt$message), NA, opt$message)
    res.hes[[i]] <- opt$hessian
  }
  sumfuncPercentile(res.val)
  table(res.con)
  fail <- which(res.con == 1)
  res.cou[fail]
  res.val[fail]
  res.par[fail]
  res.hes[fail]
  mapply(res.hes[fail], res.par[fail], FUN = function(x, y) {
    round(sqrt(diag(solve(x)))/abs(y)*100, 2)
  })
# -----------------------------------------------------------------------------
# Below is for meeting with RU, not applicable to every run of ga
  times <- seq(0, 24, by = 0.1)
  dv <- pred.sumexp(p, times)
  data <- data.frame(time = times, dv = dv)
  p1 <- NULL
  p1 <- ggplot(aes(x = time, y = dv), data = data)
  p1 <- p1 + geom_point(shape = 1)
  # 1 is good
  res.par[[1]]
  round(sqrt(diag(solve(res.hes[[1]])))/abs(res.par[[1]])*100, 2)
  p1 + geom_vline(xintercept = res.par[[1]], colour = "green4", linetype = "dashed")
  # 4, 10 for bad
  res.par[[4]]
  round(sqrt(diag(solve(res.hes[[4]])))/abs(res.par[[4]])*100, 2)
  p1 + geom_vline(xintercept = res.par[[4]], colour = "green4", linetype = "dashed")
  res.par[[10]]
  round(sqrt(diag(solve(res.hes[[10]])))/abs(res.par[[10]])*100, 2)
  p1 + geom_vline(xintercept = res.par[[10]], colour = "green4", linetype = "dashed")
  # 17 for NaN
  res.par[[17]]
  round(sqrt(diag(solve(res.hes[[17]])))/abs(res.par[[17]])*100, 2)
  diag(solve(res.hes[[17]]))
  p1 + geom_vline(xintercept = res.par[[17]], colour = "green4", linetype = "dashed")
  # 45, 46, 70 for singular hessian
  solve(res.hes[[45]])
  solve(res.hes[[46]])
  solve(res.hes[[70]])
  solve(res.hes[[71]])
  solve(res.hes[[87]])
