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
    stat6 <-  quantile(x, probs = 0.80, na.rm = T, names = F)
    stat7 <-  quantile(x, probs = 0.85, na.rm = T, names = F)
    stat8 <-  quantile(x, probs = 0.90, na.rm = T, names = F)
    stat9 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    result <- c("05perct" = stat1, "10perct" = stat2,
      "25perct" = stat3, "50perct" = stat4, "75perct" = stat5,
      "80perct" = stat6, "85perct" = stat7, "90perct" = stat8, "95perct" = stat9)
    result
  }

# -----------------------------------------------------------------------------
# New functions
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

  # optim.interv.ga <- function(fit.par, tmax = NULL) {
  #   absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
  #   x <- times[order(times)]
  #   init.par <- x[-c(1, length(x))]
  #   if (!is.null(tmax)) init.par <- init.par[-length(init.par)]
  #
  #   gares <- ga("real-valued",
  #     mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
  #     min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
  #     max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
  #     selection = gareal_lrSelection,
  #     crossover = gareal_spCrossover,
  #     mutation = gareal_raMutation,
  #     maxiter = 50,
  #     popSize = 250,
  #     monitor = F
  #   )
  #
  #   res <- optim(
  #     init.par,
  #     err.interv,
  #     method = "L-BFGS-B", control = c(maxit = 500),
  #     lower = min(x), upper = max(x),
  #     exp.par = fit.par, tfirst = min(x) + 0.01, tlast = max(x) - 0.01, tmax = tmax,
  #     a = absorp
  #   )
  #   return(res)
  # }

  optim.interv.div <- function(fit.par, tlast, nobs, div = c(0.25, 0.5, 1, 2, 4), tmax = F) {
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    best.res <- list(NULL)
    tmax.val <- tmax.sumexp(fit.par, tlast = tlast)
    if (tmax) tmax.opt <- tmax.val
    else tmax.opt <- NULL
    for (i in 1:length(div)) {
      x <- exp(seq(log(div[i]*tmax.val), log(tlast), length.out = nobs - 1))
      init.par <- x[-length(x)]
      if (!is.null(tmax.opt)) init.par <- init.par[-length(init.par)]
      res <- optim(
        init.par,
        err.interv,
        method = "L-BFGS-B", control = c(maxit = 500),
        lower = min(x), upper = max(x),
        exp.par = fit.par, tfirst = min(x) + 0.01, tlast = max(x) - 0.01, tmax = tmax,
        a = absorp
      )
      if (res$convergence == 0) {
        if (is.null(best.res[[1]])) best.res <- res
        else if ((res$value + 0.1) < best.res$value) best.res <- res
      }
    }
    return(best.res)
  }

# Set solution seems to work, and quickly!
  optim.interv.div(p, 24, 9)
# -----------------------------------------------------------------------------
# Testing iterative solution -> random
  t1 <- seq(0, 24, length.out = 9)
  t2 <- c(0, 1, 2, 3, 4, 5, 6, 7, 24)
  p <- c(-0.8243497, -0.8205869, -1.2830615, 3.6208785, 4.5177880)

  rep.par <- c(1:5)
  sumfunc.list <- list(NULL)
  table.list <- list(NULL)
  for (j in 1:5) {
    res.val <- c(NULL)
    res.con <- c(NULL)
    for (i in 1:1000) {
      opt <- optim.interv.rep(p, rep.par[j])
      res.val[i] <- opt$value
      res.con[i] <- opt$convergence
    }
    sumfunc.list[[j]] <- sumfuncPercentile(res.val)
    table.list[[j]] <- table(res.con)
  }

  rep.par <- c(1:5)
  sumfunc.list <- list(NULL)
  table.list <- list(NULL)
  for (j in 5) {
    res.val <- c(NULL)
    res.con <- c(NULL)
    for (i in 1:1000) {
      opt <- optim.interv.rep(p, rep.par[j], tmax = 1)
      res.val[i] <- opt$value
      res.con[i] <- opt$convergence
    }
    sumfunc.list[j] <- sumfuncPercentile(res.val)
    table.list[j] <- table(res.con)
  }

# -----------------------------------------------------------------------------
