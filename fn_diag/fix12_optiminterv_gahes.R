# Trialling using hessian matrices to determine if initial parameters should
# be resampled. Interested in the number of resamples required on average.
# -----------------------------------------------------------------------------
# Load libraries
  library(GA)
  library(plyr)
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
    npar <- length(times) - 2
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    flag <- 1
    repeat {
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
      vc_mat <- try(solve(res$hessian))
      if(class(vc_mat) != "try-error") {
        se <- sqrt(diag(vc_mat))
        if (!any(is.nan(se))) {
          res$flag <- flag
          break
        }
      }
      flag <- flag + 1
    }
    return(res)
  }

# -----------------------------------------------------------------------------
# Set up example (from #9)
  nsamp <- 9
  xmin <- 0
  xmax <- 24
  p <- c(-0.8243497, -0.8205869, -1.2830615, 3.6208785, 4.5177880)
  t1 <- c(rep(xmin, nsamp - 1), xmax)

# Run for err.interv using time
  set.seed(1989)
  res.par <- list(NULL)
  res.val <- c(NULL)
  res.cou <- list(NULL)
  res.con <- c(NULL)
  res.mes <- list(NULL)
  res.hes <- list(NULL)
  res.flag <- c(NULL)
  for (i in 1:1000) {
    opt <- optim.interv.ga(p, t1)
    res.par[[i]] <- opt$par
    res.val[i] <- opt$value
    res.cou[[i]] <- opt$counts
    res.con[i] <- opt$convergence
    res.mes[[i]] <- ifelse(is.null(opt$message), NA, opt$message)
    res.hes[[i]] <- opt$hessian
    res.flag[i] <- opt$flag
  }
  sumfuncPercentile(res.val)
  table(res.con)

  table(res.se)
  table(res.flag)
  table(res.se, res.flag)
  table(res.se, floor(res.val))

  time.samp <- seq(0, 24, by = 0.1)
  dv <- pred.sumexp(p, time.samp)
  data <- data.frame(time = time.samp, dv = dv)
  p1 <- NULL
  p1 <- ggplot(aes(x = time, y = dv), data = data)
  p1 <- p1 + geom_point(shape = 1)
  p1 + geom_vline(xintercept = res.par[[465]], colour = "green4", linetype = "dashed")
