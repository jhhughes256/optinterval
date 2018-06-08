# Determine how robust err.interv.dt is and trial different methods of
# init.par generation
# -----------------------------------------------------------------------------
# Set this one up properly
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32", "C:/Users/hugjh001/Documents/optinterval")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "optinterval"
    } else if (getwd() == wd[2] | getwd() == wd[5]) {
      git.dir <- "C:/Users/hugjh001/Documents"
      reponame <- "optinterval"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "splines"
    }
    rm("wd")
  }
# Load libraries
  library(GA)
  library(plyr)
  library(ggplot2)

# Source scripts to set up environment
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "study_rdata.R", sep = "/"))

# Load New functions
# Chi-squared difference test
  chisq.interv <- function(fit.par, times, tmax = F) {
    nobs <- length(times)
    tlast <- times[nobs]
    if (tmax) tmax.val <- tmax.sumexp(fit.par, tlast) else tmax.val <- NULL
    ref <- optim.interv.dt(fit.par, times, tmax = tmax.val)
    test <- optim.interv.dt(fit.par, times[-1], tmax = tmax.val)
  }
  
  optim.interv.dt <- function(fit.par, times, tmax = NULL) {
    tfirst <- min(times)
    tlast <- max(times)
    npar <- length(times) - 2
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    init.par <- cumsum(rep(tlast/48, npar))
    res <- optim(
      init.par,
      err.interv.dt,
      method = "L-BFGS-B", hessian = T,
      lower = tlast/48, upper = tlast - npar*tlast/48,
      exp.par = fit.par, tfirst = tfirst, tlast = tlast, a = absorp, tmax = tmax
    )
    res$times <- sort(c(cumsum(res$par), tmax))
    return(res)
  }
  
  optim.interv.dtmax <- function(fit.par, times, tmax = NULL) {
    tfirst <- min(times)
    tlast <- max(times)
    npar <- length(times) - 2
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    init.par <- cumsum(rep(tlast/48, npar))
    res <- optim(
      init.par,
      err.interv.dtmax,
      method = "L-BFGS-B", hessian = T,
      lower = tlast/48, upper = tlast - npar*tlast/48,
      exp.par = fit.par, tfirst = tfirst, tlast = tlast, a = absorp, tmax = tmax
    )
    res$times <- sort(c(cumsum(res$par), tmax))
    return(res)
  }
  
  err.interv.dt <- function(par, exp.par, tfirst, tlast, a = F, tmax = NULL) {
    times <- c(tfirst, cumsum(par), tmax, tlast)
    times <- sort(times)
    deltat <- diff(times)
    if (a) {
      secd.all <- abs(pred.sumexp(exp.par, times, 2))
      secd.mat <- matrix(c(secd.all[-1], secd.all[-length(secd.all)]), ncol = 2)
      secd <- apply(secd.mat, 1, max)
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    return(sum(err^2))
  }
  
  err.interv.dtmax <- function(par, exp.par, tfirst, tlast, a = F, tmax = NULL) {
    times <- c(tfirst, cumsum(par), tmax, tlast)
    times <- sort(times)
    deltat <- diff(times)
    if (a) {
      times.mat <- matrix(c(times[-length(times)], times[-1]), ncol = 2)
      times.secd <- apply(times.mat, 1, tmax.sumexp, fit.par = exp.par, d = 2)
      secd <- abs(pred.sumexp(exp.par, times.secd, 2))
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    return(sum(err^2))
  }

# -----------------------------------------------------------------------------
# Data for testing is from study_rdata, therefore the datasets in order are:
#   d1b - one compartment kinetics given iv
#   d2b - two compartment kinetics given iv
#   d3b - three compartment kinetics given iv
#   d1a - one compartment kinetics given orally
#   d2a - two compartment kinetics given orally
#   d3a - three compartment kinetics given orally
# Additionally they are in the form of a matrix!
  t.ph <- c(rep(0, 8), 24)
  facetted.interval.plot <- function(x) {
    plot.dv <- get(paste0(x, ".df"))
    plot.int <- get(paste0(x, ".int.df"))
    p <- NULL
    p <- ggplot()
    p <- p + geom_point(aes(x = TIME, y = DV), data = plot.dv, shape = 1)
    p <- p + geom_vline(aes(xintercept = INT), data = plot.int, colour = "green4", linetype = "dashed")
    p <- p + scale_x_continuous(lim = c(0, 24))
    p <- p + facet_wrap(~ ITER, ncol = 4, scales = "free_y")
    p
  }

  d2b.int <- apply(d2b.p, 2, function(p, x) chisq.interv(p, x)$times, x = t.ph)
  d2b.int.df <- data.frame(
    ITER = rep(1:niter, deach = dim(d2b.int)[1]),
    INT = as.vector(d2b.int)
  )
  d2b.df <- data.frame(
    ITER = rep(1:niter, each = dim(d2b)[1]),
    TIME = time.samp,
    DV = as.vector(d2b)
  )
  facetted.interval.plot("d2b")

  d3b.int <- apply(d3b.p, 2, function(p, x) optim.interv.dt(p, x)$times, x = t.ph)
  d3b.int.df <- data.frame(
    ITER = rep(1:niter, deach = dim(d3b.int)[1]),
    INT = as.vector(d3b.int)
  )
  d3b.df <- data.frame(
    ITER = rep(1:niter, each = dim(d3b)[1]),
    TIME = time.samp,
    DV = as.vector(d3b)
  )
  facetted.interval.plot("d3b")

  d1a.int <- apply(d1a.p, 2, function(p, x) optim.interv.dt(p, x)$times, x = t.ph)
  d1a.int.df <- data.frame(
    ITER = rep(1:niter, deach = dim(d1a.int)[1]),
    INT = as.vector(d1a.int)
  )
  d1a.df <- data.frame(
    ITER = rep(1:niter, each = dim(d1a)[1]),
    TIME = time.samp,
    DV = as.vector(d1a)
  )
  facetted.interval.plot("d1a")

  d2a.int <- apply(d2a.p, 2, function(p, x) optim.interv.dt(p, x)$times, x = t.ph)
  d2a.int.df <- data.frame(
    ITER = rep(1:niter, deach = dim(d2a.int)[1]),
    INT = as.vector(d2a.int)
  )
  d2a.df <- data.frame(
    ITER = rep(1:niter, each = dim(d2a)[1]),
    TIME = time.samp,
    DV = as.vector(d2a)
  )
  facetted.interval.plot("d2a")

  d3a.int <- apply(d3a.p, 2, function(p, x) optim.interv.dt(p, x)$times, x = t.ph)
  d3a.int.df <- data.frame(
    ITER = rep(1:niter, deach = dim(d3a.int)[1]),
    INT = as.vector(d3a.int)
  )
  d3a.df <- data.frame(
    ITER = rep(1:niter, each = dim(d3a)[1]),
    TIME = time.samp,
    DV = as.vector(d3a)
  )
  facetted.interval.plot("d3a")
