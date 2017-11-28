# Function containing the framework setup in the intial study design
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "optinterval"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
      reponame <- "optinterval"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "splines"
    }
    rm("wd")
  }

# Load packages
  library(GA)
  #library(ggplot2)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

# Source scripts to set up environment
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set basic parameters
  set.seed(256256)
  t1 <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
  niter <- 1000
  time.samp <- seq(0, 24, by = 0.05)
  sdev <- 0.15
  tlast <- max(t1)
  nobs <- 7
  t0 <- c(t1[1:(nobs+1)], tlast)

  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a.p <- c(-0.2, -0.4, 4)
  d1a <- pred.d1a(time.samp, d1a.p)
  times <- seq(0, tlast, length.out = nobs)
  err <- matrix(
    1 + rnorm(n = length(times)*niter, mean = 0, sd = sdev),
    nrow = length(times), ncol = niter
  )
  subd <- d1a[which(time.samp %in% times)]*err
  fit.par <- apply(subd, 2, function(x) {
    chisq.sumexp(optim.sumexp.hes(
      data.frame(time = times, conc = x), oral = T
    ))$par
  })
  t2 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0)$times, 1), tlast)
  )
  t3 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0, tmax = T)$times, 1), tlast)
  )
  auc <- data.frame(
    true = integrate(pred.d1a, 0, tlast, p = d1a.p)$value,
    user = auc.interv(t1, d1a.p, pred.d1a),
    optint = apply(t2, 2, function(x) auc.interv(x, d1a.p, pred.d1a)),
    optintwCmax = apply(t3, 2, function(x) auc.interv(x, d1a.p, pred.d1a))
  )
  cmax <- data.frame(
    true = pred.sumexp(d1a.p, tmax.sumexp(d1a.p)),
    user = max(pred.sumexp(d1a.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d1a.p, x))),
    optintwCmax = apply(t3, 2, function(x) max(pred.sumexp(d1a.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d1a.p),
    user = t1[which(pred.sumexp(d1a.p, t1) == cmax$user[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d1a.p, z) == y)][1]
    ),
    optintwCmax = mapply(cmax$optintwCmax, data.frame(t3),
      FUN = function(y, z) z[which(pred.sumexp(d1a.p, z) == y)][1]
    )
  )

  d1a.out <- list(
    par = d1a.p, fit.par = fit.par, t2 = t2, t3 = t3,
    auc = auc, cmax = cmax, tmax = tmax
  )
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(d1a.out, "d1a-narrow-dt15.rds")

# -----------------------------------------------------------------------------
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  d2a.p <- c(-0.06, -0.4, -0.6, log(exp(4)*0.15), log(exp(4)*0.85))
  d2a <- pred.d2a(time.samp, d2a.p)
  times <- seq(0, tlast, length.out = nobs)
  err <- matrix(
    1 + rnorm(n = length(times)*niter, mean = 0, sd = sdev),
    nrow = length(times), ncol = niter
  )
  subd <- d2a[which(time.samp %in% times)]*err
  fit.par <- apply(subd, 2, function(x) {
    chisq.sumexp(optim.sumexp.hes(
      data.frame(time = times, conc = x), oral = T
    ))$par
  })
  t2 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0)$times, 1), tlast)
  )
  t3 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0, tmax = T)$times, 1), tlast)
  )
  auc <- data.frame(
    true = integrate(pred.d2a, 0, tlast, p = d2a.p)$value,
    user = auc.interv(t1, d2a.p, pred.d2a),
    optint = apply(t2, 2, function(x) auc.interv(x, d2a.p, pred.d2a)),
    optintwCmax = apply(t3, 2, function(x) auc.interv(x, d2a.p, pred.d2a))
  )
  cmax <- data.frame(
    true = pred.sumexp(d2a.p, tmax.sumexp(d2a.p)),
    user = max(pred.sumexp(d2a.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d2a.p, x))),
    optintwCmax = apply(t3, 2, function(x) max(pred.sumexp(d2a.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d2a.p),
    user = t1[which(pred.sumexp(d2a.p, t1) == cmax$user[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d2a.p, z) == y)][1]
    ),
    optintwCmax = mapply(cmax$optintwCmax, data.frame(t3),
      FUN = function(y, z) z[which(pred.sumexp(d2a.p, z) == y)][1]
    )
  )

  d2a.out <- list(
    par = d2a.p, fit.par = fit.par, t2 = t2, t3 = t3,
    auc = auc, cmax = cmax, tmax = tmax
  )
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(d2a.out, "d2a-narrow-dt15.rds")

# -----------------------------------------------------------------------------
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }
  d3a.p <- c(-0.05, -0.25, -0.6, -1, log(exp(4)*0.10), log(exp(4)*0.35), log(exp(4)*0.55))
  d3a <- pred.d3a(time.samp, d3a.p)
  times <- seq(0, tlast, length.out = nobs)
  err <- matrix(
    1 + rnorm(n = length(times)*niter, mean = 0, sd = sdev),
    nrow = length(times), ncol = niter
  )
  subd <- d3a[which(time.samp %in% times)]*err
  fit.par <- apply(subd, 2, function(x) {
    chisq.sumexp(optim.sumexp.hes(
      data.frame(time = times, conc = x), oral = T
    ))$par
  })
  t2 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0)$times, 1), tlast)
  )
  t3 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0, tmax = T)$times, 1), tlast)
  )
  auc <- data.frame(
    true = integrate(pred.d3a, 0, tlast, p = d3a.p)$value,
    user = auc.interv(t1, d3a.p, pred.d3a),
    optint = apply(t2, 2, function(x) auc.interv(x, d3a.p, pred.d3a)),
    optintwCmax = apply(t3, 2, function(x) auc.interv(x, d3a.p, pred.d3a))
  )
  cmax <- data.frame(
    true = pred.sumexp(d3a.p, tmax.sumexp(d3a.p)),
    user = max(pred.sumexp(d3a.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d3a.p, x))),
    optintwCmax = apply(t3, 2, function(x) max(pred.sumexp(d3a.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d3a.p),
    user = t1[which(pred.sumexp(d3a.p, t1) == cmax$user[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d3a.p, z) == y)][1]
    ),
    optintwCmax = mapply(cmax$optintwCmax, data.frame(t3),
      FUN = function(y, z) z[which(pred.sumexp(d3a.p, z) == y)][1]
    )
  )

  d3a.out <- list(
    par = d3a.p, fit.par = fit.par, t2 = t2, t3 = t3,
    auc = auc, cmax = cmax, tmax = tmax
  )
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(d3a.out, "d3a-narrow-dt15.rds")

# -----------------------------------------------------------------------------
  pred.d2b <- function(x, p) {
    exp(p[1]*x + p[3]) + exp(p[2]*x + p[4])
  }
  d2b.p <- c(-0.06, -0.5, 5, 6)
  d2b <- pred.d2b(time.samp, d2b.p)
  times <- seq(0, tlast, length.out = nobs)
  err <- matrix(
    1 + rnorm(n = length(times)*niter, mean = 0, sd = sdev),
    nrow = length(times), ncol = niter
  )
  subd <- d2b[which(time.samp %in% times)]*err
  fit.par <- apply(subd, 2, function(x) {
    chisq.sumexp(optim.sumexp.hes(
      data.frame(time = times, conc = x), oral = T
    ))$par
  })
  t2 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0)$times, 1), tlast)
  )
  t3 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0, tmax = T)$times, 1), tlast)
  )
  auc <- data.frame(
    true = integrate(pred.d2b, 0, tlast, p = d2b.p)$value,
    user = auc.interv(t1, d2b.p, pred.d2b),
    optint = apply(t2, 2, function(x) auc.interv(x, d2b.p, pred.d2b)),
    optintwCmax = apply(t3, 2, function(x) auc.interv(x, d2b.p, pred.d2b))
  )
  cmax <- data.frame(
    true = pred.sumexp(d2b.p, tmax.sumexp(d2b.p)),
    user = max(pred.sumexp(d2b.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d2b.p, x))),
    optintwCmax = apply(t3, 2, function(x) max(pred.sumexp(d2b.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d2b.p),
    user = t1[which(pred.sumexp(d2b.p, t1) == cmax$user[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d2b.p, z) == y)][1]
    ),
    optintwCmax = mapply(cmax$optintwCmax, data.frame(t3),
      FUN = function(y, z) z[which(pred.sumexp(d2b.p, z) == y)][1]
    )
  )

  d2b.out <- list(
    par = d2b.p, fit.par = fit.par, t2 = t2, t3 = t3,
    auc = auc, cmax = cmax, tmax = tmax
  )
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(d2b.out, "d2b-narrow-dt15.rds")

# -----------------------------------------------------------------------------
  pred.d3b <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) + exp(p[3]*x + p[6])
  }
  d3b.p <- c(-0.03, -0.15, -0.5, 3.8, 4.7, 6)
  d3b <- pred.d3b(time.samp, d3b.p)
  times <- seq(0, tlast, length.out = nobs)
  err <- matrix(
    1 + rnorm(n = length(times)*niter, mean = 0, sd = sdev),
    nrow = length(times), ncol = niter
  )
  subd <- d3b[which(time.samp %in% times)]*err
  fit.par <- apply(subd, 2, function(x) {
    chisq.sumexp(optim.sumexp.hes(
      data.frame(time = times, conc = x), oral = T
    ))$par
  })
  t2 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0)$times, 1), tlast)
  )
  t3 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv.dtmax(x, t0, tmax = T)$times, 1), tlast)
  )
  auc <- data.frame(
    true = integrate(pred.d3b, 0, tlast, p = d3b.p)$value,
    user = auc.interv(t1, d3b.p, pred.d3b),
    optint = apply(t2, 2, function(x) auc.interv(x, d3b.p, pred.d3b)),
    optintwCmax = apply(t3, 2, function(x) auc.interv(x, d3b.p, pred.d3b))
  )
  cmax <- data.frame(
    true = pred.sumexp(d3b.p, tmax.sumexp(d3b.p)),
    user = max(pred.sumexp(d3b.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d3b.p, x))),
    optintwCmax = apply(t3, 2, function(x) max(pred.sumexp(d3b.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d3b.p),
    user = t1[which(pred.sumexp(d3b.p, t1) == cmax$user[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d3b.p, z) == y)][1]
    ),
    optintwCmax = mapply(cmax$optintwCmax, data.frame(t3),
      FUN = function(y, z) z[which(pred.sumexp(d3b.p, z) == y)][1]
    )
  )

  d3b.out <- list(
    par = d3b.p, fit.par = fit.par, t2 = t2, t3 = t3,
    auc = auc, cmax = cmax, tmax = tmax
  )
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(d3b.out, "d3b-narrow-dt15.rds")
