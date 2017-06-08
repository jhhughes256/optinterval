# Function containing the framework setup in the intial study design
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    graphics.off()
    if (getwd() == "C:/Users/Jim Hughes/Documents") {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "optinterval"
    } else if (getwd() == "C:/Users/hugjh001/Documents") {
      git.dir <- getwd()
      reponame <- "optinterval"
    } else if (getwd() == "C:/Users/hugjh001/Desktop") {
      git.dir <- "E:/Hughes/Git"
      reponanme <- "splines"
    }
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
  niter <- 1000
  time.samp <- seq(0, 24, by = 0.05)
  sdev <- 0.05
  tlast <- 24
  nobs <- 8

  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a.p <- c(-0.2, -0.4, 4)
  d1a <- pred.d1a(time.samp, d1a.p)
  t1 <- seq(0, tlast, length.out = nobs+1)
  err <- matrix(
    1 + rnorm(n = length(t1)*niter, mean = 0, sd = sdev),
    nrow = length(t1), ncol = niter
  )
  subd <- d1a[which(time.samp %in% t1)]*err
  fit.par <- apply(subd, 2, function(x) {
    chisq.sumexp(optim.sumexp(
      data.frame(time = t1, conc = x), oral = T
    ))$par
  })
  t2 <- sapply(fit.par, simplify = "array",
    FUN = function(x) c(0, round(optim.interv(t1, x)$par, 1), tlast)
  )
  t3.tmax <- sapply(fit.par, tmax.sumexp)
  t3 <- mapply(fit.par, t3.tmax,
    FUN = function(x, y) {
      z <- c(0, round(optim.interv(t1, x, tmax = y)$par, 1), y, tlast)
      return(z[order(z)])
  })

  auc <- data.frame(
    true = integrate(pred.d1a, 0, tlast, p = d1a.p)$value,
    basic = auc.interv(t1, d1a.p, pred.d1a),
    optint = apply(t2, 2, function(x) auc.interv(x, d1a.p, pred.d1a)),
    optintwCmax = apply(t3, 2, function(x) auc.interv(x, d1a.p, pred.d1a))
  )
  cmax <- data.frame(
    true = pred.sumexp(d1a.p, tmax.sumexp(d1a.p)),
    basic = max(pred.sumexp(d1a.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d1a.p, x))),
    optintwCmax = apply(t3, 2, function(x) max(pred.sumexp(d1a.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d1a.p),
    basic = t1[which(pred.sumexp(d1a.p, t1) == cmax$basic[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d1a.p, z) == y)][1]
    ),
    optintwCmax = mapply(cmax$optintwCmax, data.frame(t3),
      FUN = function(y, z) z[which(pred.sumexp(d1a.p, z) == y)][1]
    )
  )
