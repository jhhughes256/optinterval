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
  niter <- 1000
  time.samp <- seq(0, 24, by = 0.05)
  sdev <- 0.05
  tlast <- 24
  nobs <- 9
  t1 <- c(0, 0.5, 1, 2, 3, 4, 5, 6, 8, 12, 16, 24)

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
    FUN = function(x) c(0, round(optim.interv.dt(x, t1)$times, 1), tlast)
  )

  auc <- data.frame(
    true = integrate(pred.d1a, 0, tlast, p = d1a.p)$value,
    user = auc.interv(t1, d1a.p, pred.d1a),
    optint = apply(t2, 2, function(x) auc.interv(x, d1a.p, pred.d1a))
  )
  cmax <- data.frame(
    true = pred.sumexp(d1a.p, tmax.sumexp(d1a.p)),
    user = max(pred.sumexp(d1a.p, t1)),
    optint = apply(t2, 2, function(x) max(pred.sumexp(d1a.p, x)))
  )
  tmax <- data.frame(
    true = tmax.sumexp(d1a.p),
    user = t1[which(pred.sumexp(d1a.p, t1) == cmax$user[1])],
    optint = mapply(cmax$optint, data.frame(t2),
      FUN = function(y, z) z[which(pred.sumexp(d1a.p, z) == y)][1]
    )
  )

  out <- list(
    par = par, fit.par = fit.par, t2 = t2,
    auc = auc, cmax = cmax, tmax = tmax
  )
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(out, "d1a-narrow-dt.rds")
