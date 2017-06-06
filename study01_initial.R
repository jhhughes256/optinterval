# Script for the Study Design to highlight the algorithm
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
  source(paste(git.dir, reponame, "sumexp_functions_ga.R", sep = "/"))
  source(paste(git.dir, reponame, "sumexp_data.R", sep = "/"))

# -----------------------------------------------------------------------------
# Comparison of interval-choosing methods
# Method 1: Basic sequence
# Method 2: Smart sequence (as a 'PMx'er would decide)
# Method 3: Interval optimisation algorithm
# Method 4: Interval optimisation algorithm (+ Cmax)

# Set basic parameters
  nobs <- 8
  tlast <- 24
  sdev <- 0.05

# Create non-algorithm time samples
  d1a.t1 <- seq(0, tlast, length.out = nobs+1)
  d1a.t2 <- c(0, 0.5, 1, 3, 5, 7, 10, 16, 24)  # hand entered

# Set up datasets
  d1a.e1 <- 1 + rnorm(n = length(d1a.t1), mean = 0, sd = sdev)
  d1a1 <- data.frame(
    time = d1a.t1,
    conc = with(onedata.abs, sumexp[time %in% d1a.t1])*d1a.e1
  )
  d1a.e2 <- 1 + rnorm(n = length(d1a.t2), mean = 0, sd = sdev)
  d1a2 <- data.frame(
    time = d1a.t2,
    conc = with(onedata.abs, sumexp[time %in% d1a.t2])*d1a.e2
  )

# Create algorithm time samples
  d1a.sumexp1 <- chisq.sumexp(optim.sumexp(d1a1, oral = T))$par
  interv.res1 <- optim.interv(d1a.t1, d1a.sumexp1)
  d1a.t3 <- c(0, round(interv.res1$par, 1), tlast)

  tmax.sumexp <- function(fit.par, tlast = 24, res = 0.1) {
    times <- seq(0, tlast, by = res)
    yhat <- pred.sumexp(fit.par, times)
    return(times[which(yhat == max(yhat))])
  }
  d1a.tmax <- tmax.sumexp(d1a.sumexp1)
  interv.res2 <- optim.interv(d1a.t1, d1a.sumexp1, tmax = d1a.tmax)
  d1a.t4 <- c(0, round(interv.res2$par, 1), d1a.tmax, tlast)
  d1a.t4 <- d1a.t4[order(d1a.t4)]

  list(
    BasicSequence = d1a.t1,
    SmartSequence = d1a.t2,
    OptIntAlgorithm = d1a.t3,
    OptIntAlgorithm_wTmax = d1a.t4
  )

# Determine true auc_{0-tlast}, c_max and t_max
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a.auc <- list(
    true = integrate(pred.d1a, 0, tlast, p = c(-0.1, -0.4, 4)),
    sumexp = integrate(pred.d1a, 0, tlast, p = d1a.sumexp1)
  )
  d1a.tmax <- list(
    true = tmax.sumexp(c(-0.1, -0.4, 4))
  )
  d1a.cmax <- list(
    true = pred.sumexp(c(-0.1, -0.4, 4), d1a.tmax$true)
  )

  auc.interv <- function(times, fit.par, log = F) {
    C <- pred.d1a(times, fit.par)
    auc <- c(NULL)
    for (i in 2:length(C)) {
      h <- times[i] - times[i-1]
      dC <- C[i-1] - C[i]
      if (log & dC > 0) auc[i-1] <- dC*h/log(C[i-1]/C[i])
      else auc[i-1] <- (C[i-1] + C[i])*h/2
    }
    return(sum(auc))
  }

  d1a.auc$t1 <- auc.interv(d1a.t1, c(-0.1, -0.4, 4))
  d1a.auc$t2 <- auc.interv(d1a.t2, c(-0.1, -0.4, 4))
  d1a.auc$t3 <- auc.interv(d1a.t3, c(-0.1, -0.4, 4))
  d1a.auc$t4 <- auc.interv(d1a.t4, c(-0.1, -0.4, 4))

  d1a.auc$lnt1 <- auc.interv(d1a.t1, c(-0.1, -0.4, 4), log = T)
  d1a.auc$lnt2 <- auc.interv(d1a.t2, c(-0.1, -0.4, 4), log = T)
  d1a.auc$lnt3 <- auc.interv(d1a.t3, c(-0.1, -0.4, 4), log = T)
  d1a.auc$lnt4 <- auc.interv(d1a.t4, c(-0.1, -0.4, 4), log = T)

  d1a.cmax$t1 <- max(pred.sumexp(c(-0.1, -0.4, 4), d1a.t1))
  d1a.cmax$t2 <- max(pred.sumexp(c(-0.1, -0.4, 4), d1a.t2))
  d1a.cmax$t3 <- max(pred.sumexp(c(-0.1, -0.4, 4), d1a.t3))
  d1a.cmax$t4 <- max(pred.sumexp(c(-0.1, -0.4, 4), d1a.t4))

  d1a.tmax$t1 <- d1a.t1[which(pred.sumexp(c(-0.1, -0.4, 4), d1a.t1) == d1a.cmax.t1)]
  d1a.tmax$t2 <- d1a.t2[which(pred.sumexp(c(-0.1, -0.4, 4), d1a.t2) == d1a.cmax.t2)]
  d1a.tmax$t3 <- d1a.t3[which(pred.sumexp(c(-0.1, -0.4, 4), d1a.t3) == d1a.cmax.t3)]
  d1a.tmax$t4 <- d1a.t4[which(pred.sumexp(c(-0.1, -0.4, 4), d1a.t4) == d1a.cmax.t4)]
