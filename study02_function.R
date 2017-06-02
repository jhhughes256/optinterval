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
  source(paste(git.dir, reponame, "study_data.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set basic parameters
  nobs <- 8

  function(data, par, fn, nobs, t2, sdev = 0.05, tlast = 24) {
    t1 <- seq(0, tlast, length.out = nobs+1)
    e1 <- 1 + rnorm(n = length(t1), mean = 0, sd = sdev)
    subd <- data.frame(
      time = t1,
      conc = with(data, conc[time %in% t1])*e1
    )
    exp.par <- chisq.sumexp(optim.sumexp(subd, oral = T))$par
    int.t3 <- optim.interv(t1, exp.par)
    t3 <- c(0, round(int.t3$par, 1), tlast)

    t4.tmax <- tmax.sumexp(exp.par)
    int.t4 <- optim.interv(t1, exp.par, tmax = t4.tmax)
    nr.t4 <- c(0, round(int.t4$par, 1), t4.tmax, tlast)
    t4 <- nr.t4[order(nr.t4)]

    auc <- data.frame(
      true = integrate(fn, 0, tlast, p = par),
      sumexp = ,
      basic = ,
      user = ,
      optint = ,
      optintwCmax =
    )
    tmax <- data.frame(
      true = ,
      basic = ,
      user = ,
      optint = ,
      optintwCmax =
    )
    cmax <- data.frame(
      true = ,
      basic = ,
      user = ,
      optint = ,
      optintwCmax =
    )
  }



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
    true = pred.sumexp(c(-0.1, -0.4, 4), tmax.true)
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
