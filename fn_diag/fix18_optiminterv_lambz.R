# Creating optimal time samples that provide accurate terminal slopes
# Two methods
# - optimise one less time point within optim.interv, non-optimised time point
#  would be placed between the penultimate and last parameters.
#   + placed in the dead centre
#   + placed in the log centre
#   + placed so that the AUC is optimised between the 3 points
#   + placed so that the lambdaz is optimised between the 3 points
# - fix the two points like with tmax, optimised so that when combined with tlast
#   they provide an accurate lambdaz (might just trend to lower boundary....)
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32", "C:/Users/Jim Hughes/Documents/GitRepos/optinterval")

    graphics.off()
    if (getwd() == wd[1] | getwd() == wd[5]) {
      git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
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

# Load libraries
  library(GA)
  library(ggplot2)

# Source scripts to set up environment
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set up
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }
  fit.par <- c(-0.05, -0.25, -0.6, -1, log(exp(4)*0.10), log(exp(4)*0.35), log(exp(4)*0.55))
  time.samp <- seq(0, 48, by = 0.1)
  obsdata <- data.frame(
    times = time.samp,
    dv = pred.sumexp(fit.par, time.samp)
  )
  trueauc <- integrate(pred.d3a, 0, 48, p = fit.par)$value + with(obsdata, tail(dv, 1))/0.05
  t1 <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)
  tlast <- 24
  nobs <- 7
  t0 <- c(t1[1:(nobs+1)], tlast)

# Typical situation
  t2 <- c(optim.interv.dtmax(fit.par, t0)$times, tlast)
  predata <- data.frame(
    times = t2,
    dv = pred.sumexp(fit.par, t2)
  )
  lmres2 <- lm(log(dv) ~ times, tail(predata, 3))
  lmresc2 <- unname(lmres2$coefficients)
  lmrest2 <- seq(round(tail(predata, 3)$times[1], 0), 48, by = 0.1)
  lmdata2 <- data.frame(
    times = lmrest2,
    dv = exp(predict(lmres2, data.frame(times = lmrest2)))
  )
  auc2 <- auc.interv.sumexp(t2, fit.par, log = T) + with(predata, tail(dv, 1))/-lmresc2[2]
  prop2 <- auc2/trueauc

  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_line(aes(x = times, y = dv),
    data = obsdata, size = 0.8, colour = "blue")
  p2 <- p2 + geom_point(aes(x = times, y = dv),
    data = predata, size = 1.5, colour = "red")
  p2 <- p2 + geom_line(aes(x = times, y = dv),
    data = lmdata2, size = 0.8, colour = "green4", linetype = "dashed")
  p2 <- p2 + scale_y_log10(lim = c(0.1, 20))
  p2

# New methods
  nobs <- 6
  t0 <- c(t1[1:(nobs+1)], tlast)
  topt <- c(optim.interv.dtmax(fit.par, t0)$times, tlast)

# Mean of final two
  t3 <- c(head(topt), tail(topt, 2)[1], mean(tail(topt, 2)), tail(topt, 1))
  predata <- data.frame(
    times = t3,
    dv = pred.sumexp(fit.par, t3)
  )
  lmres3 <- lm(log(dv) ~ times, tail(predata, 3))
  lmresc3 <- unname(lmres3$coefficients)
  lmrest3 <- seq(round(tail(predata, 3)$times[1], 0), 48, by = 0.1)
  lmdata3 <- data.frame(
    times = lmrest3,
    dv = exp(predict(lmres3, data.frame(times = lmrest3)))
  )
  auc3 <- auc.interv.sumexp(t3, fit.par, log = T) + with(predata, tail(dv, 1))/-lmresc3[2]
  prop3 <- auc3/trueauc

  p3 <- NULL
  p3 <- ggplot()
  p3 <- p3 + geom_line(aes(x = times, y = dv),
    data = obsdata, size = 0.8, colour = "blue")
  p3 <- p3 + geom_point(aes(x = times, y = dv),
    data = predata, size = 1.5, colour = "red")
  p3 <- p3 + geom_line(aes(x = times, y = dv),
    data = lmdata3, size = 0.8, colour = "green4", linetype = "dashed")
  p3 <- p3 + scale_y_log10(lim = c(0.1, 20))
  p3

# Geomean of final two
  t4 <- c(head(topt), tail(topt, 2)[1], exp(mean(log(tail(topt, 2)))), tail(topt, 1))
  predata <- data.frame(
    times = t4,
    dv = pred.sumexp(fit.par, t4)
  )
  lmres4 <- lm(log(dv) ~ times, tail(predata, 3))
  lmresc4 <- unname(lmres4$coefficients)
  lmrest4 <- seq(round(tail(predata, 3)$times[1], 0), 48, by = 0.1)
  lmdata4 <- data.frame(
    times = lmrest4,
    dv = exp(predict(lmres4, data.frame(times = lmrest4)))
  )
  auc4 <- auc.interv.sumexp(t4, fit.par, log = T) + with(predata, tail(dv, 1))/-lmresc4[2]
  prop4 <- auc4/trueauc

  p4 <- NULL
  p4 <- ggplot()
  p4 <- p4 + geom_line(aes(x = times, y = dv),
    data = obsdata, size = 0.8, colour = "blue")
  p4 <- p4 + geom_point(aes(x = times, y = dv),
    data = predata, size = 1.5, colour = "red")
  p4 <- p4 + geom_line(aes(x = times, y = dv),
    data = lmdata4, size = 0.8, colour = "green4", linetype = "dashed")
  p4 <- p4 + scale_y_log10(lim = c(0.1, 20))
  p4

# optimised AUC
  tail.par <- c(tail(topt, 2)[1], mean(tail(topt, 2)), tail(topt, 1))
  t5 <- c(head(topt), optim.interv.dtmax(fit.par, tail.par)$times)
  
  predata <- data.frame(
    times = t5,
    dv = pred.sumexp(fit.par, t5)
  )
  lmres5 <- lm(log(dv) ~ times, tail(predata, 3))
  lmresc5 <- unname(lmres5$coefficients)
  lmrest5 <- seq(round(tail(predata, 3)$times[1], 0), 48, by = 0.1)
  lmdata5 <- data.frame(
    times = lmrest5,
    dv = exp(predict(lmres5, data.frame(times = lmrest5)))
  )
  auc5 <- auc.interv.sumexp(t5, fit.par, log = T) + with(predata, tail(dv, 1))/-lmresc5[2]
  prop5 <- auc5/trueauc
  
  p5 <- NULL
  p5 <- ggplot()
  p5 <- p5 + geom_line(aes(x = times, y = dv),
                       data = obsdata, size = 0.8, colour = "blue")
  p5 <- p5 + geom_point(aes(x = times, y = dv),
                        data = predata, size = 1.5, colour = "red")
  p5 <- p5 + geom_line(aes(x = times, y = dv),
                       data = lmdata5, size = 0.8, colour = "green4", linetype = "dashed")
  p5 <- p5 + scale_y_log10(lim = c(0.1, 20))
  p5
  
# optimised lambdaz
  topt.lamb <- optim(
    tail.par[2],
    function(p, tfirst, tlast, fit) {
      m <- fit[1:ceiling(length(fit)/2)]
      times <- c(tfirst, p, tlast)
      pred <- data.frame(
        time = times,
        dv = pred.sumexp(fit.par, times)
      )
      lmres <- lm(log(dv) ~ times, pred)$coefficients
      err <- sqrt(diff(c(max(m), lmres[2]))^2)
      return(err)
    }, 
    method = "L-BFGS-B", hessian = T,
    lower = tail.par[1], upper = tail.par[3],
    tfirst = tail.par[1], tlast = tail.par, fit = fit.par
  )$par
  t6 <- c(head(topt), tail.par[1], topt.lamb, tail.par[3])
  
  predata <- data.frame(
    times = t6,
    dv = pred.sumexp(fit.par, t6)
  )
  lmres6 <- lm(log(dv) ~ times, tail(predata, 3))
  lmresc6 <- unname(lmres6$coefficients)
  lmrest6 <- seq(round(tail(predata, 3)$times[1], 0), 48, by = 0.1)
  lmdata6 <- data.frame(
    times = lmrest6,
    dv = exp(predict(lmres6, data.frame(times = lmrest6)))
  )
  auc6 <- auc.interv.sumexp(t6, fit.par, log = T) + with(predata, tail(dv, 1))/-lmresc6[2]
  prop6 <- auc6/trueauc
  
  p6 <- NULL
  p6 <- ggplot()
  p6 <- p6 + geom_line(aes(x = times, y = dv),
                       data = obsdata, size = 0.8, colour = "blue")
  p6 <- p6 + geom_point(aes(x = times, y = dv),
                        data = predata, size = 1.5, colour = "red")
  p6 <- p6 + geom_line(aes(x = times, y = dv),
                       data = lmdata6, size = 0.8, colour = "green4", linetype = "dashed")
  p6 <- p6 + scale_y_log10(lim = c(0.1, 20))
  p6
