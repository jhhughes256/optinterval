# Determining the best fit out of the best sum of exponentials
# -----------------------------------------------------------------------------
# How to determine which is the best?
# Will AIC be useful?
# what magnitude of improvement is important?
# -----------------------------------------------------------------------------
# Setup environment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    graphics.off()
    #git.dir <- "E:/Hughes/Git"
    #git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "optinterval"
  }

# Load and configure libraries
  library(plyr)
  library(ggplot2)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Setup functions
  source(paste(git.dir, reponame, "sumexp_functions.R", sep = "/"))

  plot.sumexp <- function(res, data) {
    plotdata <- ldply(res$par, function(x) {
      data.frame(
        nexp = as.factor(ceiling(length(x)/2)),
        time = data$time,
        cobs = data$conc,
        pred = pred.sumexp(x, data$time, 0),
        objv = res$value[[ceiling(length(x)/2)]]
      )
    })
    levels(plotdata$nexp) <- paste(levels(plotdata$nexp), "exponent(s)")
    ylim <- c(0, 1.1*max(plotdata$cobs))
    xlim <- c(0, max(plotdata$time))

    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Comparison of Number of Exponents Used")
    plotobj <- plotobj + geom_point(aes(x = time, y = cobs))
    plotobj <- plotobj + geom_line(aes(x = time, y = pred), colour = "red")
    plotobj <- plotobj + geom_text(aes(x = median(time), y = max(cobs), label = objv))
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n", lim = ylim)
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
    plotobj <- plotobj + facet_wrap(~nexp)
    plotobj
  }
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Absorption Curve
  time.samp <- seq(0, 48, by = 1)
  absdata <- data.frame(
    time = time.samp,
    line1 = -0.2*time.samp + 4,
    line2 = -0.1*time.samp + 4
  )
  absdata$sumexp <- exp(absdata$line2) - exp(absdata$line1)
  #with(absdata, plot(time, sumexp))

# 2 Compartment Curve
  twodata <- data.frame(
    time = time.samp,
    line1 = -0.5*time.samp + 6,
    line2 = -0.05*time.samp + 5
  )
  twodata$sumexp <- exp(twodata$line1) + exp(twodata$line2)
  #with(twodata, plot(time, log(sumexp)))

# 2 Compartment Absorption Curve
  abstwodata <- data.frame(
    time = time.samp,
    line1 = -0.15*time.samp + log(exp(4)*0.9),
    line2 = -0.02*time.samp + log(exp(4)*0.1),
    line3 = -0.2*time.samp + 4
  )
  abstwodata$sumexp <- exp(abstwodata$line1) + exp(abstwodata$line2) - exp(abstwodata$line3)
  with(abstwodata, plot(time, sumexp))

# 3 Compartment Curve
  threedata <- data.frame(
    time = time.samp,
    line1 = -1*time.samp + 7,
    line2 = -0.2*time.samp + 5,
    line3 = -0.01*time.samp + 3
  )
  threedata$sumexp <- exp(threedata$line1) + exp(threedata$line2) + exp(threedata$line3)
  #with(threedata, plot(time, log(sumexp)))

# Create datasets
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.2)
  data1 <- data.frame(
    time = time.samp,
    conc = absdata$sumexp*err
  )
  res1 <- optim.sumexp(data1, oral = T, nexp = 4)
  data2 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp*err
  )
  res2 <- optim.sumexp(data2, oral = F, nexp = 4)
  data3 <- data.frame(
    time = time.samp,
    conc = abstwodata$sumexp#*err
  )
  res3 <- optim.sumexp(data3, oral = T, nexp = 4)
  data4 <- data.frame(
    time = time.samp,
    conc = threedata$sumexp*err
  )
  res4 <- optim.sumexp(data4, oral = F, nexp = 4)
  plot.sumexp(res1, data1)
  plot.sumexp(res2, data2)
  plot.sumexp(res3, data3)
  plot.sumexp(res4, data4)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Workflow for functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# First function that runs is optim.sumexp
# This gives a list telling us the parameters and their corresponding fitting
# values and whether they converged successfully
# This will then need to have the best number of exponentials chosen
# But how do you determine this?
  temp <- unlist(res1$value)
  which(temp == min(temp))
  # This gives the "best" curve, but four exponentials is unnecessary with such
  # little improvement!
  which(signif(temp, 5) == min(signif(temp, 5)))
  # Now we get multiple answers, but at least it contains the best one
  min(which(signif(temp, 5) == min(signif(temp, 5))))
  # This statement provides us with the "correct" answer

# So does this work with our other example? YES
  temp <- unlist(res2$value)
  min(which(signif(temp, 5) == min(signif(temp, 5))))

# This is not a very good way of determining the best model though...
# So lets use AIC?
# AIC = 2k - 2ln(objv)
# where k = number of parameters
  AIC.sumexp <- function(x) {
    AIC <- c(NULL)
    for (i in 1:length(x$par)) {
      k <- length(x$par[[i]])
      objv <- x$value[[i]]
      AIC[i] <- 2*k - 2*log(objv)
    }
    return(AIC)
  }
  c(AIC.sumexp(res1), AIC.sumexp(res2), AIC.sumexp(res3), AIC.sumexp(res4))
  # AIC is too strict for differentiating between the sums of exponentials

# How about BIC?
# BIC = kln(n) - 2ln(objv)
# where k = number of parameters, n = number of data points
  BIC.sumexp <- function(x, data) {
    n <- length(data[,1])
    BIC <- c(NULL)
    for (i in 1:length(x$par)) {
      k <- length(x$par[[i]])
      objv <- x$value[[i]]
      BIC[i] <- k*log(n) - 2*log(objv)
    }
    return(BIC)
  }
  BIC.sumexp(res1, data1)
  # Also too strict

# As these are nested models... how about using chi distributions?
# Covariate models are nested models.. surely this is similar?
  chisq.sumexp <- function(opt) {
    i <- 1
    for (j in 2:length(opt$par)) {
      degf <- length(opt$par[[j]]) - length(opt$par[[i]])
      x <- opt$value[[i]] - opt$value[[j]]
      p <- pchisq(x, degf, lower.tail = F)
      if (p < 0.01) {
        i <- i + 1
      }
    }
    return(sapply(opt, function(x) x[i]))
  }
  chisq.sumexp(res1)
  # Seems to work well for all cases
  #chisq.sumexp(res2)
  #chisq.sumexp(res3)
  #chisq.sumexp(res4)

# For chi distributions the better model needs to be x amount of points better
  optim(9.5, function(x) abs(0.01 - pchisq(x, 2, lower.tail = F)),
    method = "Brent", lower = 9, upper = 10)
  # with 2 degrees of freedom
  # for every exponential added, two parameters are added

  optim(6.5, function(x) abs(0.01 - pchisq(x, 1, lower.tail = F)),
    method = "Brent", lower = 6, upper = 7)
  # with 1 degree of freedom
  # only occurs for absorption models comparing 1 exp to 2 exp
