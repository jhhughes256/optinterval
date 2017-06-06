# Testing the use of genetic algorithm for optimisation
# -----------------------------------------------------------------------------
# Genetic algorithms explained @ the following links
# Basics: http://www.obitko.com/tutorials/genetic-algorithms/index.php
# GA Package: https://www.jstatsoft.org/article/view/v053i04/v53i04.pdf
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
  library(ggplot2)
  library(GA)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Load in functions
# Sum of exponentials predicted concentrations function
  pred.sumexp <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    for (i in 1:n) {
      if (i == 1) y <- x[i]^d*exp(x[i]*t + x[n+i])
      else if (i != n | a == 0) y <- y + x[i]^d*exp(x[i]*t + x[n+i])
      else if (a == 1) y <- y - x[i]^d*exp(x[i]*t)*sum(exp(x[(n+1):(2*n-1)]))
    }
    return(y)
  }

# Maximum likelihood estimation function for parameter optimisation
# altered to work with GA
  mle.sumexp <- function(par, x, y) {
    yhat <- pred.sumexp(par[-length(par)], x)
    sigma <- par[length(par)]
    loglik <- dnorm(y, yhat, abs(yhat*sigma), log = T)
    return(sum(loglik))
  }

  plot.sumexp <- function(res, data) {
    plotdata <- data.frame(
      time = data$time,
      cobs = data$conc,
      pred = pred.sumexp(res[-length(res)], data$time, 0)
    )
    ylim <- c(0, 1.1*max(plotdata$cobs))
    xlim <- c(0, max(plotdata$time))

    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Predicted and Observed vs. Time")
    plotobj <- plotobj + geom_point(aes(x = time, y = cobs))
    plotobj <- plotobj + geom_line(aes(x = time, y = pred), colour = "red")
    plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n", lim = ylim)
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
    plotobj
  }

# Load in data
  source(paste(git.dir, reponame, "sumexp_data.R", sep = "/"))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create datasets
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.05)
  data1 <- data.frame(
    time = time.samp,
    conc = onedata.abs$sumexp*err
  )
  data2 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp*err
  )
  data3 <- data.frame(
    time = time.samp,
    conc = twodata.abs$sumexp*err
  )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initial parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  x <- data1[which(data1[, 2] != 0), 1]
  y <- data1[which(data1[, 2] != 0), 2]

  res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y,  # maximum likelihood estimation function
    min = c(rep(-1, 2), 3, 0.001),
    max = c(rep(-0.01, 2), 5, 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation
  )

  x <- data2[which(data2[, 2] != 0), 1]
  y <- data2[which(data2[, 2] != 0), 2]

  res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y,  # maximum likelihood estimation function
    min = c(rep(-1, 2), rep(3, 2), 0.001),
    max = c(rep(-0.01, 2), rep(7, 2), 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation
  )

  x <- data3[which(data3[, 2] != 0), 1]
  y <- data3[which(data3[, 2] != 0), 2]

  res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y,  # maximum likelihood estimation function
    min = c(rep(-1, 3), rep(1, 2), 0.001),
    max = c(rep(-0.01, 3), rep(5, 2), 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation
  )

  all.res <- double(0)
  for (i in 1:50) {
    res <- ga("real-valued", # type of GA to use
      mle.sumexp, x = x, y = y,  # maximum likelihood estimation function
      min = c(rep(-1, 2), 3, 0.001),
      max = c(rep(-0.01, 2), 5, 1),
      selection = gareal_lrSelection,
      crossover = gareal_spCrossover,
      mutation = gareal_raMutation,
      maxiter = 300
    )@fitnessValue
    all.res[i] <- res
  }

  por.res <- double(0)
  for (i in 1:50) {
    res <- -Inf
    for (j in 1:5) {
      val <- ga("real-valued", # type of GA to use
        mle.sumexp, x = x, y = y,  # maximum likelihood estimation function
        min = c(rep(-1, 2), 3, 0.001),
        max = c(rep(-0.01, 2), 5, 1),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 60
      )@fitnessValue
      if (val > res) res <- val
    }
    por.res[i] <- res
  }

  c(mean.all = mean(all.res), sd.all = sd(all.res),
    mean.por = mean(por.res), sd.all = sd(all.res))

# Best for this type of problem seem to be
# lrSelection, spCrossover, raMutation
# GA operator cheatsheet
# selection: lr - linear rank | nlr - non-linear rank | rw - roulette wheel
# crossover: sp - single point | wa - whole arithmetic | la - local arithmetic |
#            blx - blend | laplace - laplacian arithmetic?
# mutation: ra - uniform random | nra - nonuniform random | pow - power? |
#           rs - random (around solution)

# ALSO better to do multiple smaller ga runs, than one big one
# Or a combination could be used?
