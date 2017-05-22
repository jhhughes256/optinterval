# Get better initial parameters using GA
# -----------------------------------------------------------------------------
# The use of genetic algorithms to determine initial parameters
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
# altered to work with GA AND optim, alter the ga argument when working with GA
  mle.sumexp <- function(par, x, y, ga = F) {
    z <- ifelse(ga, 1, -1)
    yhat <- pred.sumexp(par[-length(par)], x)
    sigma <- par[length(par)]
    loglik <- dnorm(y, yhat, abs(sigma), log = T)
    return(z*sum(loglik))
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

  lm1 <- unname(lm(log(y) ~ x)$coefficients)
  init1 <- c(lm1[2], lm1[1])

  ga.res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
    min = c(rep(lm1[2]*50, 2), lm1[1]-2, 0.001),
    max = c(rep(-0.01, 2), lm1[1]+2, 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation,
    maxiter = 100
  )
  plot.sumexp(ga.res@solution, data1)
  bfgs.res <- optim(
    ga.res@solution,
    mle.sumexp,
    method = "BFGS",
    x = x, y = y
  )
  plot.sumexp(bfgs.res$par, data1)

  plot.sumexp(ga.res@solution, data1)

  # does this method of choosing boundaries work for all cases?
  # made better by using n starting populations, but same total iterations?
  # made better by bigger pop, smaller iterations?

  ga.res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
    min = c(rep(lm1[2]*50, 2), lm1[1]-2, 0.001),
    max = c(rep(-0.01, 2), lm1[1]+2, 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation,
    popSize = 250,
    maxiter = 20
  )
  # testing done using same number of total patients (popSize * maxiter)
  # by increasing popSize and reducing maxiter, optimisation is more reliable

  for (i in 1:100) {
    pop.res[i] <- ga("real-valued", # type of GA to use
      mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
      min = c(rep(lm1[2]*50, 2), lm1[1]-2, 0.001),
      max = c(rep(-0.01, 2), lm1[1]+2, 1),
      selection = gareal_lrSelection,
      crossover = gareal_spCrossover,
      mutation = gareal_raMutation,
      maxiter = 25,
      popSize = 250
    )@fitnessValue
  }
  c(mean.res = mean(pop.res), sd.res = sd(pop.res),
  min.res = min(pop.res), median.res = median(pop.res), max.res = max(pop.res))
  # data1 100 runs
  # maxiter popSize  mean    sd   min median  max benchmark
  #     100      50  -103   110         -59           42.72
  #      50     100   -67    52         -49
  #      40     125   -51    30         -41
  #      25     200   -46    24         -39
  #      20     250   -48    27  -194   -41   -15
  #      10     500   -44    14   -95   -42   -14
  #       8     625   -50    18  -116   -46   -17
  #       5    1000   -61    20  -109   -61   -18
  #       4    1250   -66    21  -129   -67   -21

# -----------------------------------------------------------------------------
  # examples with other datasets

  x <- data2[which(data2[, 2] != 0), 1]
  y <- data2[which(data2[, 2] != 0), 2]

  lm1 <- unname(lm(log(y) ~ x)$coefficients)
  init1 <- c(lm1[2], lm1[1])

  ga.res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
    min = c(rep(lm1[2]*50, 2), rep(lm1[1]-2, 2), 0.001),
    max = c(rep(-0.01, 2), rep(lm1[1]+2, 2), 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation,
    maxiter = 100
  )
  plot.sumexp(ga.res@solution, data2)
  bfgs.res <- optim(
    ga.res@solution,
    mle.sumexp,
    method = "BFGS",
    x = x, y = y
  )
  plot.sumexp(bfgs.res$par, data2)

  x <- data3[which(data3[, 2] != 0), 1]
  y <- data3[which(data3[, 2] != 0), 2]

  lm1 <- unname(lm(log(y) ~ x)$coefficients)
  init1 <- c(lm1[2], lm1[1])

  ga.res <- ga("real-valued", # type of GA to use
    mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
    min = c(rep(lm1[2]*50, 3), rep(lm1[1]-2, 2), 0.001),
    max = c(rep(-0.01, 3), rep(lm1[1]+2, 2), 1),
    selection = gareal_lrSelection,
    crossover = gareal_spCrossover,
    mutation = gareal_raMutation,
    maxiter = 100
  )
  plot.sumexp(ga.res@solution, data3)
  bfgs.res <- optim(
    ga.res@solution,
    mle.sumexp,
    method = "BFGS",
    x = x, y = y
  )
  plot.sumexp(bfgs.res$par, data3)
