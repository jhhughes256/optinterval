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
    git.dir <- "E:/Hughes/Git"
    #git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    #git.dir <- "C:/Users/hugjh001/Documents"
    setwd(git.dir)
    reponame <- "optinterval"
  }

# Load and configure libraries
  library(ggplot2)
  library(GA)
  library(plyr)
  library(rbenchmark)
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
  data4 <- data.frame(
    time = time.samp,
    conc = onedata$sumexp*err
  )
  data5 <- data.frame(
    time = time.samp,
    conc = thrdata$sumexp*err
  )
  data6 <- data.frame(
    time = time.samp,
    conc = thrdata.abs$sumexp*err
  )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load in data from GA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # loop1, tested maxiter and popsize using data1 through 6
  loop.res1 <- readRDS(paste(reponame, "gatest_maxiter_popsize.rds", sep="/"))

  # loop2, tested      "       "      using data1 with varying sampling regimens
  # 1:25, 1:12*2, 1:13, c(2:7, 9, 11, 13, 17, 25), c(2:5, 7, 9, 25), (1:6*4-2)
  loop.res2 <- readRDS(paste(reponame, "gatest_sampling_regimen.rds"))

  # loop3, tested maxiter for popsize 250
  # using data from loop2
  loop.res3 <- readRDS(paste(reponame, "gatest_maxiter_fixedpop.rds"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Analyse
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# loop.res1
  loop.res3 <- llply(loop.res, function(x) {
    colnames(x) <- c("maxiter", "popSize", "mean", "sd", "min", "med", "max", "bench")
    x
  })





  ########## WARNING #########
  # THIS TAKES AN AGE TO RUN #
  # Written as a function so you don't accidentally run it! :O
    loop.iter <- c(10, 20, 40, 60, 80, 100)
    loop.res <- list(NULL)
    for (i in 1:6) {
      if (i == 1) data <- data1
      else if (i == 2) data <- data1[c(1:12*2),]
      else if (i == 3) data <- data1[c(1:13),]
      else if (i == 4) data <- data1[c(2:7, 9, 11, 13, 17, 25),]
      else if (i == 5) data <- data1[c(2:5, 7, 9, 25),]
      else if (i == 6) data <- data1[c(1:6*4-2),]
      x <- data[which(data[, 2] != 0), 1]
      y <- data[which(data[, 2] != 0), 2]
      lm1 <- unname(lm(log(y) ~ x)$coefficients)
      i.res <- matrix(nrow = 9, ncol = 8)
      for (j in 1:6) {
        j.iter <- loop.iter[j]
        j.pops <- 250
        j.res <- double(0)
        for (k in 1:1000) {
          j.res[k] <- ga("real-valued", # type of GA to use
            mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
            min = c(rep(lm1[2]*50, 2), lm1[1]-2, 0.001),
            max = c(rep(-0.01, 2), lm1[1]+2, 1),
            selection = gareal_lrSelection,
            crossover = gareal_spCrossover,
            mutation = gareal_raMutation,
            maxiter = j.iter,
            popSize = j.pops
          )@fitnessValue
        }
        j.bench <- benchmark({
          ga("real-valued", # type of GA to use
            mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
            min = c(rep(lm1[2]*50, 2), lm1[1]-2, 0.001),
            max = c(rep(-0.01, 2), lm1[1]+2, 1),
            selection = gareal_lrSelection,
            crossover = gareal_spCrossover,
            mutation = gareal_raMutation,
            maxiter = j.iter,
            popSize = j.pops
          )
        })
        i.res[j,] <- c(j.iter, j.pops, mean(j.res), sd(j.res),
          min(j.res), median(j.res), max(j.res), j.bench$elapsed)
      }
      loop.res[[i]] <- i.res
    }
