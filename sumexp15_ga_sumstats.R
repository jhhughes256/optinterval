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

  ########## WARNING #########
  # THIS TAKES AN AGE TO RUN #
  loop.iter <- c(100, 50, 40, 25, 20, 10, 8, 5, 4)
  loop.pops <- c(50, 100, 125, 200, 250, 500, 625, 1000, 1250)
  loop.res <- list(NULL)
  for (i in 1:6) {
    if (i == 1) data <- data4
    else if (i == 2) data <- data1
    else if (i == 3) data <- data2
    else if (i == 4) data <- data3
    else if (i == 5) data <- data5
    else if (i == 6) data <- data6
    x <- data[which(data[, 2] != 0), 1]
    y <- data[which(data[, 2] != 0), 2]
    lm1 <- unname(lm(log(y) ~ x)$coefficients)
    i.res <- matrix(nrow = 9, ncol = 8)
    for (j in 1:9) {
      j.iter <- loop.iter[j]
      j.pops <- loop.pops[j]
      j.res <- double(0)
      for (k in 1:1000) {
        j.res[k] <- ga("real-valued", # type of GA to use
          mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
          min = c(rep(lm1[2]*50, ceiling((i+1)/2)), rep(lm1[1]-2, floor((i+1)/2)), 0.001),
          max = c(rep(-0.01, ceiling((i+1)/2)), rep(lm1[1]+2, floor((i+1)/2)), 1),
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
          min = c(rep(lm1[2]*50, ceiling((i+1)/2)), rep(lm1[1]-2, floor((i+1)/2)), 0.001),
          max = c(rep(-0.01, ceiling((i+1)/2)), rep(lm1[1]+2, floor((i+1)/2)), 1),
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
