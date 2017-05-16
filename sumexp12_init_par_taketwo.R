# Determining better methods of choosing intial parameters
# -----------------------------------------------------------------------------
# addresses issue #1
# optim was trending towards non-finite parameters resulting in error
# idea to determine initial parameters using sum of exponentials!
# -----------------------------------------------------------------------------
# Setup environment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    graphics.off()
    #git.dir <- "E:/Hughes/Git"
    git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    #git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "optinterval"
  }

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
  mle.sumexp <- function(par, x, y) {
    yhat <- pred.sumexp(par, x)
    err <- y - yhat
    sigma <- 0.01
    loglik <- dnorm(y, yhat, sigma, log = T)
    return(-1*sum(loglik))
  }

# Load in data
  source(paste(git.dir, reponame, "sumexp_data.R", sep = "/"))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create datasets
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data1 <- data.frame(
    time = time.samp,
    conc = onedata.abs$sumexp#*err
  )
  data2 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp#*err
  )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initial parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Determine initial parameters for one exponential
  x <- data1[which(data1[, 2] != 0), 1]
  y <- data1[which(data1[, 2] != 0), 2]
  sub.y <- which(y == max(y)):length(y)
  lm.par1 <- unname(lm(log(y[sub.y]) ~ x[sub.y])$coefficients)
  init.par1 <- c(lm.par1[2], lm.par1[1])

# Determine for two exponential
# Use init.par1 because optim result will be suboptimal for absorption
  temp <- pred.sumexp(init.par1, time.samp) - data1[,2]
  temp[temp <= 0.1] <- NA
  lm.par2 <- unname(lm(log(temp[!is.na(temp)]) ~ x[!is.na(temp)])$coefficients)
  init.par2 <- c(lm.par1[2], lm.par2[2], lm.par1[1])

# Determine for three exponential
# Use optim result now
  optim(
    init.par2,
    mle.sumexp,
    method = "L-BFGS-B",
    lower = rep(-Inf, 3),
    upper = c(rep(-1/10^10, 2), rep(Inf, 1)),
    control = list(maxit = 500),
    x = x, y = y
  )
