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
    #git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "optinterval"
  }

# Load and configure libraries
  library(ggplot2)
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
  mle.sumexp <- function(par, x, y) {
    yhat <- pred.sumexp(par[-length(par)], x)
    sigma <- abs(par[length(par)])
    loglik <- dnorm(y, yhat, sigma, log = T)
    return(-1*sum(loglik))
  }

  plot.sumexp <- function(res, data) {
    plotdata <- data.frame(
      time = data$time,
      cobs = data$conc,
      pred = pred.sumexp(res$par[-length(res$par)], data$time, 0),
      objv = res$value
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
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.1)
  data1 <- data.frame(
    time = time.samp,
    conc = onedata.abs$sumexp*err
  )
  data2 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp*err
  )
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.1)
  data3 <- data.frame(
    time = time.samp,
    conc = twodata.abs$sumexp*err
  )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initial parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tail.fit <- function(x, y) {
    all.n <- double(0)
    i <- 1
    max.bad <- 1
    repeat {
      n <- 3
      b.r2 <- 0
      bad <- 0
      repeat {
        r2 <- summary(lm(log(tail(y, n)) ~ tail(x, n)))$adj.r.squared
        if (b.r2 < r2) {
          b.r2 <- r2
          bad <- 0
        } else {
          bad <- bad + 1
        }
        if (bad == max.bad) {
          n <- n - bad
          break
        }
        n <- n + 1
      }
      if (i == 1) {
        all.n[i] <- n
        i <- i + 1
      } else {
        if (all.n[i-1] != n) {
          all.n[i] <- n
          i <- i + 1
        }
      }
      if ((max.bad + all.n[i-1]) >= length(y)) break
      max.bad <- max.bad + 1
    }
    all.n
  }
# Determine initial parameters for one exponential
  x <- data1[which(data1[, 2] != 0), 1]
  y <- data1[which(data1[, 2] != 0), 2]
  n1 <- tail.fit(x, y)[1]
  lm.par1 <- unname(lm(tail(log(y), n1) ~ tail(x, n1))$coefficients)
  init.par1 <- c(lm.par1[2], lm.par1[1])

# Determine for two exponential
# Use init.par1 because optim result will be suboptimal for absorption
  temp <- pred.sumexp(init.par1, time.samp) - data1[,2]
  if(any(temp <= 0.1)) temp <- temp[1:(which(temp <= 0.1)[1] - 1)]
  lm.par2 <- unname(lm(log(temp) ~ x[1:length(temp)])$coefficients)
  init.par2 <- c(lm.par1[2], lm.par2[2], lm.par1[1])

# Determine for three exponential
# Use optim result now
  optres1 <- optim(
    c(init.par2, 0.01),
    mle.sumexp,
    method = "L-BFGS-B",
    lower = c(rep(-Inf, 3), 0.0001),
    upper = c(rep(-1/10^10, 2), rep(Inf, 2)),
    control = list(maxit = 500),
    x = x, y = y
  )
  #plot.sumexp(optres2, data3)

  t1 <- length(x) - n1
  n2 <- tail.fit(x[1:t1], y[1:t1])
  lm.par3 <- unname(lm(log(tail(y[1:t1], n2[2])) ~ tail(x[1:t1], n2[2]))$coefficients)
  abs.y <- log(-(data3[,2] - pred.sumexp(c(lm.par3[2], lm.par3[1]), time.samp) - pred.sumexp(c(lm.par4[2], lm.par4[1]), time.samp)))
  lm.par5 <- unname(lm(abs.y ~ c(0, x))$coefficients)
  init.par3 <- c(lm.par3[2], lm.par4[2], optres1$par[2], lm.par3[1], lm.par4[1])

  optres2 <- optim(
    c(init.par3, 0.01),
    mle.sumexp,
    method = "L-BFGS-B",
    lower = c(rep(-Inf, 5), 0.0001),
    upper = c(rep(-1/10^10, 3), rep(Inf, 3)),
    control = list(maxit = 500),
    x = x, y = y
  )

  plot.sumexp(optres2, data3)
