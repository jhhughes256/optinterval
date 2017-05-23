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
    git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    #git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "optinterval"
  }

# Load and configure libraries
  library(ggplot2)
  library(GA)
  library(plyr)
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

  plot.sumexp <- function(res, data, oral) {
    plotdata <- ldply(res$par, function(x) {
      data.frame(
        nexp = as.factor(ceiling((length(x)-1)/2)),
        time = data$time,
        cobs = data$conc,
        pred = pred.sumexp(x[-length(x)], data$time, 0),
        objv = res$value[[ceiling((length(x)-1)/2)-oral]]
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
    plotobj <- plotobj + facet_wrap(~nexp, ncol = 1)
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
# Constructing function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  optim.sumexp <- function(data, oral = F, nexp = 3) {
    x <- data[which(data[, 2] != 0), 1]
    y <- data[which(data[, 2] != 0), 2]
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    lmres <- unname(lm(log(y) ~ x)$coefficients)

    for (i in 1:nexp) {
      if (i == 1 & !oral) {
        sigres <- optim(
          0.1,
          function(z) mle.sumexp(unname(c(lmres[2], lmres[1], z)), x, y),
          method = "Brent", lower = 0.001, upper = 10^10
        )
        optres <- list(
          par = unname(c(lmres[2], lmres[1], sigres$par)),
          value = sigres$value,
          counts = NULL, convergence = 0, message = NULL
        )
      } else {
        gares <- ga("real-valued", # type of GA to use
          mle.sumexp, x = x, y = y, ga = T, # maximum likelihood estimation function
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i), 0.001),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i), 1),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 20,
          popSize = 250
        )
        optres <- optim(
          gares@solution,
          mle.sumexp,
          method = "BFGS",
          x = x, y = y
        )
      }
      opt.par[[i]] <- optres$par
      opt.val[[i]] <- optres$value
      opt.gra[[i]] <- optres$counts
      opt.con[[i]] <- optres$convergence
      opt.mes[[i]] <- ifelse(is.null(optres$message),
        "NULL", optres$message)
    }

    res <- list(par = opt.par, value = opt.val, counts = opt.gra,
      convergence = opt.con, message = opt.mes)
    res
  }

  optim.res1 <- optim.sumexp(data1, oral = T)
  plot.sumexp(optim.res1, data1, oral = T)

  optim.res2 <- optim.sumexp(data2, oral = F)
  plot.sumexp(optim.res2, data2, oral = F)

  optim.res3 <- optim.sumexp(data3, oral = T)
  plot.sumexp(optim.res3, data3, oral = T)

  optim.res4 <- optim.sumexp(data4, oral = F)
  plot.sumexp(optim.res4, data4, oral = F)

  optim.res5 <- optim.sumexp(data5, oral = F)
  plot.sumexp(optim.res5, data5, oral = F)

  optim.res6 <- optim.sumexp(data6, oral = T)
  plot.sumexp(optim.res6, data6, oral = T)
