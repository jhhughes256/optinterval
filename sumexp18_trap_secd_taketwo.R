# New method to determine secd in the error formula
# Aim: Design an algorithm that selects the correct theta always
# -----------------------------------------------------------------------------
# Incorporates ga functions from sumexp16
# Question still remains:
#  - what about linear-up log-down?
# -----------------------------------------------------------------------------
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
  #library(ggplot2)
  library(GA)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

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

# Trapezoidal error function for interval optimisation
  err.interv <- function(par, exp.par, tmin, tmax, a = F) {
    times <- c(tmin, par, tmax)
    deltat <- diff(times)
    if (a) {
      all.secd <- abs(pred.sumexp(exp.par, times, 2))
      secd <- c(NULL)
      for (i in 1:(length(times)-1)) {
        secd[i] <- all.secd[which(all.secd == max(all.secd[c(i, i + 1)]))]
      }
    } else {
      secd <- pred.sumexp(exp.par, times[-length(times)], 2)
    }
    err <- deltat^3*secd/12
    sum(err^2)
  }

# Interval optimising function
  optim.interv <- function(times, fit.par) {
    x <- times[order(times)]
    init.par <- x[-c(1, length(x))]
    xmin <- min(x)
    xmax <- max(x)
    absorp <- ifelse((length(fit.par) %% 2) != 0, T, F)
    res <- optim(
      init.par,
      err.interv,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = fit.par, tmin = xmin, tmax = xmax, a = absorp
    )
    return(res)
  }

# -----------------------------------------------------------------------------
# Set parameters
  nobs <- 12
  tmin <- 0
  tmax <- 24

# Set time sequence that is not far from the truth
# -2.5 is an arbitrary number, but does this relate in anyway to the exponentials?
  time.seq <- c(0, exp(seq(-2.5, 0, length.out = nobs))*tmax)

# Optimise
  optim.interv(time.seq, c(-0.5, -0.05, 6, 5))  # 2 comp example
  optim.interv(time.seq, c(-0.1, -0.2, 4))      # absorp example
