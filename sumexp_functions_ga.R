# Sum of exponential functions for sourcing
# -----------------------------------------------------------------------------
# The functions in order are:
#   pred.sumexp - gives sum of exponentials for given set of parameters
#   mle.sumexp - maximum likelihood estimation function to be minimised
#   optim.sumexp - determine parameters for differing numbers of exponentials
#   err.interv - trapezoidal error function to be minimised
#   optim.interv - determine intervals that give the smallest error
# -----------------------------------------------------------------------------
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
  mle.sumexp <- function(par, x, y, ga = F) {
    z <- ifelse(ga, 1, -1)
    yhat <- pred.sumexp(par[-length(par)], x)
    sigma <- par[length(par)]
    loglik <- dnorm(y, yhat, abs(sigma), log = T)
    return(z*sum(loglik))
  }

# Fit sum of exponentials to curve for different numbers of exponentials
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
        gares <- ga("real-valued",
          mle.sumexp, x = x, y = y, ga = T,
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i), 0.001),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i), 1),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 50,
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

# Chi-squared difference test
# Takes a list of optim results and gives the best optim result
  chisq.sumexp <- function(opt) {
    i <- 1
    for (j in 2:length(opt$par)) {
      degf <- length(opt$par[[j]]) - length(opt$par[[i]])
      x <- opt$value[[i]] - opt$value[[j]]
      p <- pchisq(x, degf, lower.tail = F)
      if (p < 0.05) {
        i <- i + 1
      }
    }
    return(sapply(opt, function(x) x[i]))
  }

# Trapezoidal error function for interval optimisation
  err.interv <- function(par, exp.par, tmin, tmax, theta) {
    times <- c(tmin, par, tmax)
    deltat <- diff(times)
    secd <- pred.sumexp(exp.par, theta, 2)  # d = 2 (second derivative)
    err <- abs(deltat^3*secd/12)
    sum(err)
  }

# Interval optimising function
  optim.interv <- function(times, fit.par) {
    x <- times[order(times)]
    init.par <- x[-c(1, length(x))]
    xmin <- min(x)
    xmax <- max(x)
    theta <- c(NULL)

    for (i in 1:(length(times)-1)) {
      theta[i] <- optim(
        mean(c(times[i], times[i+1])),
        function(t, x) -abs(pred.sumexp(x, t, 2)),
        method = "L-BFGS-B", control = c(maxit = 500),
        lower = times[i], upper = times[i+1], x = fit.par
      )$par
    }

    res <- optim(
      init.par,
      err.interv,
      method = "L-BFGS-B", control = c(maxit = 500),
      lower = xmin, upper = xmax,
      exp.par = fit.par, tmin = xmin, tmax = xmax, theta = theta
    )
  # Useful values to return are $par
  # $value is a sum of errors, errors would be more interesting when separated
    return(res)
  }
