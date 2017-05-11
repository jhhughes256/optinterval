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
  mle.sumexp <- function(par, x, y) {
    yhat <- pred.sumexp(par, x)
    err <- y - yhat
    sigma <- sd(err)
    loglik <- dnorm(y, yhat, sigma, log = T)
    return(-1*sum(loglik))
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
    sub.y <- which(y == max(y)):length(y)
    lm.par <- unname(lm(log(y[sub.y]) ~ x[sub.y])$coefficients)

    for (i in 1:nexp) {
      if (i == 1) {
        init.par <- c(lm.par[2], lm.par[1])
      } else if (i == 2) {
        if (oral) {
          init.par <- c(lm.par[2]*c(0.8, 1.2), lm.par[1])
        } else {
          init.par <- c(lm.par[2]*c(0.8, 1.2), lm.par[1]*c(1, 1.2))
        }
      } else {
        if (oral) {
          s <- seq(1-(i-(oral+1))*0.05, 1+(i-(oral+1))*0.05, length.out = i-oral)
          init.par <- c(
            mean(optres$par[1:(i-1)]), optres$par[1:(i-1)],
            log(sum(exp(optres$par[i:(2*i-3)]))/(i-1)*s)
          )
        } else {
          init.par <- c(
            mean(optres$par[1:(i-1)]), optres$par[1:(i-1)],
            mean(optres$par[i:(2*i-2)]), optres$par[i:(2*i-2)]
          )
        }
      }
      optres <- optim(
        init.par,
        mle.sumexp,
        method = "L-BFGS-B",
        lower = rep(-Inf, i),
        upper = c(rep(-1/10^10, i), rep(Inf, i)),
        control = list(maxit = 500),
        x = x, y = y
      )
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
      if (p < 0.01) {
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
