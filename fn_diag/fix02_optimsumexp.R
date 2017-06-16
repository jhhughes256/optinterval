# Fixing predsumexp as per Issue #8
# -----------------------------------------------------------------------------
# Set up libraries
  library(GA)
# -----------------------------------------------------------------------------
# Set up the functions
# pred.sumexp with for loop, changed for issue #8
# Identified as the best option in fix01_predsumexp.R
  pred.sumexp <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    m <- x[1:n]  # new
    ord <- order(m, decreasing = T)  # new
    p <- c(m[ord], x[(n+1):l])  # new
    for (i in 1:n) {
      if (i == 1) y <- p[i]^d*exp(p[i]*t + p[n+i])
      else if (i != n | a == 0) y <- y + p[i]^d*exp(p[i]*t + p[n+i])
      else if (a == 1) y <- y - p[i]^d*exp(p[i]*t)*sum(exp(p[(n+1):(2*n-1)]))
    }
    return(y)
  }

# Unchanged maximum likelihood estimation function for parameter optimisation
  mle.sumexp <- function(par, x, y, sigma, ga = F) {
    z <- ifelse(ga, 1, -1)
    yhat <- pred.sumexp(par, x)
    loglik <- dnorm(y, yhat, abs(yhat)*sigma, log = T)
    return(z*sum(loglik))
  }

# Updated optim.sumexp (using browser)
  optim.sumexp <- function(data, oral = F, nexp = 3) {
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
    opt.par <- list(NULL)
    opt.val <- list(NULL)
    opt.gra <- list(NULL)
    opt.con <- list(NULL)
    opt.mes <- list(NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    lmres <- unname(lm(log(y[lm.sub]) ~ x[lm.sub])$coefficients)
    for (i in 1:nexp) {
      if (i == 1 & !oral) {
        optres <- list(
          par = c(lmres[2], lmres[1]),
          value = mle.sumexp(unname(c(lmres[2], lmres[1])), x, y, 0.01),
          counts = NULL, convergence = 0, message = NULL
        )
      } else {
        gares <- ga("real-valued",
          mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
          min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
          max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
          selection = gareal_lrSelection,
          crossover = gareal_spCrossover,
          mutation = gareal_raMutation,
          maxiter = 50,
          popSize = 250,
          monitor = F
        )

        optres <- optim(
          gares@solution[1, ],
          mle.sumexp,
          method = "BFGS",
          x = x, y = y, sigma = 0.01
        )
      }
      slope.par <- optres$par[1:(i+oral)]
      slope.ord <- order(slope.par, decreasing = T)
      par.ord <- unname(c(slope.par[slope.ord], optres$par[(i+oral+1):length(optres$par)]))
      opt.par[[i]] <- par.ord
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

# -----------------------------------------------------------------------------
# Set up data for basic testing
  time.samp <- seq(0, 24, length.out = 9)
  d1a1.p <- c(-0.1, -0.4, 4)
  d1a1 <- data.frame(time = time.samp, conc = pred.sumexp(d1a1.p, time.samp))
  d3a1.p <- c(-0.001, -0.08, -0.5, -0.8, log(exp(4)*0.15), log(exp(4)*0.25), log(exp(4)*0.6))
  d3a1 <- data.frame(time = time.samp, conc = pred.sumexp(d3a1.p, time.samp))
  d3b1.p <- c(-0.01, -0.1, -0.4, 4.1, 4.7, 6)
  d3b1 <- data.frame(time = time.samp, conc = pred.sumexp(d3b1.p, time.samp))

  optim.sumexp(d1a1, oral = T)
  optim.sumexp(d3a1, oral = T)
  optim.sumexp(d3b1, oral = F)
  chisq.sumexp(optim.sumexp(d1a1, oral = T))
  chisq.sumexp(optim.sumexp(d3a1, oral = T))
  chisq.sumexp(optim.sumexp(d3b1, oral = F))

# -----------------------------------------------------------------------------
# Set up data for testing in regards to issue #8
  d1a2.p <- c(-0.2, -0.1, 4)
  d1a2 <- data.frame(time = time.samp, conc = pred.sumexp(d1a2.p, time.samp))
  d2a2.p <- c(-0.8, -0.01, -0.4, log(exp(4)*0.8), log(exp(4)*0.2))
  d2a2 <- data.frame(time = time.samp, conc = pred.sumexp(d2a2.p, time.samp))
  d3a2.p <- c(-0.001, -0.8, -0.08, -0.5, log(exp(4)*0.15), log(exp(4)*0.25), log(exp(4)*0.6))
  d3a2 <- data.frame(time = time.samp, conc = pred.sumexp(d3a2.p, time.samp))

  optim.sumexp(d1a2, oral = T)
  optim.sumexp(d2a2, oral = T)
  optim.sumexp(d3a2, oral = T)
  chisq.sumexp(optim.sumexp(d1a2, oral = T))
  chisq.sumexp(optim.sumexp(d2a2, oral = T))
  chisq.sumexp(optim.sumexp(d3a2, oral = T))
