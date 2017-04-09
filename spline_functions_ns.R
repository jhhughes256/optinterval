mlespline.fit <- function(par, x, y) {
  err <- y
  for (i in 1:length(par)) {
    if (i < length(par) - 1) {
      err <- err - par[i]*x[, i]
    } else if (i == length(par) - 1) {
      err <- err - par[i]
    }
  }
  loglik <- dnorm(err, 0, tail(par, 1), log = T)
  return(-1*sum(loglik))
}

mlespline.opt <- function(data, K) {
  ndata <- data.frame(x = data[,1], y = data[,2])
  lm.mod <- lm(y ~ ns(x, knots = c(K)), data = ndata)
  init.par <- unname(lm.mod$coefficients[c(2:(length(knots) + 2), 1)])
  browser()
  suppressWarnings(
    optim(
      c(init.par, 1),
      mlespline.fit,
      method = "BFGS",
      x = ns(ndata$x, knots = c(K)), y = ndata$y
    )
  )
}

mlespline.predict <- function(x, par, K) {
  ns.mat <- ns(x, knots = c(K))
  y <- double(length(x))
  for (i in 1:length(par)) {
    if (i < length(par) - 1) {
      y <- y + par[i]*ns.mat[, i]
    } else if (i == length(par) - 1) {
      y <- y + par[i]
    }
  }
  out <- list(
    time = x,
    pred = y,
    knots = K,
    coefficients = par[-length(par)],
    sigma = par[length(par)]
  )
  names(out$coefficients) <- paste0("beta", 1:(length(par) - 1))
  return(out)
}

mlespline <- function(data, K, predict.times = NULL) {
  if (is.null(predict.times)) predict.times <- data[,1]
  time <- data[,1]
  conc <- data[,2]
  spline <- mlespline.opt(data, K)
  return(mlespline.predict(predict.times, spline$par, K))
}

knot.func <- function(x) {
  y <- double(length(x))
  for (i in 1:length(x)) {
    y[i] <- sum(x[1:i])
  }
  return(y)
}

mlespline.knotfit <- function(par, data) {
  Chat <- mlespline(data, knot.func(par))$pred
  Cobs <- data[, 2]
  err <- Cobs - Chat
  sigma <- sd(err)
  loglik <- dnorm(Cobs, Chat, sigma, log = T)
  return(-1*sum(loglik))
}
