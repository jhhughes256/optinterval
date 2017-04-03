mlespline.fit <- function(par, x, y) {
  err <- y
  for (i in 1:length(par)) {
    if (i < length(par) - 1) {
      err <- err - par[i]*x[,i]
    } else if (i == length(par) - 1) {
      err <- err - par[i]
    }
  }
  loglik <- dnorm(err, 0, tail(par, 1), log = T)
  return(-1*sum(loglik))
}

mlespline.opt <- function(data, K, n) {
  ndata <- data.frame(x = data[,1], y = data[,2])
  lm.mod <- lm(y ~ bs(x, knots = c(K), degree = 2), data = ndata)
  init.par <- unname(pk.modlm$coefficients[c(2:5, 1)])
  suppressWarnings(
    optim(
      c(init.par, 1),
      mlespline.fit,
      method = "BFGS",
      x = bs(ndata$x, knots = c(K), degree = 2), y = ndata$y
    )
  )
}

mlespline.predict <- function(x, par, K, n) {
  bs.mat <- bs(x, knots = c(K), degree = n)
  y <- double(length(x))
  for (i in 1:length(par)) {
    if (i < length(par) - 1) {
      y <- y + par[i]*bs.mat[,i]
    } else if (i == length(par) - 1) {
      y <- y + par[i]
    }
  }
  out <- list(
    time = x,
    pred = y,
    knots = K,
    degree = n,
    coefficients = par[-length(par)],
    sigma = par[length(par)]
  )
  names(out$coefficients) <- paste0("beta", 1:(length(par) - 1))
  return(out)
}

mlespline <- function(data, K, n, predict.times = NULL) {
  if (is.null(predict.times)) predict.times <- data[,1]
  time <- data[,1]
  conc <- data[,2]
  spline <- mlespline.opt(data, K, n)
  return(mlespline.predict(predict.times, spline$par, K, n))
}

knot.func <- function(x) {
  y <- double(length(x))
  for (i in 1:length(x)) {
    y[i] <- sum(x[1:i])
  }
  return(y)
}

mlespline.knotfit <- function(par, data, n) {
  Chat <- mlespline(data, knot.func(par), n)$pred
  err <- Cobs-Chat
  sigma <- sd(err)
  loglik <- dnorm(Cobs, Chat, sigma, log = T)
  return(-1*sum(loglik))
}
