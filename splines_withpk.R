# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  K <- c(5, 10)  # knots
  # NOT onecomp - > PK! WIP
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  ERR <- 0.3	# Standard deviation of error
  dose <- 50	# mg
  times <- seq(from = 0,to = 24,by = 0.25)	# Time sequence for simulating concentrations
  sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)	# Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))*(1+rnorm(n = length(times),mean = 0,sd = ERR))
  sample.conc <- conc[times %in% sample.times]
  onecomp <- data.frame(time = sample.times, conc = sample.conc)
  plot(onecomp)

# ------------------------------------------------------------------------------
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
  u <- seq(0, 24, by = 0.1)  # test model at these times
  B <- data.frame(time = u)  # data.frame(times)
  Y <- predict(reg, newdata = B)  # predicted conc for the desired times
  lines(u, Y, lwd = 2, col = "red")

#  i.e. we have the following (nice) picture
#  But if we look at the output of the regression we get this
  summary(reg)

#  So, what can we do with those numbers ? First, assume know that we
#  consider only one knot (we have to start somewhere), and we
#  consider a b-spline interpolation of degree 1 (i.e. linear by parts).
  K <- c(5)
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 1), data = onecomp)
  u <- seq(0, 24, by = 0.1)
  B <- data.frame(time = u)
  Y <- predict(reg, newdata = B)
  dev.off()
  plot(onecomp)
  lines(u, Y, lwd = 2, col = "red")
  summary(reg)

#  Recall that b-splines work like that:
#  given knots: 0 <= t[1] <= ... <= t[m-1] <= 1
#  (we define splines on the unit interval)
#  then: b[j,0](t) = 1 if t[j] <= t <= t[j+1] otherwise b[j,0](t) = 0
#  for all: j = 0, 1, ..., m-1, m-2

#                           t - t[j]                       t[j+n+1] - t
#  while: b[j,n](t) =  ----------------- b[j,n-1](t) + --------------------- b[j+1,n-1](t)
#                        t[j+n] - t[j]                   t[j+n+1] - t[j+1]

#  or: b[j,n](t) = a1 + a2
#  for all: j = 0, ..., m-n-2
#  The code is something like that:
  Bfun <- function(x, j, n, K) {
    b <- 0
    a1 <- 0
    a2 <- 0
    if(((K[j+n+1] > K[j+1]) & (j+n<=length(K)) & (n > 0)) == TRUE) {
      a2 <- (K[j+n+1] - x)/(K[j+n+1] - K[j+1]) * Bfun(x, j+1, n-1, K)
    }
    if(((K[j+n] > K[j]) & (n > 0)) == TRUE) {
      a1 <- (x - K[j])/(K[j+n] - K[j]) * Bfun(x, j, n-1, K)
    }
    if(n==0) {
      b <- ((x > K[j]) & (x <= K[j+1]))*1  # uncertain how this line works
    }
    if(n>0) {
      b <- a1 + a2
    }
    return(b)
  }

#  So, for instance, for splines of degree 1, we have
  u <- seq(0, 1, by = 0.01)
  plot(u, Bfun(u, 1, 1, c(0, 0.4, 1, 1)),
		lwd = 2, col = "red", type = "l", ylim = c(0, 1))
  lines(u, Bfun(u, 2, 1, c(0, 0.4, 1, 1)), lwd = 2, col = "blue")
  abline(v = c(0, 0.4, 1), lty = 2)

#  and for splines of degree 2, the basis is
  u <- seq(0, 1, by = 0.01)
  plot(u, Bfun(u, 1, 2, c(0 ,0, 0.4, 1, 1, 1)),
		lwd = 2, col = "red", type = "l", ylim = c(0, 1))
  lines(u, Bfun(u, 2, 2, c(0, 0, 0.4, 1, 1, 1)), lwd = 2, col = "blue")
  lines(u, Bfun(u, 3, 2, c(0, 0, 0.4, 1, 1, 1)), lwd = 2, col = "green")
  abline(v = c(0, 0.4, 1), lty = 2)

#  …etc. Note that I need to duplicate sometimes starting and end points (but it
#  should be possible to fix it in the function).
#  So, how do we use that, here ? Actually, there are two steps:

# 1. We get from [x[min],x[max]] to the unit interval (using a simple affine
#    transformation)
# 2. We consider Yhat[i] = Bhat[0] + Bhat[1]b[1,1](X[i]) + Bhat[2]b[2,1](X[i])

# Here, based on the graph above (with the basis function), note that we can use
	plot(onecomp)
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[2]*u0 + reg$coefficients[1]
  x1 <- seq(min(onecomp$time), K, length = length(u0))
  lines(x1, v, col = "green", lwd = 2)
  u0 <- seq(0, 1, by = 0.01)
  v <- (reg$coefficients[3] - reg$coefficients[2])*u0 +
		reg$coefficients[1] +
		reg$coefficients[2]
  x2 <- seq(K, max(onecomp$time), length = length(u0))
  lines(x2, v, col = "blue", lwd = 2)
#	which gives us exactly the graph we obtained previously.
# But we can also consider
  plot(onecomp)
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 1), data = onecomp)
  k <- (K - min(onecomp$time))/(max(onecomp$time) - min(onecomp$time))
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[1] +
    reg$coefficients[2]*Bfun(u0, 1, 1, c(0, k, 1, 1)) +
    reg$coefficients[3]*Bfun(u0, 2, 1, c(0, k, 1, 1))
  lines(x = min(onecomp$time) + u0*(max(onecomp$time) - min(onecomp$time)), y = v,
	 	col = "purple", lwd = 2)
  abline(v = K, lty = 2, col = "red")

# So, we should be able to try with two knots (but we keep it linear, so far)
  K <- c(5, 10)
  plot(onecomp)
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 1), data = onecomp)
  u <- seq(0, 24, by = 0.1)
  B <- data.frame(time = u)
  Y <- predict(reg, newdata = B)
  lines(u, Y, lwd = 2, col = "red")
  abline(v = K, lty = 2, col = "red")
# First, we can plot our basis functions, with two knots,
  u <- seq(0, 1, by = 0.01)
  plot(u, Bfun(u, 1, 1, c(0, 0.4, 0.7, 1)),
		lwd = 2, col = "red", type = "l", ylim = c(0, 1))
  lines(u, Bfun(u, 2, 1, c(0, 0.4, 0.7, 1, 1)), lwd = 2, col = "blue")
  lines(u, Bfun(u, 3, 1, c(0, 0.4, 0.7, 1, 1)), lwd = 2, col = "green")
  abline(v = c(0, 0.4, 0.7, 1), lty = 2)

# so we can use those functions here, as we did before,
  plot(onecomp)
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 1), data = onecomp)
  k <- (K - min(onecomp$time))/(max(onecomp$time) - min(onecomp$time))
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[1] +
    reg$coefficients[2]*Bfun(u0, 1, 1, c(0, k, 1, 1)) +
    reg$coefficients[3]*Bfun(u0, 2, 1, c(0, k, 1, 1)) +
    reg$coefficients[4]*Bfun(u0, 3, 1, c(0, k, 1, 1))
  lines(x = min(onecomp$time)+u0*(max(onecomp$time)-min(onecomp$time)), y = v,
		col = "red", lwd = 2)
  abline(v = K, lty = 2, col = "red")

# Great, it looks promising…. Let us look finally at the case we have two knots,
# and some quadratic splines. Here, with two knots, the basis is
  u <- seq(0, 1, by = 0.01)
  plot(u, Bfun(u, 1, 2, c(0, 0, 0.4, 0.7, 1, 1, 1)),
		lwd = 2, col = "red", type = "l", ylim = c(0, 1))
  lines(u, Bfun(u, 2, 2, c(0, 0, 0.4, 0.7, 1, 1, 1)), lwd = 2, col = "blue")
  lines(u, Bfun(u, 3, 2, c(0, 0, 0.4, 0.7, 1, 1, 1)), lwd = 2, col = "green")
  lines(u, Bfun(u, 4, 2, c(0, 0, 0.4, 0.7, 1, 1, 1)), lwd = 2, col = "orange")
  abline(v = c(0, 0.4, 0.7, 1), lty = 2)

# so if we just rewrite our previous function, we have
  plot(onecomp)
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
  k <- (K - min(onecomp$time))/(max(onecomp$time) - min(onecomp$time))
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[1] +
    reg$coefficients[2]*Bfun(u0, 1, 2, c(0, 0, k, 1, 1, 1)) +
    reg$coefficients[3]*Bfun(u0, 2, 2, c(0, 0, k, 1, 1, 1)) +
    reg$coefficients[4]*Bfun(u0, 3, 2, c(0, 0, k, 1, 1, 1)) +
    reg$coefficients[5]*Bfun(u0, 4, 2, c(0, 0, k, 1, 1, 1))
  lines(x = min(onecomp$time) + u0*(max(onecomp$time) - min(onecomp$time)), y = v,
	  col = "purple", lwd = 2)
  abline(v = K, lty = 2, col = "red")

# And just one final comment: how do I optimally
# choose my knots ?… A simple idea can be to consider all possible
# knots, and consider the model which gives us the residuals with the
# smaller variance,
  vk <- seq(0.05, 0.95, by = 0.05)
  SSR <- matrix(NA, length(vk))
  for(i in 1:(length(vk))){
    k <- vk[i]
    K <- min(onecomp$time) + k*(max(onecomp$time) - min(onecomp$time))
    reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
    SSR[i] <- sum(residuals(reg)^2)
  }
  plot(vk, SSR, type = "b", col = "blue")

# Here, the best model is obtained when we split 3/4-1/4…
