# Splines - Opening the Black Box
# by Arthur Charpentier (2010)
# URL: https://www.r-bloggers.com/splines-opening-the-black-box/

# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  Splines in regression is something which looks like a black box (or
#  maybe like some dishes you get when you travel away from home: it tastes
#  good, but you don’t what’s inside… even if you might have some clues,
#  you never know for sure*). With splines, it is the same: there are knots,
#  then we consider polynomial interpolations on parts between knots, and we
#  make sure that there is no discontinuity (on the prediction, but on the
#  derivative as well).

#  That sounds nice, but when you look at the output of the
#  regression… you got figures, but you barely see how to interpret
#  them… So let us have a look at the box, and I mean what is inside that
#  box…

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  K <- c(14, 20)  # knots
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 2), data = cars)
  u <- seq(4, 25, by = 0.1)  # test model at these speeds
  B <- data.frame(speed = u)  # data.frame(speeds)
  Y <- predict(reg, newdata = B)  # predicted dist for the desired speeds
  lines(u, Y, lwd = 2, col = "red")

#  i.e. we have the following (nice) picture
#  But if we look at the output of the regression we get this
  summary(reg)

#  So, what can we do with those numbers ? First, assume know that we
#  consider only one knot (we have to start somewhere), and we
#  consider a b-spline interpolation of degree 1 (i.e. linear by parts).
  K <- c(14)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 1), data = cars)
  u <- seq(4, 25, by = 0.1)
  B <- data.frame(speed = u)
  Y <- predict(reg, newdata = B)
  dev.off()
  plot(cars)
  lines(u, Y, lwd = 2, col = "red")
  summary(reg)

# Reading done from: https://cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf
# Usually variables for splines are:
#   m = the order of the spline
#   n = the degree of the spline
#   n = m - 1                         # fairly certain on this
#   N = number of internal knots
#   t = {t[i]|i ∈ Z}                    # the set of t[i] such that i is an integer
#   t = knot sequence, with individual terms being knots
#   number of knots in augmented knot set = N + 2m
#   augmented knot set =  internal knots + boundary knots * m
#   t[-(m-1)] = ... = t[0] <= t[1] <= ... <= t[N] <= t[N+1] = ... = t[N+m]
#   t[-(m-1)] through to t[0] are the lower boundary knots
#   t[N+1] through to t[N+m] are the upper boundary knots

# Unfortunately the terminology used in the code below does not match the
# mathematical functions provided by the blog, both are replicated as seen in
# the blog. The author uses j and n instead of i and j and appears to use t
# without subscript instead of x. It's all rather confusing... Bfun() appears
# to work though...

#  Recall that b-splines work like that:
#  given knots: 0 <= t[1] <= ... <= t[m-1] <= 1
#  (we define splines on the unit interval)
#  then: b[j,0](t) = 1 if t[j] <= t <= t[j+1] otherwise b[j,0](t) = 0
#  for all: j = 0, 1, ..., m-3, m-2

#                           t - t[j]                       t[j+n+1] - t
#  while: b[j,n](t) =  ----------------- b[j,n-1](t) + --------------------- b[j+1,n-1](t)
#                        t[j+n] - t[j]                   t[j+n+1] - t[j+1]

#  or: b[j,n](t) = a1 + a2
#  for all: j = 0, ..., m-n-2   # this doesn't look right to me
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
      b <- ((x > K[j]) & (x <= K[j+1]))*1  # this is described in line 59, uses the TRUE or FALSE statement gained by determining t[j] <= t <= t[j+1]
    }
    if(n>0) {
      b <- a1 + a2
    }
    return(b)
  }
# This is where
#   x = a time series (like 0 to 1 in 0.01 increments)
#   j = a sequence traditionally from 0, ..., n, however this is represented as 1, ..., n+1 as you can't get the zero index of a vector
#   n = degree of spline used (3 = cubic etc.)
#   K = augmented knot set
#   K = interior knots + boundary knots * if(lower.boundary) n else if (upper.boundary) n+1
#       it appears the first knot in the augmented knot set is missing for some
#       reason, as in equations typically there are n+1 boundary knots due to
#       the recursive nature of splines

#  So, for instance, for splines of degree 1, we have
  u <- seq(0, 1, by = 0.01)
  plot(u, Bfun(u, 1, 1, c(0, 0.4, 1, 1)),
    lwd = 2, col = "red", type = "l", ylim = c(0, 1))
  lines(u, Bfun(u, 2, 1, c(0, 0.4, 1, 1)), lwd = 2, col = "blue")
  abline(v = c(0, 0.4, 1), lty = 2)

# So here we have one internal knot 0.4, with our boundary knots 0 and 1, with an appended knot (none on the lower boundary due to degree of 1)

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
  plot(cars)
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[2]*u0 + reg$coefficients[1]
  x1 <- seq(min(cars$speed), K, length = length(u0))
  lines(x1, v, col = "green", lwd = 2)
  u0 <- seq(0, 1, by = 0.01)
  v <- (reg$coefficients[3] - reg$coefficients[2])*u0 +
    reg$coefficients[1] +
    reg$coefficients[2]
  x2 <- seq(K, max(cars$speed), length = length(u0))
  lines(x2, v, col = "blue", lwd = 2)
#	which gives us exactly the graph we obtained previously.
# But we can also consider
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 1), data = cars)
  k <- (K - min(cars$speed))/(max(cars$speed) - min(cars$speed))
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[1] +
    reg$coefficients[2]*Bfun(u0, 1, 1, c(0, k, 1, 1)) +
    reg$coefficients[3]*Bfun(u0, 2, 1, c(0, k, 1, 1))
  lines(x = min(cars$speed) + u0*(max(cars$speed) - min(cars$speed)), y = v,
     col = "purple", lwd = 2)
  abline(v = K, lty = 2, col = "red")

# So, we should be able to try with two knots (but we keep it linear, so far)
  K <- c(14, 20)
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 1), data = cars)
  u <- seq(4, 25, by = 0.1)
  B <- data.frame(speed = u)
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
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 1), data = cars)
  k <- (K - min(cars$speed))/(max(cars$speed) - min(cars$speed))
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[1] +
    reg$coefficients[2]*Bfun(u0, 1, 1, c(0, k, 1, 1)) +
    reg$coefficients[3]*Bfun(u0, 2, 1, c(0, k, 1, 1)) +
    reg$coefficients[4]*Bfun(u0, 3, 1, c(0, k, 1, 1))
  lines(x = min(cars$speed)+u0*(max(cars$speed)-min(cars$speed)), y = v,
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
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 2), data = cars)
  k <- (K - min(cars$speed))/(max(cars$speed) - min(cars$speed))
  u0 <- seq(0, 1, by = 0.01)
  v <- reg$coefficients[1] +
    reg$coefficients[2]*Bfun(u0, 1, 2, c(0, 0, k, 1, 1, 1)) +
    reg$coefficients[3]*Bfun(u0, 2, 2, c(0, 0, k, 1, 1, 1)) +
    reg$coefficients[4]*Bfun(u0, 3, 2, c(0, 0, k, 1, 1, 1)) +
    reg$coefficients[5]*Bfun(u0, 4, 2, c(0, 0, k, 1, 1, 1))
  lines(x = min(cars$speed) + u0*(max(cars$speed) - min(cars$speed)), y = v,
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
    K <- min(cars$speed) + k*(max(cars$speed) - min(cars$speed))
    reg <- lm(dist ~ bs(speed, knots = c(K), degree = 2), data = cars)
    SSR[i] <- sum(residuals(reg)^2)
  }
  plot(vk, SSR, type = "b", col = "blue")

# Here, the best model is obtained when we split 3/4-1/4…
