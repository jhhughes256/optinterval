# Splines - Opening the Black Box
# by Arthur Charpentier (2010)
# URL: https://www.r-bloggers.com/splines-opening-the-black-box/

# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	set.seed(123)

#  Splines in regression is something which looks like a black box (or
#  maybe like some dishes you get when you travel away from home: it tastes
#  good, but you don’t what’s inside… even if you might have some clues,
#  you never know for sure*). With
#  splines, it is the same: there are knots, then we consider polynomial
#  interpolations on parts between knots, and we make sure that there is
#  no discontinuity (on the prediction, but on the derivative as well).

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
  B <- function(x, j, n, K) {
    b <- 0
    a1 <- 0
    a2 <- 0
    if(((K[j+n+1] > K[j+1]) & (j+n<=length(K)) & (n > 0)) == TRUE) {
      a2 <- (K[j+n+1] - x)/(K[j+n+1] - K[j+1]) * B(x, j+1, n-1, K)
    }
    if(((K[j+n] > K[j]) & (n > 0)) == TRUE) {
      a1 <- (x - K[j])/(K[j+n] - K[j]) * B(x, j, n-1, K)
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
  plot(u, B(u, 1, 1, c(0, 0.4, 1, 1)), lwd = 2, col = "red", type = "l", ylim = c(0, 1))
  lines(u, B(u, 2, 1, c(0, 0.4, 1, 1)), lwd = 2, col = "blue")
  abline(v = c(0, 0.4, 1), lty = 2)

#  and for splines of degree 2, the basis is
	u <- seq(0, 1, by = 0.01)
	plot(u, B(u, 1, 2, c(0 ,0, 0.4, 1, 1, 1)), lwd = 2, col = "red", type = "l", ylim = c(0, 1))
	lines(u, B(u, 2, 2, c(0, 0, 0.4, 1, 1, 1)), lwd = 2, col = "blue")
	lines(u, B(u, 3, 2, c(0, 0, 0.4, 1, 1, 1)), lwd = 2, col = "green")
	abline(v = c(0, 0.4, 1), lty = 2)

#  …etc. Note that I need to duplicate sometimes starting and end points (but it
#  should be possible to fix it in the function).
#  So, how do we use that, here ? Actually, there are two steps:
