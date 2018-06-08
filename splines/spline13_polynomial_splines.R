# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  K <- c(14)  # knots
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 1), data = cars)
  u <- seq(4, 25, by = 0.1)  # test model at these speeds
  B <- data.frame(speed = u)  # data.frame(speeds)
  Y <- predict(reg, newdata = B)  # predicted dist for the desired speeds
  lines(u, Y, lwd = 2, col = "red")

#  i.e. we have the following (nice) picture
#  But if we look at the output of the regression we get this
  summary(reg)

# Well we can replicate this
  plot(cars)
  u0 <- seq(0, 1, by = 0.01)  # splines run from 0 -> 1 (with the resolution of the curve dictated by how many points you ask for)
  v1 <- reg$coefficients[2]*u0 + reg$coefficients[1]  # here we calculate y = mx + c m = coefficient 2, c = coefficient 1
  x1 <- seq(min(cars$speed), K[1], length = length(u0))  # create times for plotting from first time: min(time), to the first knot: K[1]
  lines(x1, v1, col = "green", lwd = 2)  # draw the line
  u0 <- seq(0, 1, by = 0.01)  # set up the second section of spline
  v2 <- (reg$coefficients[3] - reg$coefficients[2])*u0 +
    reg$coefficients[1] +
    reg$coefficients[2]  # here we calculate m (as coefficient 3 - coefficient 2) and c (coefficient 1 + coefficient 2)
  x2 <- seq(K[1], max(cars$speed), length = length(u0))  # create times for plotting from the first knot: K[1], to the final time  first time: max(time)
  lines(x2, v2, col = "blue", lwd = 2)  # draw the second line

# But what about when things get more complex? Let's say with two knots!
  K <- c(14, 20)  # knots
  plot(cars)
  reg <- lm(dist ~ bs(speed, knots = c(K), degree = 1), data = cars)
  u <- seq(4, 25, by = 0.1)  # test model at these speeds
  B <- data.frame(speed = u)  # data.frame(speeds)
  Y <- predict(reg, newdata = B)  # predicted dist for the desired speeds
  lines(u, Y, lwd = 2, col = "red")

# This can be replicated as well!
  plot(cars)
### this part is all the same as with the one knot
  u0 <- seq(0, 1, by = 0.01)
  v1 <- reg$coefficients[2]*u0 + reg$coefficients[1]
  x1 <- seq(min(cars$speed), K[1], length = length(u0))
  lines(x1, v1, col = "green", lwd = 2)
  v2 <- (reg$coefficients[3] - reg$coefficients[2])*u0 +
    reg$coefficients[1] +
    reg$coefficients[2]
  x2 <- seq(K[1], K[2], length = length(u0))  # except here we calculate times from the first knot to the second knot!
  lines(x2, v2, col = "blue", lwd = 2)
### this part is all the same as with the one knot
# then we need to calculate for the extra bit of linear regression we did!
  v3 <- (reg$coefficients[4] - reg$coefficients[3])*u0 +
    reg$coefficients[1] +
    reg$coefficients[3]
  x3 <- seq(K[2], max(cars$speed), length = length(u0))
  lines(x3, v3, col = "purple", lwd = 2)

# Therefore the output of summary(model) is:
# coefficient 1 <- intercept
# coefficient 2 <- slope 1
# coefficient 3 <- slope 2
# ... coefficient x <- slope x

# For the first linear spline:
# m = coefficient 2
# c = coefficient 1
# For the nth linear spline:
# m = coefficient n+1 - coefficient n
# c = coefficient 1 + coefficient n

# Therefore we can use a function to extract this information!
  interpLinearLoop <- function(x, y, K) {
  # function to provide slope and intercept parameters for a linear interpolation
    require(splines)
    mod <- lm(y ~ bs(x, knots = c(K), degree = 1))  # fit linear spline
    p <- unname(mod$coefficients)  # create object for model coefficients
    m <- c(NULL)  # create placeholder vectors to store values in later
    b <- c(NULL)  # use b instead of c to represent intercept
    for (i in 1:(length(K)+1)) {  # for all the required knots
      if (i == 1) {
        m <- c(m, p[2])
        b <- c(b, p[1])
      } else {
        m <- c(m, p[i+1] - p[i])
        b <- c(b, p[1] + p[i])
      }
    }
    data.frame(
      K = paste(c(min(x), K), c(K, max(x)), sep = "-"),
      m = m,
      c = b
    )
  }
  interpLinearLoop(cars$speed, cars$dist, c(14, 20))

# Can also be done without for loops! # This is the faster function, so use this one
  interpLinearOld <- function(x, y, K) {
  # function to provide slope and intercept parameters for a linear interpolation
    require(splines)
    mod <- lm(y ~ bs(x, knots = c(K), degree = 1))  # fit linear spline
    p <- unname(mod$coefficients)  # create object for model coefficients
    p1 <- tail(p, length(K))
    p2 <- head(p[-1], length(K))
    data.frame(
      K = paste(c(min(x), K), c(K, max(x)), sep = "-"),
      m = c(p[2], p1 - p2),
      c = c(p[1], p[1] + p2)
    )
  }
  interpLinearOld(cars$speed, cars$dist, c(14, 20))


# how about something that looks like concentrations
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  ERR <- 0.05	# Standard deviation of error
  dose <- 50	# mg
  times <- seq(from = 0,to = 24,by = 0.5)	# Time sequence for simulating concentrations
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))*(1+rnorm(n = length(times),mean = 0,sd = ERR))
  plot(conc ~ times)
  knots <- tail(head(times, -1), -1)  # the knots mustn't include the final or start time!

  out <- interpLinearOld(times, conc, knots)

# These values are not linear interpolation values, but in fact linear splines still
# To convert them to linear slopes more needs to be done!
  y1 <- out$m[1]*times + out$c[1]
  y2 <- out$m[2]*times + out$c[2]
  y3 <- out$m[3]*times + out$c[3]
  plot(conc ~ times)
  lines(y1 ~ times)
  lines(y2 ~ times)
  lines(y3 ~ times)

# The slopes are correct if the difference between each knot is equal to 1
# But in reality this is not true at all, and could be 0.5, or 10 etc.
# So you must divide the spline slopes by the difference between independent variable
# for time as an independent variable this includes 0 and the final time, which we removed when specifying knots
  out$m2 <- out$m/diff(times)

# The intercept is correct if you start at time = 0 each time, but you aren't
# you start at the end of the previous line, so must draw the slope back to zero
# So we subract the slope multplied by how much time has passed from the concentration
  out$c2 <- out$c - head(times, -1)*out$m2

# Now we are good
  fixed1 <- out$m2[1]*times + out$c2[1]
  fixed2 <- out$m2[2]*times + out$c2[2]
  fixed3 <- out$m2[3]*times + out$c2[3]
  plot(conc ~ times)
  lines(fixed1 ~ times)
  lines(fixed2 ~ times)
  lines(fixed3 ~ times)

# Last function
  interpLinearNew <- function(x, y, K) {
  # function to provide slope and intercept parameters for a linear interpolation
    require(splines)
    mod <- lm(y ~ bs(x, knots = c(K), degree = 1))  # fit linear spline
    p <- unname(mod$coefficients)  # create object for model coefficients
    p1 <- tail(p, length(K))
    p2 <- head(p[-1], length(K))
    all.K <- c(head(x, 1), K, tail(x, 1))
    m <- c(p[2], p1 - p2)/diff(all.K)
    data.frame(
      K = paste(c(min(x), K), c(K, max(x)), sep = "-"),
      m = m,
      c = c(p[1], p[1] + p2) - head(all.K, -1)*m
    )
  }
  out2 <- interpLinearNew(cars$speed, cars$dist, c(14, 20))

  fixed1 <- out2$m[1]*cars$speed + out2$c[1]
  fixed2 <- out2$m[2]*cars$speed + out2$c[2]
  fixed3 <- out2$m[3]*cars$speed + out2$c[3]
  with(cars, plot(dist ~ speed))
  with(cars, lines(fixed1 ~ speed, col = "green"))
  with(cars, lines(fixed2 ~ speed, col = "blue"))
  with(cars, lines(fixed3 ~ speed, col = "purple"))
  
