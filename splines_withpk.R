# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

#  So, consider the following dataset, with the following spline regression,
  library(splines)
  library(ggplot2)

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

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
  reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
  u <- seq(0, 24, by = 0.1)  # test model at these times
  B <- data.frame(time = u)  # data.frame(times)
  Y <- predict(reg, newdata = B)  # predicted conc for the desired times
  plotobj1 <- ggplot()
  plotobj1 <- plotobj1 + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj1 <- plotobj1 + geom_line(aes(x = u, y = Y), size = 1, colour = "red")
  plotobj1

  vk <- seq(0.05, 0.95, by = 0.05)
  SSR <- matrix(NA, length(vk))
  for(i in 1:(length(vk))){
    k <- vk[i]
    K <- min(onecomp$time) + k*(max(onecomp$time) - min(onecomp$time))
    reg <- lm(conc ~ bs(time, knots = c(K), degree = 2), data = onecomp)
    SSR[i] <- sum(residuals(reg)^2)
  }
  dat <- data.frame(vk, c(SSR))
  plotobj2 <- ggplot(dat, aes(x = vk, y = SSR))
  plotobj2 <- plotobj2 + geom_point(shape = 21, size = 2, colour = "blue")
  plotobj2 <- plotobj2 + geom_line(linetype = 2, colour = "blue")
  plotobj2
