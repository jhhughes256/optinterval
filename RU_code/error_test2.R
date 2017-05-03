# Some notes on achieving constant numerical integration error

  time <- seq(0,120,by=10)

  comp1abs <- function(time)
  {
    Dose <- 100
    V <- 2
    ka <- 0.05
    kel <- 0.02
    Cp <- (Dose/V)*ka/(ka-kel)*(exp(-kel*time)-exp(-ka*time))
  }

  #Calculate a concentration - time course (representing observed data)
  testdf <- data.frame(time,Cp=comp1abs(time))
  plot(Cp~time, data=testdf)

  #Calculate the first derivative (slope) - deltay/deltax
  testdf$firstd <- c(0,diff(testdf$Cp))/c(0,diff(testdf$time))
  testdf$firstd[1] <- 0

  #Calculate the second derivative (slope of slope)
  testdf$secondd <- c(0,diff(testdf$firstd))/c(0,diff(testdf$time))
  testdf$secondd[1] <- 0

  #As direction of gradient isn't important, use abs()
  plot(abs(secondd)~time, data=testdf)

  #delta t to get constant integration error
  k <- 0.002
  testdf$deltat <- (1/(k*abs(testdf$secondd)))^(1/3)

  testdf

# deltat means that for for the preceding time-interval and this error
# you could have samples every deltat min?

# needs to be framed as: given x pk samples what interval apart should they be?
# Some sort of iterative process that minimises the intergration error?
# Need spline function to estimate derivatives at times based on observed data?
