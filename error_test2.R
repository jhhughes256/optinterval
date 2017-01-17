#Some notes on achieving constant numerical integration error

  time <- c(0,0.25,0.5,1,2,4,8,12,24)

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
  k <- 0.083
  testdf$deltat <- (1/(k*abs(testdf$secondd)))^(1/3)

  nobs <- length(time)-1
  ldply(seq_len(nobs), function(i) {
    t <- time[c(i,i+1)]
    abs.err <- -(t[2] - t[1])^3*testdf$secondd[i+1]/12
    c(a = t[1],
      b = t[2],
      abs = abs(abs.err)
    )
  })

  testdf

  #deltat means that for for the preceding time-interval and this error, you could have samples every deltat min?

  #needs to be framed as: given x pk samples what interval apart should they be?  Some sort of iterative process that minimises the intergration error?  Needs a spline function to estimate derivatives at any time based on observed data?
