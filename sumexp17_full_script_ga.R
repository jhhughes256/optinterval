# Combination of all functions in one process
# -----------------------------------------------------------------------------
# How will the function be used?
# How will this translate into the use of the functions provided?
# An overarching function should be available to run the underlying functions
# Should do so in a reliable way
# -----------------------------------------------------------------------------
# Setup environment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify git directory and remove previous objects if not already done
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    graphics.off()
    #git.dir <- "E:/Hughes/Git"
    #git.dir <- "C:/Users/Jim Hughes/Documents/GitRepos"
    git.dir <- "C:/Users/hugjh001/Documents"
    reponame <- "optinterval"
  }
  library(GA)
  library(ggplot2)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# Setup directory
  source(paste(git.dir, reponame, "sumexp_functions_ga.R", sep = "/"))
  source(paste(git.dir, reponame, "sumexp_data.R", sep = "/"))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create datasets
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data1 <- data.frame(
    time = time.samp,
    conc = onedata.abs$sumexp*err
  )
  err <- 1 + rnorm(n = length(time.samp), mean = 0, sd = 0.3)
  data2 <- data.frame(
    time = time.samp,
    conc = twodata$sumexp*err
  )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Output Function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create a function that presents the data
  res.out <- function(fitpar, interv, data = NULL, plotdata = F) {
    if (plotdata) {
      if (is.null(data)) stop("'data' argument not provided")
      plotdata <- data.frame(
        time = data$time,
        cobs = data$conc,
        pred = pred.sumexp(fitpar, data[, 1])
      )
      ylim <- c(0, 1.1*max(plotdata$cobs))
      xlim <- c(0, max(plotdata$time))
      plotobj <- NULL
      plotobj <- ggplot(data = plotdata)
      plotobj <- plotobj + ggtitle("Predicted and Observed vs. Time")
      plotobj <- plotobj + geom_point(aes(x = time, y = cobs))
      plotobj <- plotobj + geom_line(aes(x = time, y = pred), colour = "red")
      plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
      plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n", lim = ylim)
      plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
      print(plotobj)
    }
    out <- list(
      sumexp.par = unname(best.sumexp1$par),
      time.interv = c(0, interv1$par, tlast)
    )
    return(out)
  }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Workflow for functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# First function that runs is
  list.sumexp1 <- optim.sumexp(data1, oral = T)

# Then we determine the best model using chi-squares
  best.sumexp1 <- chisq.sumexp(list.sumexp1)

# Now we use these optim results along with our time sequence
# Originally 49 samples were taken (so many!), but now we only want to use 12
  nobs <- 12
  tlast <- 48
  time.seq <- c(0, exp(seq(-2.5, 0, length.out = nobs))*tlast)
  interv1 <- optim.interv(time.seq, best.sumexp1$par)

# Lastly we need to provide the information that the output
  res.out(best.sumexp1$par, c(0, interv1$par, tlast), data1, T)

# First function that runs is
  list.sumexp2 <- optim.sumexp(data2, oral = F)

# Then we determine the best model using chi-squares
  best.sumexp2 <- chisq.sumexp(list.sumexp2)

# Now we use these optim results along with our time sequence
# Originally 49 samples were taken (so many!), but now we only want to use 12
  nobs <- 6
  tlast <- 24
  time.seq <- c(0, exp(seq(-2.5, 0, length.out = nobs))*tlast)
  interv2 <- optim.interv(time.seq, best.sumexp2$par)

# Lastly we need to provide the information that the output
  output2 <- list(
    sumexp.par = best.sumexp2$par,
    time.interv = c(0, interv2$par, tlast)
  )
