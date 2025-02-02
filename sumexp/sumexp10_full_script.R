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

# Setup directory
  source(paste(git.dir, reponame, "sumexp_functions.R", sep = "/"))
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
# Workflow for functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# First function that runs is
  list.sumexp1 <- optim.sumexp(data1, oral = T, nexp = 4)

# Then we determine the best model using chi-squares
  best.sumexp1 <- chisq.sumexp(list.sumexp1)

# Now we use these optim results along with our time sequence
# Originally 49 samples were taken (so many!), but now we only want to use 12
  nobs <- 12
  tlast <- 24
  time.seq <- c(0, exp(seq(-2.5, 0, length.out = nobs))*tlast)
  interv1 <- optim.interv(time.seq, best.sumexp1$par)

# Lastly we need to provide the information that the output
  output <- list(
    sumexp.par = best.sumexp1$par,
    time.interv = c(0, interv1$par, tlast)
  )
