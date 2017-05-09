# Combination of all functions in one process
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Simulate data
# Absorption Curve
  time.samp <- seq(0, 48, by = 1)
  absdata <- data.frame(
    time = time.samp,
    line1 = -0.2*time.samp + 4,
    line2 = -0.1*time.samp + 4
  )
  absdata$sumexp <- exp(absdata$line2) - exp(absdata$line1)
  #with(absdata, plot(time, sumexp))

# 2 Compartment Curve
  twodata <- data.frame(
    time = time.samp,
    line1 = -0.5*time.samp + 6,
    line2 = -0.05*time.samp + 5
  )
  twodata$sumexp <- exp(twodata$line1) + exp(twodata$line2)
  #with(twodata, plot(time, log(sumexp)))

# Create datasets
  sub.time <- c(0, 1, 2, 3, 4, 8, 12, 24)

  #err <- 1 + rnorm(n = length(sub.time), mean = 0, sd = 0.2)
  data1 <- data.frame(
    time = time.samp[which(time.samp %in% sub.time)],
    conc = absdata$sumexp[which(time.samp %in% sub.time)]  # *err
  )
  #err <- 1 + rnorm(n = length(sub.time), mean = 0, sd = 0.2)
  data2 <- data.frame(
    time = time.samp[which(time.samp %in% sub.time)],
    conc = twodata$sumexp[which(time.samp %in% sub.time)]  # *err
  )
# -----------------------------------------------------------------------------
# Workflow for functions
  optim.sumexp(data1, oral = T, nexp = 4)
