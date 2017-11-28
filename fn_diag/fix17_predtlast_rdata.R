# Trying to determine away to predict appropriate tlast
# Using the rule that AUC0-tlast must be 80% of AUC0-inf
# -----------------------------------------------------------------------------
# Set up directory
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "optinterval"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
      reponame <- "optinterval"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "splines"
    }
    rm("wd")
  }

# Load libraries
  library(GA)

# Source functions
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "study_rdata.R", sep = "/"))

# -----------------------------------------------------------------------------
  dataset <- "d3a"
  par <- get(paste0(dataset, ".p"))
  obsdata <- get(dataset)
  optres <- apply(obsdata, 2, function(x) {
    chisq.sumexp(optim.sumexp.hes(
      data.frame(time = time.samp, conc = x), oral = T
    ))$par
  })

  tlast <- lapply(optres, function(x) {
    pred.tlast(x, 24)
  })

  pred.tlast <- function(fit.par, tlast) {
    i <- round(tlast/12, 0)
    perc.term <- 1
    timeloop <- seq(0, i*12, by = i*12/120)
    predloop <- pred.sumexp(fit.par, timeloop)
    tmax <- timeloop[which(predloop == max(predloop))]
    while(perc.term > 0.2) {
      if (exists("init")) {
        repeat {
          i <- i + 1
          timesloop <- seq(0, i*12, by = i*12/120)
          predloop <- pred.sumexp(fit.par, timesloop)
          clast <- tail(predloop, 1)
          if (clast < cterm) break
        }
      }
      clast <- tail(predloop, 1)
      auclast <- auc.interv.sumexp(timeloop, fit.par, log = T)
      lambz <- max(head(fit.par, ceiling(length(fit.par)/2)))
      aucinf <- clast/-lambz
      perc.term <- aucinf/(auclast+aucinf)
      cterm <- clast*(0.18/perc.term)
      init <- 1
    }
    return(c(i*12, 1-perc.term))
  }

  pred.tlast(optres$par, 24)



  # i == 5

  auc060 <- auc.interv.sumexp(timesloop, fit.par, log = T)
  clast60 <- tail(predloop, 1)
  auc60inf <- clast60/-lambz
  auc60inf/(auc060+auc60inf)
  # less than 10% we have gone too far

  times24 <- seq(0, 2*tlast, by = 2*tlast/120)
  pred24 <- pred.sumexp(fit.par, times24)
  whichtmax <- which(data24$dv == max(data24$dv))
  any(data24$dv[whichtmax:121] < newclast)

  times36 <- seq(0, 36, by = 0.3)
  pred24 <- pred.sumexp(fit.par, times24)
  whichtmax <- which(data24$dv == max(data24$dv))
  any(data24$dv[whichtmax:121] < newclast)
