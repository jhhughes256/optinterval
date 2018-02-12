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

# -----------------------------------------------------------------------------
  pred.lambdaz <- function(dv, t) {
    if (t[1] == 0) dv[1] <- 0
    i <- 2
    j <- 0
    bestR2 <- -1
    bestk <- 0
    mdv <- which(dv <= 0)
    rem <- matrix(mdv)
    terminal <- which(dv == max(dv)):length(dv)
    if (length(terminal) >= i) {
      repeat {
        for (l in 1:ncol(rem)) {
          mod <- suppressWarnings(lm(log(tail(dv[-rem[,l]], i)) ~ tail(unique(t)[-rem[,l]], i)))
          k <- unname(mod$coefficients)
          R2 <- suppressWarnings(as.numeric(summary(mod)["adj.r.squared"]))
          if (is.na(k[2])) browser()
          if (is.nan(R2)) R2 <- suppressWarnings(as.numeric(summary(mod)["r.squared"]))
          if (k[2] < 0) {
            if (R2 > bestR2) {
              if (i > 2) bestR2 <- R2
              bestk <- k
            } else {
              break  # break out of for loop
            }
          }
        }
        if (R2 < bestR2) break  # break out of repeat loop
        if (length(terminal) == j + 3) {  # last ditch effort (intended for simulation study only)
          bestk <- c(max(dv), log(2)/56)
          break  # break out of repeat loop
        }
        if (i == length(terminal) & bestk[2] == 0) {  # if all models are invalid, trial removing random points
          j <- j + 1
          i <- j + 2
          comb <- length(dv) - combn(i + 1, j) + 1
          rem <- rbind(matrix(rep(mdv, ncol(comb)), ncol = ncol(comb)), comb)
        }
        i <- i + 1
      }
    }
    bestk
  }

  obs.tlast.lam <- function(obs) {
    lambz <- -pred.lambdaz(obs[, 2], obs[, 1])[2]
    hl <- log(2)/lambz
    ceiling(hl/4)*12
  }

  auc.interv.lam <- function(par, t) {
    dv <- pred.sumexp(par, unique(t))
    lambz <- -pred.lambdaz(dv, t)[2]
    tail(dv, 1)/lambz
  }

  par <- c(-0.06, -0.4, -0.6, log(exp(4)*0.15), log(exp(4)*0.85))
  tlast <- 24
  i <- round(tlast/12, 0)
  obstimes <- seq(0, i*12, by = i*12/12)
  obsdata <- data.frame(
    times = obstimes,
    dv = pred.sumexp(par, obstimes)
  )
  optres <- chisq.sumexp(optim.sumexp(obsdata, oral = T))
  pred.tlast(optres$par, 12)
  pred.tlast.lam(optres$par)
  obs.tlast.lam(obsdata)
