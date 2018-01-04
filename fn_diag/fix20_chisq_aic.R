# Incorporate sigma into model fitting
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
            "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
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

# Load packages
  library(GA)

# Source scripts to set up environment
  set.seed(256256)
  niter <- 100
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "study_rdata.R", sep = "/"))

# -----------------------------------------------------------------------------
  chisq.sumexp.aic <- function(opt) {
    x <- unlist(opt$value)
    k <- unlist(lapply(opt$par, length))
    aic <- x + 2*k
    return(sapply(opt, function(x) x[which(aic == min(aic))]))
  }

  chisq.sumexp.bic <- function(opt, nobs) {
    x <- unlist(opt$value)
    k <- unlist(lapply(opt$par, length))
    aic <- x + log(nobs)*k
    return(sapply(opt, function(x) x[which(aic == min(aic))]))
  }

# -----------------------------------------------------------------------------
# Set up
  data <- d2b
  par <- d2b.p
  t1 <- d2b.t

  niter <- dim(data)[2]
  absorp <- ifelse((dim(par)[1] %% 2) != 0, T, F)
  err <- matrix(
    1 + rnorm(n = length(t1)*niter, mean = 0, sd = 0.05),
    nrow = length(t1), ncol = niter
  )
  subd <- data[which(time.samp %in% t1),]*err
  res.sumexp <- apply(subd, 2, function(x) {
    optim.sumexp.hes(
      data.frame(time = t1, conc = x), oral = absorp
    )
  })
  auc.tlast.true <- apply(par, 2, function(x) pred.tlast.lam(x))

# Original
  lrt.sumexp <- lapply(res.sumexp, chisq.sumexp)
  fit.par <- lapply(lrt.sumexp, function(x) x$par)
  auc.tlast.lrt <- unlist(lapply(fit.par, function(x) pred.tlast.lam(x)))

# AIC
  aic.sumexp <- lapply(res.sumexp, chisq.sumexp.aic)
  fit.par <- lapply(aic.sumexp, function(x) x$par)
  auc.tlast.aic <- unlist(lapply(fit.par, function(x) pred.tlast.lam(x)))

# BIC
  bic.sumexp <- lapply(res.sumexp, chisq.sumexp.bic, nobs = length(t1))
  fit.par <- lapply(bic.sumexp, function(x) x$par)
  auc.tlast.bic <- unlist(lapply(fit.par, function(x) pred.tlast.lam(x)))
