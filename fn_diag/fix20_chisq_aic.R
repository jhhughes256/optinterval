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
  auc.tlast.true <- apply(par, 2, function(x) pred.tlast.lam(x))

# Original
  res.sumexp <- apply(subd, 2, function(x) {
    out <- optim.sumexp.hes(
      data.frame(time = t1, conc = x), oral = absorp
    )
    chisq.sumexp(out)
  })
  fit.par <- lapply(res.sumexp, FUN = function(x) {
    x$par
  })
  auc.tlast.lrt <- lapply(fit.par, function(x) pred.tlast.lam(x))

# AIC
  res.sumexp <- apply(subd, 2, function(x) {
    out <- optim.sumexp.hes(
      data.frame(time = t1, conc = x), oral = absorp
    )
    chisq.sumexp.aic(out)
  })
  fit.par <- lapply(res.sumexp, FUN = function(x) {
    x$par
  })
  auc.tlast.aic <- lapply(fit.par, function(x) pred.tlast.lam(x))
