# Create a function that is able to parse a matrix containing alot of data
# -----------------------------------------------------------------------------
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
  # library(splines)
  #library(ggplot2)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

# Source scripts to set up environment
  set.seed(256256)
  niter <- 1000
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "study_rdata.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set basic parameters

  data.names <- paste0("d", as.vector(outer(1:3, c("b", "a"), paste0)))[-1]
  par.names <- paste(data.names, "p", sep = ".")
  fn.names <- paste("pred", data.names, sep = ".")
  t1.names <- paste(data.names, "t", sep = ".")

  study.fn <- function(data, par, fn, nobs, t1, sdev = 0.05, tlast = 24, logauc = F) {
    niter <- dim(data)[2]
    absorp <- ifelse((dim(par)[1] %% 2) != 0, T, F)
    err <- matrix(
      1 + rnorm(n = length(t1)*niter, mean = 0, sd = sdev),
      nrow = length(t1), ncol = niter
    )
    subd <- data[which(time.samp %in% t1),]*err
    res.sumexp <- apply(subd, 2, function(x) {
      chisq.sumexp(optim.sumexp.hes(
        data.frame(time = t1, conc = x), oral = absorp
      ))
    })
    fit.par <- lapply(res.sumexp, FUN = function(x) {
      x$par
    })
    res.interv <- lapply(fit.par,
      FUN = function(x) optim.interv.dtmax(x, t1)
    )
    t2 <- sapply(res.interv, FUN = function(x) {
      c(0, round(x$times, 2), tlast)
    })
    res.interv.tmax <- lapply(fit.par,
      FUN = function(x) optim.interv.dtmax(x, t1, tmax = T)
    )
    t3 <- sapply(res.interv.tmax, FUN = function(x) {
      c(0, round(x$times, 2), tlast)
    })
    auc <- data.frame(
      true = apply(par, 2, function(x) integrate(fn, 0, tlast, p = x)$value),
      basic = apply(par, 2, function(x) auc.interv(t1, x, fn)),
      optint = mapply(data.frame(par), data.frame(t2),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      optintwCmax = mapply(data.frame(par), data.frame(t3),
        FUN = function(x, y) auc.interv(y, x, fn)
      )
    )
    cmax <- data.frame(
      true = apply(par, 2, function(x) pred.sumexp(x, tmax.sumexp(x))),
      basic = apply(par, 2, function(x) max(pred.sumexp(x, t1))),
      optint = mapply(data.frame(par), data.frame(t2),
        FUN = function(x, y) max(pred.sumexp(x, y))
      ),
      optintwCmax = mapply(data.frame(par), data.frame(t3),
        FUN = function(x, y) max(pred.sumexp(x, y))
      )
    )
    tmax <- data.frame(
      true = apply(par, 2, tmax.sumexp),
      basic = mapply(data.frame(par), cmax$basic,
        FUN = function(x, y) t1[which(pred.sumexp(x, t1) == y)][1]
      ),
      optint = mapply(data.frame(par), cmax$optint, data.frame(t2),
        FUN = function(x, y, z) z[which(pred.sumexp(x, z) == y)][1]
      ),
      optintwCmax = mapply(data.frame(par), cmax$optintwCmax, data.frame(t3),
        FUN = function(x, y, z) z[which(pred.sumexp(x, z) == y)][1]
      )
    )
    return(list(par = par, fit.par = fit.par, t2 = t2, t3 = t3, sumexp = res.sumexp, interv = res.interv, interv.tmax = res.interv.tmax, auc = auc, cmax = cmax, tmax = tmax))
  }

# -----------------------------------------------------------------------------

  fin.res <- list(NULL)
  for (i in 1:5) {
    fin.res[[i]] <- list(
      data = data.names[i],
      result = study.fn(get(data.names[i]),
        par = get(par.names[i]), fn = get(fn.names[i]),
        t1 = get(t1.names[i]), nobs = 9
      )  # study.fn
    )  # list
    print(paste0(i, "done"))
  }  # for loop
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(fin.res[[1]]$result, "d2b-narrow-dt05-RMSE.rds")
  saveRDS(fin.res[[2]]$result, "d3b-narrow-dt05-RMSE.rds")
  saveRDS(fin.res[[3]]$result, "d1a-narrow-dt05-RMSE.rds")
  saveRDS(fin.res[[4]]$result, "d2a-narrow-dt05-RMSE.rds")
  saveRDS(fin.res[[5]]$result, "d3a-narrow-dt05-RMSE.rds")
