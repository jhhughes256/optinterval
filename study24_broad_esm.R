# Create a function that is able to parse a matrix containing alot of data
# -----------------------------------------------------------------------------
# Set up directories
  git.dir <- "E:/Hughes/Git"
  reponame <- "splines"  #"optinterval"

# Load packages
  library(GA)
  # library(splines)
  #library(ggplot2)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

# Source scripts to set up environment
  set.seed(256256)
  niter <- 1000
  sdev <- 4
  source(paste(git.dir, reponame, "fn_diag/fix_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "fn_diag/study_rdata_resamp.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set basic parameters
  data.names <- paste0("d", as.vector(outer(1:3, c("b", "a"), paste0)))[-1]
  par.names <- paste(data.names, "p", sep = ".")
  fn.names <- paste("pred", data.names, sep = ".")
  t1.names <- paste(data.names, "t", sep = ".")

  study.fn <- function(data, par, fn, nobs, t1, tlast = 24, logauc = F) {
    niter <- dim(data)[2]
    absorp <- ifelse((dim(par)[1] %% 2) != 0, T, F)
    if (absorp) data[1, ] <- 0
    all.sumexp <- apply(data, 2, function(x) {
      optim.sumexp.new(
        data.frame(time = t1, conc = x), oral = absorp
      )
    })
    print("sumexp done")

    res.sumexp <- lapply(all.sumexp, best.sumexp.aic)
    fit.par <- lapply(res.sumexp, function(x) x$sumexp)
    true.tlast <- apply(par, 2, function(x) {
      list(seq(0, pred.tlast.lam(x), length.out = length(t1)))
    })
    obs.tlast.mat <- apply(data, 2, function(x) {
      out <- try(seq(0, obs.tlast.lam(data.frame(t1, x)), length.out = length(t1)))
      if (class(out) == "try-error") browser()
      out
    })
    obs.tlast <- split(t(obs.tlast.mat), seq(NROW(t(obs.tlast.mat))))
    print("tlast done")

    res.interv <- mapply(fit.par, obs.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t)
    })
    res.times <- sapply(res.interv, FUN = function(x) {
      x$times
    })
    print("intervals done")

    auc24 <- data.frame(
      true = apply(par, 2, function(x) integrate(fn, 0, 24, p = x)$value),
      basic = apply(par, 2, function(x) auc.interv(t1, x, fn)),
      optint = mapply(data.frame(par), data.frame(res.times),
        FUN = function(x, y) auc.interv(y, x, fn)
      )
    )

    auctlast <- data.frame(
      true = apply(par, 2, function(x) integrate(fn, 0, 168, p = x)$value),
      basic = apply(par, 2, function(x) auc.interv(t1, x, fn)),
      optint = mapply(data.frame(par), data.frame(res.times),
        FUN = function(x, y) auc.interv(y, x, fn)
      )
    )

    aucinf <- try(data.frame(
      true = apply(par, 2, function(x) {
        auc <- integrate(fn, 0, 168, p = x)$value
        inf <- fn(168, x)/abs(max(x[ceiling(length(x)/2)]))
        auc + inf
      }),
      basic = apply(par, 2, function(x) {
        auc <- auc.interv(t1, x, fn)
        inf <- auc.interv.lam(x, t1)
        auc + inf
      }),
      optint = mapply(data.frame(par), data.frame(res.times), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      })
    ))
    if (class(aucinf) == "try-error") browser()

    test <- try(apply(par, 2, function(x) pred.sumexp(x, tmax.sumexp(x))))
    if (class(test) == "try-error") browser()
    cmax <- data.frame(
      true = apply(par, 2, function(x) pred.sumexp(x, tmax.sumexp(x))),
      basic = apply(par, 2, function(x) max(pred.sumexp(x, t1))),
      optint = mapply(data.frame(par), data.frame(res.times),
        FUN = function(x, t) max(pred.sumexp(x, t))
      )
    )

    tmax <- data.frame(
      true = apply(par, 2, tmax.sumexp),
      basic = mapply(data.frame(par), cmax$basic,
        FUN = function(x, cmax) t1[which(pred.sumexp(x, t1) == cmax)][1]
      ),
      optint = mapply(data.frame(par), cmax$optint, data.frame(res.times),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      )
    )
    tlast <- data.frame(
      true = sapply(true.tlast, function(x) tail(unlist(x), 1)),
      obs = sapply(obs.tlast, function(x) tail(unlist(x), 1))
    )
    return(list(par = par, fit.par = fit.par, tlast = tlast, sumexp = res.sumexp,
      tbas = t1, optint = res.times, interv.optint = res.interv, auc24 = auc24,
      auctlast = auctlast, aucinf = aucinf, cmax = cmax, tmax = tmax
    ))
  }

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
  # Runs
  # ARabcd
  # a - max number of exponentials (2, 3)
  # b - include sigma in mle (0 - off, 1 - on)
  # c - ofv comparitive criterion (1 - lrt, 2 - aic, 3 - bic)
  # d - error model (1 - loprop, 2- hiprop, 3 - loboth, 4- hiboth)
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(fin.res[[1]]$result, "d2b-broad-AR3024.rds")
  saveRDS(fin.res[[2]]$result, "d3b-broad-AR3024.rds")
  saveRDS(fin.res[[3]]$result, "d1a-broad-AR3024.rds")
  saveRDS(fin.res[[4]]$result, "d2a-broad-AR3024.rds")
  saveRDS(fin.res[[5]]$result, "d3a-broad-AR3024.rds")
