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
  sdev <- 1
  source(paste(git.dir, reponame, "fn_diag/fix_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "fn_diag/study_rdata_resamp.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set basic parameters
  data.names <- paste0("d", as.vector(outer(1:3, c("b", "a"), paste0)))[-1]
  par.names <- paste(data.names, "p", sep = ".")
  fn.names <- paste("pred", data.names, sep = ".")
  t1.names <- paste(data.names, "t", sep = ".")

  study.fn <- function(data, par, fn, nobs, t1, tlast = 24, logauc = F) {  # sdev = 1:4
    niter <- dim(data)[2]
    absorp <- ifelse((dim(par)[1] %% 2) != 0, T, F)
    if (absorp) data[1, ] <- 0
    all.sumexp <- apply(data, 2, function(x) {
      optim.sumexp.new(
      # optim.sumexp.sig(
        data.frame(time = t1, conc = x), oral = absorp
        # , nexp = 2
      )
    })
    print("sumexp done")
    res.sumexp <- lapply(all.sumexp, best.sumexp.aic)  # ".lrt)", ".aic)", ".bic, nobs = length(t1))"
    fit.par <- lapply(res.sumexp, function(x) x$sumexp)
    true.tlast <- apply(par, 2, function(x) {
      list(seq(0, pred.tlast.lam(x), length.out = length(t1)))
    })
    auc.tlast <- lapply(fit.par, function(x) {
      out <- pred.tlast(x, 12)[1]
      if (out == 0) browser()
      seq(0, out, length.out = length(t1))
    })
    lapply(auc.tlast, function(x) {
      if (tail(x, 1) == 0) 0
      else 1
    })
    lam.tlast <- lapply(fit.par, function(x) {
      seq(0, pred.tlast.lam(x), length.out = length(t1))
    })
    obs.tlast.mat <- apply(data, 2, function(x) {
      out <- try(seq(0, obs.tlast.lam(data.frame(t1, x)), length.out = length(t1)))
      if (class(out) == "try-error") browser()
      out
    })
    obs.tlast <- split(t(obs.tlast.mat), seq(NROW(t(obs.tlast.mat))))
    print("tlast done")
  # Explanation of Option Naming
  #         a b c
  #       t 0 0 0
  # a - fixed tmax (0 - off, 1 - on)
  # b - variable tlast (0 - off, 1 - 80% auc, 2 - three half-lives, 3 - three half-lives using observed data)
  # c - optimal lambdaz (0 - off, 1 - geomean, 2 - optimise)

    t000.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t1)
    })
    t001.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t1[-(nobs-1)])$times
    })
    t010.res <- mapply(fit.par, auc.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t)
    })
    t011.res <- mapply(fit.par, auc.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t[-(nobs-1)])$times
    })
    t020.res <- mapply(fit.par, lam.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t)
    })
    t021.res <- mapply(fit.par, lam.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t[-(nobs-1)])$times
    })
    t030.res <- mapply(fit.par, obs.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t)
    })
    t031.res <- mapply(fit.par, obs.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t[-(nobs-1)])$times
    })
    t100.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t1, tmax = T)
    })
    t101.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t1[-(nobs-1)], tmax = T)$times
    })
    t110.res <- mapply(fit.par, auc.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t, tmax = T)
    })
    t111.res <- mapply(fit.par, auc.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t[-(nobs-1)], tmax = T)$times
    })
    t120.res <- mapply(fit.par, lam.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t, tmax = T)
    })
    t121.res <- mapply(fit.par, lam.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t[-(nobs-1)], tmax = T)$times
    })
    t130.res <- mapply(fit.par, obs.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t, tmax = T)
    })
    t131.res <- mapply(fit.par, obs.tlast, SIMPLIFY = F, FUN = function(x, t) {
      optim.interv.dtmax(x, t[-(nobs-1)], tmax = T)$times
    })
    print("intervals done")
    t000 <- sapply(t000.res, FUN = function(x) {
      x$times
    })
    t001 <- sapply(t001.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t002 <- mapply(t001.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t010 <- sapply(t010.res, FUN = function(x) {
      x$times
    })
    t011 <- sapply(t011.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t012 <- mapply(t011.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t020 <- sapply(t020.res, FUN = function(x) {
      x$times
    })
    t021 <- sapply(t021.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t022 <- mapply(t021.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t030 <- sapply(t030.res, FUN = function(x) {
      x$times
    })
    t031 <- sapply(t031.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t032 <- mapply(t031.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t100 <- sapply(t100.res, FUN = function(x) {
      x$times
    })
    t101 <- sapply(t101.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t102 <- mapply(t101.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t110 <- sapply(t110.res, FUN = function(x) {
      x$times
    })
    t111 <- sapply(t111.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t112 <- mapply(t111.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t120 <- sapply(t120.res, FUN = function(x) {
      x$times
    })
    t121 <- sapply(t121.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t122 <- mapply(t121.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    t130 <- sapply(t130.res, FUN = function(x) {
      x$times
    })
    t131 <- sapply(t131.res, FUN = function(x) {
      geomean <- exp(mean(log(tail(x, 2))))
      c(head(x, nobs-2), geomean, tail(x, 1))
    })
    t132 <- mapply(t131.res, fit.par, FUN = function(x, fit) {
      tail.par <- tail(unique(x), 2)
      optres <- try(optim(
        mean(tail.par),
        function(p, tfirst, tlast, fit) {
          m <- fit[1:ceiling(length(fit)/2)]
          times <- c(tfirst, p, tlast)
          pred <- data.frame(
            time = times,
            dv = pred.sumexp(fit, times)
          )
          lmres <- lm(log(dv) ~ times, pred)$coefficients
          err <- sqrt(diff(c(max(m), lmres[2]))^2)
          return(err)
        },
        method = "L-BFGS-B", hessian = T,
        lower = tail.par[1], upper = tail.par[2],
        tfirst = tail.par[1], tlast = tail.par[2], fit = fit
      ))
      if (class(optres) == "try-error") {
        sort(c(head(x, nobs-2), exp(mean(log(tail.par))), tail(x, 1)))
      } else {
        sort(c(head(x, nobs-2), optres$par, tail(x, 1)))
      }
    })
    print("times done")
    auc24 <- data.frame(
      true = apply(par, 2, function(x) integrate(fn, 0, 24, p = x)$value),
      basic = apply(par, 2, function(x) auc.interv(t1, x, fn)),
      t000 = mapply(data.frame(par), data.frame(t000),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t001 = mapply(data.frame(par), data.frame(t001),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t002 = mapply(data.frame(par), data.frame(t002),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t010 = mapply(data.frame(par), data.frame(t010),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t011 = mapply(data.frame(par), data.frame(t011),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t012 = mapply(data.frame(par), data.frame(t012),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t020 = mapply(data.frame(par), data.frame(t020),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t021 = mapply(data.frame(par), data.frame(t021),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t022 = mapply(data.frame(par), data.frame(t022),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t030 = mapply(data.frame(par), data.frame(t030),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t031 = mapply(data.frame(par), data.frame(t031),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t032 = mapply(data.frame(par), data.frame(t032),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t100 = mapply(data.frame(par), data.frame(t100),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t101 = mapply(data.frame(par), data.frame(t101),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t102 = mapply(data.frame(par), data.frame(t102),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t110 = mapply(data.frame(par), data.frame(t110),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t111 = mapply(data.frame(par), data.frame(t111),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t112 = mapply(data.frame(par), data.frame(t112),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t120 = mapply(data.frame(par), data.frame(t120),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t121 = mapply(data.frame(par), data.frame(t121),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t122 = mapply(data.frame(par), data.frame(t122),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t130 = mapply(data.frame(par), data.frame(t130),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t131 = mapply(data.frame(par), data.frame(t131),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t132 = mapply(data.frame(par), data.frame(t132),
        FUN = function(x, y) auc.interv(y, x, fn)
      )
    )

    auctlast <- data.frame(
      true = apply(par, 2, function(x) integrate(fn, 0, 168, p = x)$value),
      basic = apply(par, 2, function(x) auc.interv(t1, x, fn)),
      t000 = mapply(data.frame(par), data.frame(t000),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t001 = mapply(data.frame(par), data.frame(t001),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t002 = mapply(data.frame(par), data.frame(t002),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t010 = mapply(data.frame(par), data.frame(t010),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t011 = mapply(data.frame(par), data.frame(t011),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t012 = mapply(data.frame(par), data.frame(t012),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t020 = mapply(data.frame(par), data.frame(t020),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t021 = mapply(data.frame(par), data.frame(t021),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t022 = mapply(data.frame(par), data.frame(t022),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t030 = mapply(data.frame(par), data.frame(t030),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t031 = mapply(data.frame(par), data.frame(t031),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t032 = mapply(data.frame(par), data.frame(t032),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t100 = mapply(data.frame(par), data.frame(t100),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t101 = mapply(data.frame(par), data.frame(t101),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t102 = mapply(data.frame(par), data.frame(t102),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t110 = mapply(data.frame(par), data.frame(t110),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t111 = mapply(data.frame(par), data.frame(t111),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t112 = mapply(data.frame(par), data.frame(t112),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t120 = mapply(data.frame(par), data.frame(t120),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t121 = mapply(data.frame(par), data.frame(t121),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t122 = mapply(data.frame(par), data.frame(t122),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t130 = mapply(data.frame(par), data.frame(t130),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t131 = mapply(data.frame(par), data.frame(t131),
        FUN = function(x, y) auc.interv(y, x, fn)
      ),
      t132 = mapply(data.frame(par), data.frame(t132),
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
      t000 = mapply(data.frame(par), data.frame(t000), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t001 = mapply(data.frame(par), data.frame(t001), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t002 = mapply(data.frame(par), data.frame(t002), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t010 = mapply(data.frame(par), data.frame(t010), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t011 = mapply(data.frame(par), data.frame(t011), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t012 = mapply(data.frame(par), data.frame(t012), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t020 = mapply(data.frame(par), data.frame(t020), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t021 = mapply(data.frame(par), data.frame(t021), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t022 = mapply(data.frame(par), data.frame(t022), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t030 = mapply(data.frame(par), data.frame(t030), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t031 = mapply(data.frame(par), data.frame(t031), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t032 = mapply(data.frame(par), data.frame(t032), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t100 = mapply(data.frame(par), data.frame(t100), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t101 = mapply(data.frame(par), data.frame(t101), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t102 = mapply(data.frame(par), data.frame(t102), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t110 = mapply(data.frame(par), data.frame(t110), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t111 = mapply(data.frame(par), data.frame(t111), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t112 = mapply(data.frame(par), data.frame(t112), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t120 = mapply(data.frame(par), data.frame(t120), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t121 = mapply(data.frame(par), data.frame(t121), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t122 = mapply(data.frame(par), data.frame(t122), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t130 = mapply(data.frame(par), data.frame(t130), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t131 = mapply(data.frame(par), data.frame(t131), FUN = function(x, t) {
        auc <- auc.interv(t, x, fn)
        inf <- auc.interv.lam(x, t)
        auc + inf
      }),
      t132 = mapply(data.frame(par), data.frame(t132), FUN = function(x, t) {
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
      t000 = mapply(data.frame(par), data.frame(t000),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t001 = mapply(data.frame(par), data.frame(t001),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t002 = mapply(data.frame(par), data.frame(t002),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t010 = mapply(data.frame(par), data.frame(t010),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t011 = mapply(data.frame(par), data.frame(t011),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t012 = mapply(data.frame(par), data.frame(t012),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t020 = mapply(data.frame(par), data.frame(t020),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t021 = mapply(data.frame(par), data.frame(t021),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t022 = mapply(data.frame(par), data.frame(t022),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t030 = mapply(data.frame(par), data.frame(t030),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t031 = mapply(data.frame(par), data.frame(t031),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t032 = mapply(data.frame(par), data.frame(t032),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t100 = mapply(data.frame(par), data.frame(t100),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t101 = mapply(data.frame(par), data.frame(t101),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t102 = mapply(data.frame(par), data.frame(t102),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t110 = mapply(data.frame(par), data.frame(t110),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t111 = mapply(data.frame(par), data.frame(t111),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t112 = mapply(data.frame(par), data.frame(t112),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t120 = mapply(data.frame(par), data.frame(t120),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t121 = mapply(data.frame(par), data.frame(t121),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t122 = mapply(data.frame(par), data.frame(t122),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t130 = mapply(data.frame(par), data.frame(t130),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t131 = mapply(data.frame(par), data.frame(t131),
        FUN = function(x, t) max(pred.sumexp(x, t))
      ),
      t132 = mapply(data.frame(par), data.frame(t132),
        FUN = function(x, t) max(pred.sumexp(x, t))
      )
    )

    tmax <- data.frame(
      true = apply(par, 2, tmax.sumexp),
      basic = mapply(data.frame(par), cmax$basic,
        FUN = function(x, cmax) t1[which(pred.sumexp(x, t1) == cmax)][1]
      ),
      t000 = mapply(data.frame(par), cmax$t000, data.frame(t000),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t001 = mapply(data.frame(par), cmax$t001, data.frame(t001),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t002 = mapply(data.frame(par), cmax$t002, data.frame(t002),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t010 = mapply(data.frame(par), cmax$t010, data.frame(t010),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t011 = mapply(data.frame(par), cmax$t011, data.frame(t011),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t012 = mapply(data.frame(par), cmax$t012, data.frame(t012),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t020 = mapply(data.frame(par), cmax$t020, data.frame(t020),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t021 = mapply(data.frame(par), cmax$t021, data.frame(t021),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t022 = mapply(data.frame(par), cmax$t022, data.frame(t022),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t030 = mapply(data.frame(par), cmax$t030, data.frame(t030),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t031 = mapply(data.frame(par), cmax$t031, data.frame(t031),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t032 = mapply(data.frame(par), cmax$t032, data.frame(t032),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t100 = mapply(data.frame(par), cmax$t100, data.frame(t100),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t101 = mapply(data.frame(par), cmax$t101, data.frame(t101),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t102 = mapply(data.frame(par), cmax$t102, data.frame(t102),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t110 = mapply(data.frame(par), cmax$t110, data.frame(t110),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t111 = mapply(data.frame(par), cmax$t111, data.frame(t111),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t112 = mapply(data.frame(par), cmax$t112, data.frame(t112),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t120 = mapply(data.frame(par), cmax$t120, data.frame(t120),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t121 = mapply(data.frame(par), cmax$t121, data.frame(t121),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t122 = mapply(data.frame(par), cmax$t122, data.frame(t122),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t130 = mapply(data.frame(par), cmax$t130, data.frame(t130),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t131 = mapply(data.frame(par), cmax$t131, data.frame(t131),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      ),
      t132 = mapply(data.frame(par), cmax$t132, data.frame(t132),
        FUN = function(x, cmax, t) t[which(pred.sumexp(x, t) == cmax)][1]
      )
    )
    tlast <- data.frame(
      true = sapply(true.tlast, function(x) tail(unlist(x), 1)),
      auc = sapply(auc.tlast, function(x) tail(unlist(x), 1)),
      lam = sapply(lam.tlast, function(x) tail(unlist(x), 1)),
      obs = sapply(obs.tlast, function(x) tail(unlist(x), 1))
    )
    return(list(par = par, fit.par = fit.par, tlast = tlast, sumexp = res.sumexp, tbas = t1,
      t000 = t000, t001 = t001, t002 = t002, t010 = t010, t011 = t011, t012 = t012,
      t020 = t020, t021 = t021, t022 = t022, t030 = t030, t031 = t031, t032 = t032,
      t100 = t100, t101 = t101, t102 = t102, t110 = t110, t111 = t111, t112 = t112,
      t120 = t120, t121 = t121, t122 = t122, t130 = t130, t131 = t131, t132 = t132,
      interv.t000 = t000.res, interv.t001 = t001.res, interv.t010 = t010.res,
      interv.t011 = t011.res, interv.t020 = t020.res, interv.t021 = t021.res,
      interv.t030 = t030.res, interv.t031 = t031.res, interv.t100 = t100.res,
      interv.t101 = t101.res, interv.t110 = t110.res, interv.t111 = t111.res,
      interv.t120 = t120.res, interv.t121 = t121.res, interv.t130 = t130.res,
      interv.t131 = t131.res,
      auc24 = auc24, auctlast = auctlast, aucinf = aucinf, cmax = cmax, tmax = tmax
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
  saveRDS(fin.res[[1]]$result, "d2b-gaiter100-AR3021.rds")
  saveRDS(fin.res[[2]]$result, "d3b-gaiter100-AR3021.rds")
  saveRDS(fin.res[[3]]$result, "d1a-gaiter100-AR3021.rds")
  saveRDS(fin.res[[4]]$result, "d2a-gaiter100-AR3021.rds")
  saveRDS(fin.res[[5]]$result, "d3a-gaiter100-AR3021.rds")

  # source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  # fin.res <- list(NULL)
  # for (i in 5) {
  #   fin.res[[i]] <- list(
  #     data = data.names[i],
  #     result = study.fn(get(data.names[i]),
  #       par = get(par.names[i]), fn = get(fn.names[i]),
  #       t1 = get(t1.names[i]), nobs = 9
  #     )  # study.fn
  #   )  # list
  #   print(paste0(i, "done"))
  # }  # for loop
  # # saveRDS(fin.res[[4]]$result, "d2a-broad-AR2021.rds")
  # saveRDS(fin.res[[5]]$result, "d3a-broad-AR3022.rds")
