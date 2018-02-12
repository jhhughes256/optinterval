# Function containing the framework setup in the intial study design
# -----------------------------------------------------------------------------
# Set up directories
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

# Load packages
  library(GA)
  #library(ggplot2)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

# Source scripts to set up environment
  set.seed(256256)
  niter <- 1000
  sdev <- 4
  source(paste(git.dir, reponame, "fn_diag/fix_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "fn_diag/study_data.R", sep = "/"))

# Set basic parameters
  data.names <- paste0("d", as.vector(outer(1:3, c("b", "a"), paste0)))[-1]
  par.names <- paste(data.names, "p", sep = ".")
  fn.names <- paste("pred", data.names, sep = ".")
  t0.names <- paste(data.names, "t", sep = ".")
  t1 <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24)

# -----------------------------------------------------------------------------
  study.fn <- function(data, par, fn, nobs, t0, tlast = 24, logauc = F) {  # sdev = 1:4
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    if (absorp) data[1] <- 0
    all.sumexp <- apply(data, 2, function(x) {
      optim.sumexp.new(
      # optim.sumexp.sig(
        data.frame(time = t0, conc = x), oral = absorp
        # , nexp = 2
      )
    })
    print("sumexp done")
    res.sumexp <- lapply(all.sumexp, best.sumexp.aic)  # ".lrt)", ".aic)", ".bic, nobs = length(t1))"
    fit.par <- lapply(res.sumexp, function(x) x$sumexp)
    true.tlast <- rep(list(
      seq(0, pred.tlast.lam(par), length.out = nobs)
    ), niter)
    auc.tlast <- lapply(fit.par, function(x) {
      c(0, exp(seq(log(t0[2]), log(pred.tlast(x, 12)[1]), length.out = nobs-1)))
    })
    lam.tlast <- lapply(fit.par, function(x) {
      c(0, exp(seq(log(t0[2]), log(pred.tlast.lam(x)), length.out = nobs-1)))
    })
    obs.tlast.mat <- apply(data, 2, function(x) {
      out <- try(
        c(0, exp(seq(log(t0[2]), log(obs.tlast.lam(data.frame(t0, x))), length.out = nobs-1)))
      )
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
    t2 <- c(0, exp(seq(log(t0[2]), log(tail(t0, 1)), length.out = nobs-1)))
    t000.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t2)
    })
    t001.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t2[-(nobs-1)])$times
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
      optim.interv.dtmax(x, t2, tmax = T)
    })
    t101.res <- lapply(fit.par, FUN = function(x) {
      optim.interv.dtmax(x, t2[-(nobs-1)], tmax = T)$times
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
      true = rep(integrate(fn, 0, 24, p = par)$value, niter),
      basic = rep(auc.interv(t1, par, fn), niter),
      t000 = apply(t000, 2, function(x) auc.interv(x, par, fn)),
      t001 = apply(t001, 2, function(x) auc.interv(x, par, fn)),
      t002 = apply(t002, 2, function(x) auc.interv(x, par, fn)),
      t010 = apply(t010, 2, function(x) auc.interv(x, par, fn)),
      t011 = apply(t011, 2, function(x) auc.interv(x, par, fn)),
      t012 = apply(t012, 2, function(x) auc.interv(x, par, fn)),
      t020 = apply(t020, 2, function(x) auc.interv(x, par, fn)),
      t021 = apply(t021, 2, function(x) auc.interv(x, par, fn)),
      t022 = apply(t022, 2, function(x) auc.interv(x, par, fn)),
      t030 = apply(t030, 2, function(x) auc.interv(x, par, fn)),
      t031 = apply(t031, 2, function(x) auc.interv(x, par, fn)),
      t032 = apply(t032, 2, function(x) auc.interv(x, par, fn)),
      t100 = apply(t100, 2, function(x) auc.interv(x, par, fn)),
      t101 = apply(t101, 2, function(x) auc.interv(x, par, fn)),
      t102 = apply(t102, 2, function(x) auc.interv(x, par, fn)),
      t110 = apply(t110, 2, function(x) auc.interv(x, par, fn)),
      t111 = apply(t111, 2, function(x) auc.interv(x, par, fn)),
      t112 = apply(t112, 2, function(x) auc.interv(x, par, fn)),
      t120 = apply(t120, 2, function(x) auc.interv(x, par, fn)),
      t121 = apply(t121, 2, function(x) auc.interv(x, par, fn)),
      t122 = apply(t122, 2, function(x) auc.interv(x, par, fn)),
      t130 = apply(t130, 2, function(x) auc.interv(x, par, fn)),
      t131 = apply(t131, 2, function(x) auc.interv(x, par, fn)),
      t132 = apply(t132, 2, function(x) auc.interv(x, par, fn))
    )

    auctlast <- data.frame(
      true = rep(integrate(fn, 0, tail(true.tlast[[1]], 1), p = par)$value, niter),
      basic = rep(auc.interv(t1, par, fn), niter),
      t000 = apply(t000, 2, function(x) auc.interv(x, par, fn)),
      t001 = apply(t001, 2, function(x) auc.interv(x, par, fn)),
      t002 = apply(t002, 2, function(x) auc.interv(x, par, fn)),
      t010 = apply(t010, 2, function(x) auc.interv(x, par, fn)),
      t011 = apply(t011, 2, function(x) auc.interv(x, par, fn)),
      t012 = apply(t012, 2, function(x) auc.interv(x, par, fn)),
      t020 = apply(t020, 2, function(x) auc.interv(x, par, fn)),
      t021 = apply(t021, 2, function(x) auc.interv(x, par, fn)),
      t022 = apply(t022, 2, function(x) auc.interv(x, par, fn)),
      t030 = apply(t030, 2, function(x) auc.interv(x, par, fn)),
      t031 = apply(t031, 2, function(x) auc.interv(x, par, fn)),
      t032 = apply(t032, 2, function(x) auc.interv(x, par, fn)),
      t100 = apply(t100, 2, function(x) auc.interv(x, par, fn)),
      t101 = apply(t101, 2, function(x) auc.interv(x, par, fn)),
      t102 = apply(t102, 2, function(x) auc.interv(x, par, fn)),
      t110 = apply(t110, 2, function(x) auc.interv(x, par, fn)),
      t111 = apply(t111, 2, function(x) auc.interv(x, par, fn)),
      t112 = apply(t112, 2, function(x) auc.interv(x, par, fn)),
      t120 = apply(t120, 2, function(x) auc.interv(x, par, fn)),
      t121 = apply(t121, 2, function(x) auc.interv(x, par, fn)),
      t122 = apply(t122, 2, function(x) auc.interv(x, par, fn)),
      t130 = apply(t130, 2, function(x) auc.interv(x, par, fn)),
      t131 = apply(t131, 2, function(x) auc.interv(x, par, fn)),
      t132 = apply(t132, 2, function(x) auc.interv(x, par, fn))
    )

    aucinf <- try(data.frame(
      true = {
        auc <- integrate(fn, 0, 168, p = par)$value
        inf <- fn(168, par)/abs(max(par[ceiling(length(par)/2)]))
        rep(auc + inf, niter)
      },
      basic = {
        auc <- auc.interv(t1, par, fn)
        inf <- auc.interv.lam(par, t1)
        rep(auc + inf, niter)
      },
      t000 = apply(t000, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t001 = apply(t001, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t002 = apply(t002, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t010 = apply(t010, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t011 = apply(t011, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t012 = apply(t012, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t020 = apply(t020, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t021 = apply(t021, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t022 = apply(t022, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t030 = apply(t030, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t031 = apply(t031, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t032 = apply(t032, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t100 = apply(t100, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t101 = apply(t101, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t102 = apply(t102, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t110 = apply(t110, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t111 = apply(t111, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t112 = apply(t112, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t120 = apply(t120, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t121 = apply(t121, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t122 = apply(t122, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t130 = apply(t130, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t131 = apply(t131, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      }),
      t132 = apply(t132, 2, FUN = function(t) {
        auc <- auc.interv(t, par, fn)
        inf <- auc.interv.lam(par, t)
        auc + inf
      })
    ))
    if (class(aucinf) == "try-error") browser()
    test <- pred.sumexp(par, tmax.sumexp(par))
    if (class(test) == "try-error") browser()
    cmax <- data.frame(
      true = rep(pred.sumexp(par, tmax.sumexp(par)), niter),
      basic = rep(max(pred.sumexp(par, t1)), niter),
      t000 = apply(t000, 2, function(t) max(pred.sumexp(par, t))),
      t001 = apply(t001, 2, function(t) max(pred.sumexp(par, t))),
      t002 = apply(t002, 2, function(t) max(pred.sumexp(par, t))),
      t010 = apply(t010, 2, function(t) max(pred.sumexp(par, t))),
      t011 = apply(t011, 2, function(t) max(pred.sumexp(par, t))),
      t012 = apply(t012, 2, function(t) max(pred.sumexp(par, t))),
      t020 = apply(t020, 2, function(t) max(pred.sumexp(par, t))),
      t021 = apply(t021, 2, function(t) max(pred.sumexp(par, t))),
      t022 = apply(t022, 2, function(t) max(pred.sumexp(par, t))),
      t030 = apply(t030, 2, function(t) max(pred.sumexp(par, t))),
      t031 = apply(t031, 2, function(t) max(pred.sumexp(par, t))),
      t032 = apply(t032, 2, function(t) max(pred.sumexp(par, t))),
      t100 = apply(t100, 2, function(t) max(pred.sumexp(par, t))),
      t101 = apply(t101, 2, function(t) max(pred.sumexp(par, t))),
      t102 = apply(t102, 2, function(t) max(pred.sumexp(par, t))),
      t110 = apply(t110, 2, function(t) max(pred.sumexp(par, t))),
      t111 = apply(t111, 2, function(t) max(pred.sumexp(par, t))),
      t112 = apply(t112, 2, function(t) max(pred.sumexp(par, t))),
      t120 = apply(t120, 2, function(t) max(pred.sumexp(par, t))),
      t121 = apply(t121, 2, function(t) max(pred.sumexp(par, t))),
      t122 = apply(t122, 2, function(t) max(pred.sumexp(par, t))),
      t130 = apply(t130, 2, function(t) max(pred.sumexp(par, t))),
      t131 = apply(t131, 2, function(t) max(pred.sumexp(par, t))),
      t132 = apply(t132, 2, function(t) max(pred.sumexp(par, t)))
    )

    tmax <- data.frame(
      true = rep(tmax.sumexp(par), niter),
      basic = apply(matrix(cmax$basic), 1,
        FUN = function(x) t1[which(pred.sumexp(par, t1) == x)][1]
      ),
      t000 = mapply(cmax$t000, data.frame(t000),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t001 = mapply(cmax$t001, data.frame(t001),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t002 = mapply(cmax$t002, data.frame(t002),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t010 = mapply(cmax$t010, data.frame(t010),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t011 = mapply(cmax$t011, data.frame(t011),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t012 = mapply(cmax$t012, data.frame(t012),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t020 = mapply(cmax$t020, data.frame(t020),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t021 = mapply(cmax$t021, data.frame(t021),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t022 = mapply(cmax$t022, data.frame(t022),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t030 = mapply(cmax$t030, data.frame(t030),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t031 = mapply(cmax$t031, data.frame(t031),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t032 = mapply(cmax$t032, data.frame(t032),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t100 = mapply(cmax$t100, data.frame(t100),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t101 = mapply(cmax$t101, data.frame(t101),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t102 = mapply(cmax$t102, data.frame(t102),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t110 = mapply(cmax$t110, data.frame(t110),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t111 = mapply(cmax$t111, data.frame(t111),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t112 = mapply(cmax$t112, data.frame(t112),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t120 = mapply(cmax$t120, data.frame(t120),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t121 = mapply(cmax$t121, data.frame(t121),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t122 = mapply(cmax$t122, data.frame(t122),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t130 = mapply(cmax$t130, data.frame(t130),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t131 = mapply(cmax$t131, data.frame(t131),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
      ),
      t132 = mapply(cmax$t132, data.frame(t132),
        FUN = function(x, t) t[which(pred.sumexp(par, t) == x)][1]
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
        t0 = get(t0.names[i]), nobs = 6
      )  # study.fn
    )  # list
    print(paste0(i, "done"))
  }  # for loop
  setwd("E:/Hughes/Git/splines/fn_diag")
  saveRDS(fin.res[[1]]$result, "d2b-varnobs6-AR3024.rds")
  saveRDS(fin.res[[2]]$result, "d3b-varnobs6-AR3024.rds")
  saveRDS(fin.res[[3]]$result, "d1a-varnobs6-AR3024.rds")
  saveRDS(fin.res[[4]]$result, "d2a-varnobs6-AR3024.rds")
  saveRDS(fin.res[[5]]$result, "d3a-varnobs6-AR3024.rds")
