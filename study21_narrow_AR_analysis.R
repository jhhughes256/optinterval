# Printing of result tables
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
  library(plyr)
  library(ggplot2)
  library(grid)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Source functions
  source(paste(git.dir, reponame, "sumstat_functions.R", sep = "/"))

# Source files of interest
  run.string <- "narrow6-AR3022"
  d2b <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d2b-", run.string, ".rds"), sep = "/"))
  d3b <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d3b-", run.string, ".rds"), sep = "/"))
  d1a <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d1a-", run.string, ".rds"), sep = "/"))
  d2a <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d2a-", run.string, ".rds"), sep = "/"))
  d3a <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d3a-", run.string, ".rds"), sep = "/"))

# -----------------------------------------------------------------------------
# Data structure
  data.names <- c("d2b", "d3b", "d1a", "d2a", "d3a")
  slot.names <- c("auc24", "auctlast", "aucinf", "cmax", "tmax")
  niter <- 1000

  res <- data.frame(NULL)
  plotdata <- data.frame(NULL)
  for (i in 1:length(data.names)) {
    r.out <- data.frame(NULL)
    d.out <- data.frame(NULL)
    ref.par.nexp <- ceiling(length(get(data.names[i])[["par"]])/2)
    ref.par.m <- get(data.names[i])[["par"]][1:ref.par.nexp]
    ref.par.c <- get(data.names[i])[["par"]][-(1:ref.par.nexp)]
    fit.par.nexp <- unlist(lapply(get(data.names[i])[["fit.par"]], function(x) {
      ceiling(length(x)/2)
    }))
    fit.par.m <- lapply(get(data.names[i])[["fit.par"]], function(x) {
      x[1:ceiling(length(x)/2)]
    })
    fit.par.c <- lapply(get(data.names[i])[["fit.par"]], function(x) {
      x[-(1:ceiling(length(x)/2))]
    })
    fit.par.mse <- lapply(get(data.names[i])[["sumexp"]], function(x) {
      vc_mat <- try(solve(x$hessian))
      if(class(vc_mat) != "try-error") {
        se <- sqrt(diag(vc_mat))
        se_percent <- abs(se/x$par*100)
        se_percent[1:ceiling(length(x$par)/2)]
      } else {
        rep(NA, 4)
      }
    })
    fit.par.cse <- lapply(get(data.names[i])[["sumexp"]], function(x) {
      vc_mat <- try(solve(x$hessian))
      if(class(vc_mat) != "try-error") {
        se <- sqrt(diag(vc_mat))
        se_percent <- abs(se/x$par*100)
      } else {
        rep(NA, 4)
      }
    })
    fit.par.tse <- c(
      rep(NA, niter),
      unlist(lapply(get(data.names[i])[["interv.t000"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t001"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t001"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t010"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t011"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t011"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t020"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t021"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t021"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t030"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t031"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t031"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t100"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t101"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t101"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t110"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t111"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t111"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t120"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t121"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t121"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t130"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t131"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })),
      unlist(lapply(get(data.names[i])[["interv.t131"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      }))
    )
    for (j in 1:length(slot.names)) {
      d.in <- get(data.names[i])[[slot.names[j]]]
      d.melt <- data.frame(
        id = rep(1:niter, 25),
        data = data.names[i],
        metric = slot.names[j],
        ref = rep(d.in$true, 25),
        test = c(d.in$basic, d.in$t000, d.in$t001, d.in$t002, d.in$t010,
          d.in$t011, d.in$t012, d.in$t020, d.in$t021, d.in$t022, d.in$t030,
          d.in$t031, d.in$t032, d.in$t100, d.in$t101, d.in$t102, d.in$t110,
          d.in$t111, d.in$t112, d.in$t120, d.in$t121, d.in$t122, d.in$t130,
          d.in$t131, d.in$t132),
        type = c(rep("bas", niter), rep("000", niter), rep("001", niter),
          rep("002", niter), rep("010", niter), rep("011", niter), rep("012", niter),
          rep("020", niter), rep("021", niter), rep("022", niter), rep("030", niter),
          rep("031", niter), rep("032", niter), rep("100", niter), rep("101", niter),
          rep("102", niter), rep("110", niter), rep("111", niter), rep("112", niter),
          rep("120", niter), rep("121", niter), rep("122", niter), rep("130", niter),
          rep("131", niter), rep("132", niter)),
        ref.m1 = rep(ref.par.m[1], 25*niter),
        ref.m2 = rep(ref.par.m[2], 25*niter),
        ref.m3 = if (length(ref.par.m) > 2) {
            rep(ref.par.m[3], 25*niter)
          } else {
            NA
          },
        ref.m4 = if (length(ref.par.m) > 3) {
            rep(ref.par.m[4], 25*niter)
          } else {
            NA
          },
        ref.c1 = rep(ref.par.c[1], 25*niter),
        ref.c2 = if (is.vector(ref.par.c)) {
            NA
          } else if (length(ref.par.c) > 1) {
            rep(ref.par.c[2], 25*niter)
          } else {
            NA
          },
        ref.c3 = if (is.vector(ref.par.c)) {
            NA
          } else if (length(ref.par.c) > 2) {
            rep(ref.par.c[3], 25*niter)
          } else {
            NA
          },
        test.nexp = rep(fit.par.nexp, 25),
        test.m1 = rep(unlist(lapply(fit.par.m, function(x) x[1])), 25),
        test.m2 = rep(unlist(lapply(fit.par.m, function(x) x[2])), 25),
        test.m3 = rep(unlist(lapply(fit.par.m, function(x) x[3])), 25),
        test.m4 = rep(unlist(lapply(fit.par.m, function(x) x[4])), 25),
        test.c1 = rep(unlist(lapply(fit.par.c, function(x) x[1])), 25),
        test.c2 = rep(unlist(lapply(fit.par.c, function(x) x[2])), 25),
        test.c3 = rep(unlist(lapply(fit.par.c, function(x) x[3])), 25),
        test.m1se = rep(unlist(lapply(fit.par.mse, function(x) x[1])), 25),
        test.m2se = rep(unlist(lapply(fit.par.mse, function(x) x[2])), 25),
        test.m3se = rep(unlist(lapply(fit.par.mse, function(x) x[3])), 25),
        test.m4se = rep(unlist(lapply(fit.par.mse, function(x) x[4])), 25),
        test.c1se = rep(unlist(lapply(fit.par.cse, function(x) x[1])), 25),
        test.c2se = rep(unlist(lapply(fit.par.cse, function(x) x[2])), 25),
        test.c3se = rep(unlist(lapply(fit.par.cse, function(x) x[3])), 25),
        time.1 = c(rep(get(data.names[i])[["tbas"]][1], niter),
          get(data.names[i])[["t000"]][1,], get(data.names[i])[["t001"]][1,],
          get(data.names[i])[["t002"]][1,], get(data.names[i])[["t010"]][1,],
          get(data.names[i])[["t001"]][1,], get(data.names[i])[["t012"]][1,],
          get(data.names[i])[["t020"]][1,], get(data.names[i])[["t021"]][1,],
          get(data.names[i])[["t022"]][1,], get(data.names[i])[["t030"]][1,],
          get(data.names[i])[["t031"]][1,], get(data.names[i])[["t032"]][1,],
          get(data.names[i])[["t100"]][1,], get(data.names[i])[["t101"]][1,],
          get(data.names[i])[["t102"]][1,], get(data.names[i])[["t110"]][1,],
          get(data.names[i])[["t111"]][1,], get(data.names[i])[["t112"]][1,],
          get(data.names[i])[["t120"]][1,], get(data.names[i])[["t121"]][1,],
          get(data.names[i])[["t122"]][1,], get(data.names[i])[["t130"]][1,],
          get(data.names[i])[["t131"]][1,], get(data.names[i])[["t132"]][1,]),
        time.2 = c(rep(get(data.names[i])[["tbas"]][2], niter),
          get(data.names[i])[["t000"]][2,], get(data.names[i])[["t001"]][2,],
          get(data.names[i])[["t002"]][2,], get(data.names[i])[["t010"]][2,],
          get(data.names[i])[["t001"]][2,], get(data.names[i])[["t012"]][2,],
          get(data.names[i])[["t020"]][2,], get(data.names[i])[["t021"]][2,],
          get(data.names[i])[["t022"]][2,], get(data.names[i])[["t030"]][2,],
          get(data.names[i])[["t031"]][2,], get(data.names[i])[["t032"]][2,],
          get(data.names[i])[["t100"]][2,], get(data.names[i])[["t101"]][2,],
          get(data.names[i])[["t102"]][2,], get(data.names[i])[["t110"]][2,],
          get(data.names[i])[["t111"]][2,], get(data.names[i])[["t112"]][2,],
          get(data.names[i])[["t120"]][2,], get(data.names[i])[["t121"]][2,],
          get(data.names[i])[["t122"]][2,], get(data.names[i])[["t130"]][2,],
          get(data.names[i])[["t131"]][2,], get(data.names[i])[["t132"]][2,]),
        time.3 = c(rep(get(data.names[i])[["tbas"]][3], niter),
          get(data.names[i])[["t000"]][3,], get(data.names[i])[["t001"]][3,],
          get(data.names[i])[["t002"]][3,], get(data.names[i])[["t010"]][3,],
          get(data.names[i])[["t001"]][3,], get(data.names[i])[["t012"]][3,],
          get(data.names[i])[["t020"]][3,], get(data.names[i])[["t021"]][3,],
          get(data.names[i])[["t022"]][3,], get(data.names[i])[["t030"]][3,],
          get(data.names[i])[["t031"]][3,], get(data.names[i])[["t032"]][3,],
          get(data.names[i])[["t100"]][3,], get(data.names[i])[["t101"]][3,],
          get(data.names[i])[["t102"]][3,], get(data.names[i])[["t110"]][3,],
          get(data.names[i])[["t111"]][3,], get(data.names[i])[["t112"]][3,],
          get(data.names[i])[["t120"]][3,], get(data.names[i])[["t121"]][3,],
          get(data.names[i])[["t122"]][3,], get(data.names[i])[["t130"]][3,],
          get(data.names[i])[["t131"]][3,], get(data.names[i])[["t132"]][3,]),
        time.4 = c(rep(get(data.names[i])[["tbas"]][4], niter),
          get(data.names[i])[["t000"]][4,], get(data.names[i])[["t001"]][4,],
          get(data.names[i])[["t002"]][4,], get(data.names[i])[["t010"]][4,],
          get(data.names[i])[["t001"]][4,], get(data.names[i])[["t012"]][4,],
          get(data.names[i])[["t020"]][4,], get(data.names[i])[["t021"]][4,],
          get(data.names[i])[["t022"]][4,], get(data.names[i])[["t030"]][4,],
          get(data.names[i])[["t031"]][4,], get(data.names[i])[["t032"]][4,],
          get(data.names[i])[["t100"]][4,], get(data.names[i])[["t101"]][4,],
          get(data.names[i])[["t102"]][4,], get(data.names[i])[["t110"]][4,],
          get(data.names[i])[["t111"]][4,], get(data.names[i])[["t112"]][4,],
          get(data.names[i])[["t120"]][4,], get(data.names[i])[["t121"]][4,],
          get(data.names[i])[["t122"]][4,], get(data.names[i])[["t130"]][4,],
          get(data.names[i])[["t131"]][4,], get(data.names[i])[["t132"]][4,]),
        time.5 = c(rep(get(data.names[i])[["tbas"]][5], niter),
          get(data.names[i])[["t000"]][5,], get(data.names[i])[["t001"]][5,],
          get(data.names[i])[["t002"]][5,], get(data.names[i])[["t010"]][5,],
          get(data.names[i])[["t001"]][5,], get(data.names[i])[["t012"]][5,],
          get(data.names[i])[["t020"]][5,], get(data.names[i])[["t021"]][5,],
          get(data.names[i])[["t022"]][5,], get(data.names[i])[["t030"]][5,],
          get(data.names[i])[["t031"]][5,], get(data.names[i])[["t032"]][5,],
          get(data.names[i])[["t100"]][5,], get(data.names[i])[["t101"]][5,],
          get(data.names[i])[["t102"]][5,], get(data.names[i])[["t110"]][5,],
          get(data.names[i])[["t111"]][5,], get(data.names[i])[["t112"]][5,],
          get(data.names[i])[["t120"]][5,], get(data.names[i])[["t121"]][5,],
          get(data.names[i])[["t122"]][5,], get(data.names[i])[["t130"]][5,],
          get(data.names[i])[["t131"]][5,], get(data.names[i])[["t132"]][5,]),
        time.6 = c(rep(get(data.names[i])[["tbas"]][6], niter),
          get(data.names[i])[["t000"]][6,], get(data.names[i])[["t001"]][6,],
          get(data.names[i])[["t002"]][6,], get(data.names[i])[["t010"]][6,],
          get(data.names[i])[["t001"]][6,], get(data.names[i])[["t012"]][6,],
          get(data.names[i])[["t020"]][6,], get(data.names[i])[["t021"]][6,],
          get(data.names[i])[["t022"]][6,], get(data.names[i])[["t030"]][6,],
          get(data.names[i])[["t031"]][6,], get(data.names[i])[["t032"]][6,],
          get(data.names[i])[["t100"]][6,], get(data.names[i])[["t101"]][6,],
          get(data.names[i])[["t102"]][6,], get(data.names[i])[["t110"]][6,],
          get(data.names[i])[["t111"]][6,], get(data.names[i])[["t112"]][6,],
          get(data.names[i])[["t120"]][6,], get(data.names[i])[["t121"]][6,],
          get(data.names[i])[["t122"]][6,], get(data.names[i])[["t130"]][6,],
          get(data.names[i])[["t131"]][6,], get(data.names[i])[["t132"]][6,]),
        # time.7 = c(rep(get(data.names[i])[["tbas"]][7], niter),
        #   get(data.names[i])[["t000"]][7,], get(data.names[i])[["t001"]][7,],
        #   get(data.names[i])[["t002"]][7,], get(data.names[i])[["t010"]][7,],
        #   get(data.names[i])[["t001"]][7,], get(data.names[i])[["t012"]][7,],
        #   get(data.names[i])[["t020"]][7,], get(data.names[i])[["t021"]][7,],
        #   get(data.names[i])[["t022"]][7,], get(data.names[i])[["t030"]][7,],
        #   get(data.names[i])[["t031"]][7,], get(data.names[i])[["t032"]][7,],
        #   get(data.names[i])[["t100"]][7,], get(data.names[i])[["t101"]][7,],
        #   get(data.names[i])[["t102"]][7,], get(data.names[i])[["t110"]][7,],
        #   get(data.names[i])[["t111"]][7,], get(data.names[i])[["t112"]][7,],
        #   get(data.names[i])[["t120"]][7,], get(data.names[i])[["t121"]][7,],
        #   get(data.names[i])[["t122"]][7,], get(data.names[i])[["t130"]][7,],
        #   get(data.names[i])[["t131"]][7,], get(data.names[i])[["t132"]][7,]),
        # time.8 = c(rep(get(data.names[i])[["tbas"]][8], niter),
        #   get(data.names[i])[["t000"]][8,], get(data.names[i])[["t001"]][8,],
        #   get(data.names[i])[["t002"]][8,], get(data.names[i])[["t010"]][8,],
        #   get(data.names[i])[["t001"]][8,], get(data.names[i])[["t012"]][8,],
        #   get(data.names[i])[["t020"]][8,], get(data.names[i])[["t021"]][8,],
        #   get(data.names[i])[["t022"]][8,], get(data.names[i])[["t030"]][8,],
        #   get(data.names[i])[["t031"]][8,], get(data.names[i])[["t032"]][8,],
        #   get(data.names[i])[["t100"]][8,], get(data.names[i])[["t101"]][8,],
        #   get(data.names[i])[["t102"]][8,], get(data.names[i])[["t110"]][8,],
        #   get(data.names[i])[["t111"]][8,], get(data.names[i])[["t112"]][8,],
        #   get(data.names[i])[["t120"]][8,], get(data.names[i])[["t121"]][8,],
        #   get(data.names[i])[["t122"]][8,], get(data.names[i])[["t130"]][8,],
        #   get(data.names[i])[["t131"]][8,], get(data.names[i])[["t132"]][8,]),
        # time.9 = c(rep(get(data.names[i])[["tbas"]][9], niter),
        #   get(data.names[i])[["t000"]][9,], get(data.names[i])[["t001"]][9,],
        #   get(data.names[i])[["t002"]][9,], get(data.names[i])[["t010"]][9,],
        #   get(data.names[i])[["t001"]][9,], get(data.names[i])[["t012"]][9,],
        #   get(data.names[i])[["t020"]][9,], get(data.names[i])[["t021"]][9,],
        #   get(data.names[i])[["t022"]][9,], get(data.names[i])[["t030"]][9,],
        #   get(data.names[i])[["t031"]][9,], get(data.names[i])[["t032"]][9,],
        #   get(data.names[i])[["t100"]][9,], get(data.names[i])[["t101"]][9,],
        #   get(data.names[i])[["t102"]][9,], get(data.names[i])[["t110"]][9,],
        #   get(data.names[i])[["t111"]][9,], get(data.names[i])[["t112"]][9,],
        #   get(data.names[i])[["t120"]][9,], get(data.names[i])[["t121"]][9,],
        #   get(data.names[i])[["t122"]][9,], get(data.names[i])[["t130"]][9,],
        #   get(data.names[i])[["t131"]][9,], get(data.names[i])[["t132"]][9,]),
        tse = fit.par.tse
      )
      d.melt$prop <- with(d.melt, test/ref)
      r.out <- rbind(r.out,
        ddply(d.melt, .(type), function(x) {
          c(data = data.names[i], metric = slot.names[j], sumfuncBOX(x$prop))
        })
      )
      d.mid <- ddply(d.melt, .(type), function(x) {
        q75 <- as.numeric(r.out[r.out$type == unique(x$type), ]$q75)[j]
        q25 <- as.numeric(r.out[r.out$type == unique(x$type), ]$q25)[j]
        lower <- q25 - (q75 - q25)
        upper <- q75 + (q75 - q25)
        x$outlier <- ifelse(x$prop < lower | x$prop > upper, T, F)
        x$bioq <- ifelse(x$prop < 0.95 | x$prop > 1.05, F, T)
        x
      })
      d.out <- rbind(d.out, d.mid)
    }
    res <- rbind(res, r.out)
    plotdata <- rbind(plotdata, d.out)
  }
# -----------------------------------------------------------------------------
# Table
  broad.out <- ddply(plotdata, .(data, type, metric), function(x) {
    prop <- with(x, test/ref)
    stat1 <-  mean(prop, na.rm = T)
    stat2 <-  sd(prop, na.rm = T)
    out <- data.frame(
      mean = stat1,
      sd = stat2,
      cv = stat2/stat1*100,
      median = median(prop, na.rm = T),
      CI90lo = quantile(prop, probs = 0.05, na.rm = T, names = F),
      CI90hi = quantile(prop, probs = 0.95, na.rm = T, names = F),
      n = length(na.omit(prop))
    )
  })
  # write.csv(broad.out, "broad-output.csv")

# True vs. Pred
  plotdata$bioqf <- factor(plotdata$bioq)
  lineofid <- 10^seq(
    from = log10(min(plotdata$test[plotdata$test >= 0], na.rm = T)+0.0001),
    to = log10(max(plotdata$test[plotdata$test != Inf], na.rm = T)),
    length.out = 1000
  )
  bioqhiline <- data.frame(ref = lineofid, test = lineofid*1.11)
  bioqloline <- data.frame(ref = lineofid, test = lineofid*0.9)
  plotdata2 <- plotdata
  plotdata2$typef <- factor(plotdata2$type)
  levels(plotdata2$typef) <- c("otter", "otter+glam", "otter+olam",
    "otter+auct", "otter+auct+glam", "otter+auct+olam",
    "otter+halft", "otter+halft+glam", "otter+halft+olam",
    "otter+obst", "otter+obst+glam", "otter+obst+olam",
    "otter+tmax", "otter+tmax+glam", "otter+tmax+olam",
    "otter+tmax+auct", "otter+tmax+auct+glam",
    "otter+tmax+auct+olam", "otter+tmax+halft",
    "otter+tmax+halft+glam", "otter+tmax+halft+olam", "otter+tmax+obst",
    "otter+tmax+obst+glam", "otter+tmax+obst+olam", "basic")
  statdata <- ddply(plotdata2, .(data, metric, typef), function(x) {
    data.frame(
      percbioq = paste0("bioq = ", signif(sum(x$bioq, na.rm = T)/length(x$bioq)*100, 3), "%"),
      meanratio = paste("mean =", round(mean(x$prop), 2))
    )
  })

  truevpred.plot <- function(data, metric, n.col) {
    plotsub <- plotdata2[plotdata2$data == data & plotdata2$metric == metric, ]
    statsub <- statdata[statdata$data == data & statdata$metric == metric, ]
    allval <- c(plotsub$ref)
    minval <- min(allval[allval != 0], na.rm = T)
    maxval <- max(allval[allval != Inf], na.rm = T)
    valrange <- log10(maxval) - log10(minval)
    p0 <- NULL
    p0 <- ggplot()
    p0 <- p0 + ggtitle(paste("True vs. Predicted", metric, "-", data))
    p0 <- p0 + geom_point(aes(x = ref, y = test, colour = bioqf), data = plotsub, alpha = 0.2)
    p0 <- p0 + geom_line(aes(x = ref, y = test), data = bioqhiline,
      colour = "green4", linetype = "dashed")
    p0 <- p0 + geom_line(aes(x = ref, y = test), data = bioqloline,
      colour = "green4", linetype = "dashed")
    p0 <- p0 + geom_abline(intercept = 0, slope = 1)
    p0 <- p0 + geom_text(aes(label = percbioq), data = statsub,
      x = log10(minval)+valrange*0.29, y = log10(maxval)-valrange*0.11)
    p0 <- p0 + geom_text(aes(label = meanratio), data = statsub,
      x = log10(maxval)-valrange*0.29, y = log10(minval)+valrange*0.11)
    p0 <- p0 + scale_x_log10("True AUC", lim = c(minval, maxval))
    p0 <- p0 + scale_y_log10("Predicted AUC", lim = c(minval, maxval))
    p0 <- p0 + scale_colour_manual(values = c("red", "blue"))
    p0 <- p0 + facet_wrap(~typef, ncol = n.col)
    p0
  }

  p11 <- truevpred.plot("d2b", "auc24", 6)
  p12 <- truevpred.plot("d3b", "auc24", 6)
  p13 <- truevpred.plot("d1a", "auc24", 6)
  p14 <- truevpred.plot("d2a", "auc24", 6)
  p15 <- truevpred.plot("d3a", "auc24", 6)
  p21 <- truevpred.plot("d2b", "auctlast", 6)
  p22 <- truevpred.plot("d3b", "auctlast", 6)
  p23 <- truevpred.plot("d1a", "auctlast", 6)
  p24 <- truevpred.plot("d2a", "auctlast", 6)
  p25 <- truevpred.plot("d3a", "auctlast", 6)
  p31 <- truevpred.plot("d2b", "aucinf", 6)
  p32 <- truevpred.plot("d3b", "aucinf", 6)
  p33 <- truevpred.plot("d1a", "aucinf", 6)
  p34 <- truevpred.plot("d2a", "aucinf", 6)
  p35 <- truevpred.plot("d3a", "aucinf", 6)
  p41 <- truevpred.plot("d1a", "cmax", 6)
  p42 <- truevpred.plot("d2a", "cmax", 6)
  p43 <- truevpred.plot("d3a", "cmax", 6)
  p51 <- truevpred.plot("d1a", "tmax", 6)
  p52 <- truevpred.plot("d2a", "tmax", 6)
  p53 <- truevpred.plot("d3a", "tmax", 6)

  plotdata2 <- plotdata[plotdata$type %in% c("000", "030", "100", "130", "bas"), ]
  plotdata2$typef <- factor(plotdata2$type)
  levels(plotdata2$typef) <- c("otter", "otter+obst", "otter+tmax", "otter+tmax+obst", "basic")
  statdata <- ddply(plotdata2, .(data, metric, typef), function(x) {
    data.frame(
      percbioq = paste0("bioq = ", signif(sum(x$bioq, na.rm = T)/length(x$bioq)*100, 3), "%"),
      meanratio = paste("mean =", round(mean(x$prop), 2))
    )
  })

  s11 <- truevpred.plot("d2b", "auc24", 6)
  s12 <- truevpred.plot("d3b", "auc24", 6)
  s13 <- truevpred.plot("d1a", "auc24", 6)
  s14 <- truevpred.plot("d2a", "auc24", 6)
  s15 <- truevpred.plot("d3a", "auc24", 6)
  s21 <- truevpred.plot("d2b", "auctlast", 3)
  s22 <- truevpred.plot("d3b", "auctlast", 3)
  s23 <- truevpred.plot("d1a", "auctlast", 3)
  s24 <- truevpred.plot("d2a", "auctlast", 3)
  s25 <- truevpred.plot("d3a", "auctlast", 3)
  s31 <- truevpred.plot("d2b", "aucinf", 3)
  s32 <- truevpred.plot("d3b", "aucinf", 3)
  s33 <- truevpred.plot("d1a", "aucinf", 3)
  s34 <- truevpred.plot("d2a", "aucinf", 3)
  s35 <- truevpred.plot("d3a", "aucinf", 3)
  s41 <- truevpred.plot("d1a", "cmax", 3)
  s42 <- truevpred.plot("d2a", "cmax", 3)
  s43 <- truevpred.plot("d3a", "cmax", 3)
  s51 <- truevpred.plot("d1a", "tmax", 3)
  s52 <- truevpred.plot("d2a", "tmax", 3)
  s53 <- truevpred.plot("d3a", "tmax", 3)

  # plotsub <- plotdata[plotdata$data == "d2b" & plotdata$metric == "auc24", ]
  # minval <- min(c(plotsub$test, plotsub$ref), na.rm = T)
  # maxval <- max(c(plotsub$test, plotsub$ref), na.rm = T)
  # p01 <- NULL
  # p01 <- ggplot()
  # p01 <- p01 + ggtitle("True vs. Predicted AUC (0-24) - d2b")
  # p01 <- p01 + geom_point(aes(x = ref, y = test, colour = bioqf), data = plotsub)
  # p01 <- p01 + geom_line(aes(x = ref, y = test), data = bioqhiline, colour = "green4", linetype = "dashed")
  # p01 <- p01 + geom_line(aes(x = ref, y = test), data = bioqloline, colour = "green4", linetype = "dashed")
  # p01 <- p01 + geom_abline(intercept = 0, slope = 1)
  # p01 <- p01 + scale_x_log10("True AUC", lim = c(minval, maxval))
  # p01 <- p01 + scale_y_log10("Predicted AUC", lim = c(minval, maxval))
  # p01 <- p01 + scale_colour_manual(values = c("red", "blue"))
  # p01 <- p01 + facet_wrap(~typef, ncol = 6)
  # p01

# Box plot
  subdata <- plotdata
  subdata$typef <- as.factor(subdata$type)
  levels(subdata$typef) <- c("otter", "otter+glam", "otter+olam",
    "otter+auct", "otter+auct+glam", "otter+auct+olam",
    "otter+halft", "otter+halft+glam", "otter+halft+olam",
    "otter+obst", "otter+obst+glam", "otter+obst+olam",
    "otter+tmax", "otter+tmax+glam", "otter+tmax+olam",
    "otter+tmax+auct", "otter+tmax+auct+glam",
    "otter+tmax+auct+olam", "otter+tmax+halft",
    "otter+tmax+halft+glam", "otter+tmax+halft+olam", "otter+tmax+obst",
    "otter+tmax+obst+glam", "otter+tmax+obst+olam", "basic")
  subdata <- ddply(subdata, .(data, typef, metric), function(x) {
    data.frame(
      prop = x$prop,
      CI90lo = quantile(x$prop, probs = 0.05, na.rm = T, names = F),
      CI90hi = quantile(x$prop, probs = 0.95, na.rm = T, names = F),
      CI95lo = quantile(x$prop, probs = 0.025, na.rm = T, names = F),
      CI95hi = quantile(x$prop, probs = 0.975, na.rm = T, names = F)
    )
  })
  plotauc24 <- subdata[subdata$metric == "auc24", ]
  plotauctlast <- subdata[subdata$metric == "auctlast", ]
  plotaucinf <- subdata[subdata$metric == "aucinf", ]
  plotcmax <- subdata[subdata$metric == "cmax" & subdata$data %in% c("d1a", "d2a", "d3a"), ]
  plottmax <- subdata[subdata$metric == "tmax" & subdata$data %in% c("d1a", "d2a", "d3a"), ]

  p1 <- NULL
  p1 <- ggplot(plotauc24, aes(x = typef, y = prop))
  p1 <- p1 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  p1 <- p1 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # p1 <- p1 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  p1 <- p1 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  p1 <- p1 + ggtitle("Broad Study AUC (0-24)")
  p1 <- p1 + xlab("\nMethod")
  p1 <- p1 + scale_y_continuous("Method/Reference Ratio\n", lim = c(0, 5))
  p1 <- p1 + facet_wrap(~data, nrow = 1)

  p2 <- ggplot(plotauctlast, aes(x = typef, y = prop))
  p2 <- p2 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  p2 <- p2 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # p2 <- p2 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  p2 <- p2 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  p2 <- p2 + ggtitle("Broad Study AUC (0-tlast)")
  p2 <- p2 + xlab("\nMethod")
  p2 <- p2 + scale_y_continuous("Method/Reference Ratio\n", lim = c(0, 5))
  p2 <- p2 + facet_wrap(~data, nrow = 1)

  p3 <- ggplot(plotaucinf, aes(x = typef, y = prop))
  p3 <- p3 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  p3 <- p3 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # p3 <- p3 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  p3 <- p3 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  p3 <- p3 + ggtitle("Broad Study AUC (0-Inf)")
  p3 <- p3 + xlab("\nMethod")
  p3 <- p3 + scale_y_continuous("Method/Reference Ratio\n", lim = c(0, 5))
  p3 <- p3 + facet_wrap(~data, nrow = 1)

  p4 <- ggplot(plotcmax, aes(x = typef, y = prop))
  p4 <- p4 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  p4 <- p4 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # p4 <- p4 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  p4 <- p4 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  p4 <- p4 + ggtitle("Broad Study Cmax")
  p4 <- p4 + xlab("\nMethod")
  p4 <- p4 + scale_y_continuous("Method/Reference Ratio\n", lim = c(0, 5))
  p4 <- p4 + facet_wrap(~data, nrow = 1)

  p5 <- ggplot(plottmax, aes(x = typef, y = prop))
  p5 <- p5 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  p5 <- p5 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # p5 <- p5 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  p5 <- p5 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  p5 <- p5 + ggtitle("Broad Study Tmax")
  p5 <- p5 + xlab("\nMethod")
  p5 <- p5 + scale_y_continuous("Method/Reference Ratio\n", lim = c(0, 5))
  p5 <- p5 + facet_wrap(~data, nrow = 1)

  s1 <- NULL
  s1 <- ggplot(plotauc24[plotauc24$typef %in% c("otter", "otter+obst", "otter+tmax", "otter+tmax+obst", "basic"), ], aes(x = typef, y = prop))
  s1 <- s1 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  s1 <- s1 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # s1 <- s1 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  s1 <- s1 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  s1 <- s1 + ggtitle("Broad Study AUC (0-24)")
  s1 <- s1 + xlab("\nMethod")
  s1 <- s1 + ylab("Method/Reference Ratio\n")
  s1 <- s1 + facet_wrap(~data, nrow = 1)

  s2 <- ggplot(plotauctlast[plotauctlast$typef %in% c("otter", "otter+obst", "otter+tmax", "otter+tmax+obst", "basic"), ], aes(x = typef, y = prop))
  s2 <- s2 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  s2 <- s2 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # s2 <- s2 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  s2 <- s2 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  s2 <- s2 + ggtitle("Broad Study AUC (0-tlast)")
  s2 <- s2 + xlab("\nMethod")
  s2 <- s2 + ylab("Method/Reference Ratio\n")
  s2 <- s2 + facet_wrap(~data, nrow = 1)

  s3 <- ggplot(plotaucinf[plotaucinf$typef %in% c("otter", "otter+obst", "otter+tmax", "otter+tmax+obst", "basic"), ], aes(x = typef, y = prop))
  s3 <- s3 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  s3 <- s3 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # s3 <- s3 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  s3 <- s3 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  s3 <- s3 + ggtitle("Broad Study AUC (0-Inf)")
  s3 <- s3 + xlab("\nMethod")
  s3 <- s3 + ylab("Method/Reference Ratio\n")
  s3 <- s3 + facet_wrap(~data, nrow = 1)

  s4 <- ggplot(plotcmax[plotcmax$typef %in% c("otter", "otter+obst", "otter+tmax", "otter+tmax+obst", "basic"), ], aes(x = typef, y = prop))
  s4 <- s4 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  s4 <- s4 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # s4 <- s4 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  s4 <- s4 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  s4 <- s4 + ggtitle("Broad Study Cmax")
  s4 <- s4 + xlab("\nMethod")
  s4 <- s4 + ylab("Method/Reference Ratio\n")
  s4 <- s4 + facet_wrap(~data, nrow = 1)

  s5 <- ggplot(plottmax[plottmax$typef %in% c("otter", "otter+obst", "otter+tmax", "otter+tmax+obst", "basic"), ], aes(x = typef, y = prop))
  s5 <- s5 + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
  s5 <- s5 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # s5 <- s5 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  s5 <- s5 + geom_boxplot(aes(ymin = CI95lo, ymax = CI95hi))
  s5 <- s5 + ggtitle("Broad Study Tmax")
  s5 <- s5 + xlab("\nMethod")
  s5 <- s5 + ylab("Method/Reference Ratio\n")
  s5 <- s5 + facet_wrap(~data, nrow = 1)


# Attempt at new plots using %difference from truth
  plotdata3 <- plotdata[plotdata$type %in% c("000", "030", "bas"), ]
  plotdata3$diff <- abs(1 - plotdata3$prop)
  plotdata3 <- ddply(plotdata3, .(type, metric, data), function(x) {
    x$diffrank <- rank(x$diff, ties.method = "first")/1000
    x$diffrank
    x
  })

  p6 <- NULL
  p6 <- ggplot(data = plotdata3[plotdata3$metric == "aucinf" & plotdata3$data == "d2a", ])
  p6 <- p6 + geom_line(aes(x = diffrank, y = diff*100, colour = type), size = 0.8)
  p6 <- p6 + ggtitle("Percent Difference From Truth - AUC")
  p6 <- p6 + xlab("\nRanked PK Profiles")
  p6 <- p6 + ylab("Difference in AUC (%)\n")
  p6

  p7 <- NULL
  p7 <- ggplot(data = plotdata3[plotdata3$metric == "aucinf", ])
  p7 <- p7 + geom_line(aes(x = diffrank, y = diff*100, colour = type), size = 0.8)
  p7 <- p7 + ggtitle("Percent Difference From Truth - AUC")
  p7 <- p7 + xlab("\nRanked PK Profiles")
  p7 <- p7 + ylab("Difference in AUC (%)\n")
  p7 <- p7 + facet_wrap(~data, nrow = 1)
  p7

# Attempt at new plots using difference in % difference from Truth
  basdata <- plotdata3[plotdata3$type == "bas", ]
  plotdata3 <- ddply(plotdata3, .(type, metric, data), function(x) {
    x.metric <- unique(x$metric)
    x.data <- unique(x$data)
    bassub <- basdata[basdata$metric == x.metric & basdata$data == x.data, ]
    x$diffbase <- bassub$diff - x$diff
    x$diffbaserank <- rank(x$diffbase, ties.method = "first")/1000
    if (unique(x$type) == "bas") {
      x$diffbasemed <- NA
    } else {
      x$diffbasemed <- x$diffbaserank[which(abs(x$diffbase) == min(abs(x$diffbase), na.rm = T))][1]
    }
    x
  })

  p8 <- NULL
  p8 <- ggplot(data = plotdata3[plotdata3$metric == "aucinf" & plotdata3$data == "d1a", ])
  p8 <- p8 + geom_line(aes(x = diffbaserank, y = diffbase*100, colour = type), size = 0.8)
  p8 <- p8 + ggtitle("Percent Difference From Empirical - AUC")
  p8 <- p8 + xlab("\nRanked PK Profiles")
  p8 <- p8 + ylab("Difference in AUC (%)\n")
  p8

  p9 <- NULL
  p9 <- ggplot(data = plotdata3[plotdata3$metric == "aucinf", ])
  p9 <- p9 + geom_line(aes(x = diffbaserank, y = diffbase*100, colour = type), size = 0.8)
  p9 <- p9 + geom_vline(aes(xintercept = diffbasemed, colour = type))
  p9 <- p9 + ggtitle("Percent Difference From Empirical - AUC")
  p9 <- p9 + xlab("\nRanked PK Profiles")
  p9 <- p9 + ylab("Difference in AUC (%)\n")
  p9 <- p9 + facet_wrap(~data, nrow = 1)
  p9
