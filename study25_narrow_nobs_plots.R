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
  library(cowplot)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Source functions
  source(paste(git.dir, reponame, "sumstat_functions.R", sep = "/"))

# Source files of interest
  nobs <- 6:12
  run.string <- paste0("newnobs", nobs, "-AR3024")

# -----------------------------------------------------------------------------
# Data structure
  data.names <- c("d2b", "d3b", "d1a", "d2a", "d3a")
  slot.names <- c("auc24", "auctlast", "aucinf", "cmax", "tmax")
  niter <- 1000

  alldata <- data.frame(NULL)
  for (k in 1:length(run.string)) {
    d2b <- readRDS(paste(git.dir, reponame,
      paste0("fn_diag/d2b-", run.string[k], ".rds"), sep = "/"))
    d3b <- readRDS(paste(git.dir, reponame,
      paste0("fn_diag/d3b-", run.string[k], ".rds"), sep = "/"))
    d1a <- readRDS(paste(git.dir, reponame,
      paste0("fn_diag/d1a-", run.string[k], ".rds"), sep = "/"))
    d2a <- readRDS(paste(git.dir, reponame,
      paste0("fn_diag/d2a-", run.string[k], ".rds"), sep = "/"))
    d3a <- readRDS(paste(git.dir, reponame,
      paste0("fn_diag/d3a-", run.string[k], ".rds"), sep = "/"))
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
    plotdata$bioqf <- factor(plotdata$bioq)
    plotdata$nobs <- nobs[k]
    plotdata$nobsf <- factor(nobs[k])
    alldata <- rbind(alldata, plotdata)
  }

# -----------------------------------------------------------------------------
centre_legend <- function(p, draw = T) {
# Function created from SO thread https://stackoverflow.com/questions/48000292
# Input: ggplot or grob
# Output: draws the plot (or a grob if draw = F)
# This is the last thing you do with your plot as your output will no longer
# be a ggplot object!
# Load required libraries
  require(ggplot2)
  require(gtable)
  require(grid)
# Turn plot into a grob if not already in that format
  if ("gg" %in% class(p)) g <- ggplotGrob(p)
  else if ("grob" %in% class(p)) g <- p
  else stop("ggplots & grobs are the only suitable input for this function")
# Extract legend
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
# Extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  guides <- legend$grobs[[guides_index]]
# Add extra column to space the legend, effectively making it centred
# guides$width[5] is the extra spacing from the end of the legend text
# to the end of the legend title. If we instead distribute it 50:50 on
# both sides, we get a centred legend.
  guides <- gtable_add_cols(guides, 0.5*guides$width[5], 1)
  guides$widths[6] <- guides$widths[2]
  title_index <- guides$layout$name == "title"
  guides$layout$l[title_index] <- 2
# Reconstruct legend and write back
  legend$grobs[[guides_index]] <- guides
  g$grobs[[legend_index]] <- legend
  if (draw) {
    grid.newpage()
    grid.draw(g)
  } else g
}
# Attempt at new plots using %difference from truth
  basdata <- alldata[alldata$type == "bas", ]
  basdata$diff <- abs(1 - basdata$prop)
  subdata <- alldata[alldata$type == "000", ]
  subdata$diff <- abs(1 - subdata$prop)
  subdata <- ddply(subdata, .(nobs, metric, data), function(x) {
    x$diffrank <- rank(x$diff, ties.method = "first")/1000
    x$diffrank
    x
  })

  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = id/1000, y = diff*100),
    data = basdata[basdata$metric == "aucinf", ], size = 0.8)
  p1 <- p1 + geom_line(aes(x = diffrank, y = diff*100, group = nobs, colour = nobs),
    data = subdata[subdata$metric == "aucinf", ], size = 0.8)
  p1 <- p1 + ggtitle("Percent Difference From Truth - AUC")
  p1 <- p1 + xlab("\nRanked PK Profiles")
  p1 <- p1 + ylab("Difference in AUC (%)\n")
  p1 <- p1 + coord_cartesian(ylim = c(0, 30))
  p1 <- p1 + facet_wrap(~data, nrow = 1)
  p1 <- p1 + scale_colour_gradient(trans = "log", lim = c(5.75, 12.8), breaks = 3:6*2)
  p1

  basdata$dataf <- factor(basdata$data)
  levels(basdata$dataf) <- c("2 compartment bolus", "3 compartment bolus", "1 compartment oral", "2 compartment oral", "3 compartment oral")
  subdata$dataf <- factor(subdata$data)
  levels(subdata$dataf) <- c("2 compartment bolus", "3 compartment bolus", "1 compartment oral", "2 compartment oral", "3 compartment oral")

  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + geom_vline(aes(xintercept = prop), linetype = "dashed", colour = "red",
    data = basdata[basdata$metric == "aucinf", ], size = 0.8)
  p2 <- p2 + geom_vline(xintercept = 1, size = 0.8)
  p2 <- p2 + geom_line(aes(prop, group = nobs, colour = nobs), stat = "density",
    data = subdata[subdata$metric == "aucinf", ], size = 0.8, alpha = 0.7)
  p2 <- p2 + scale_colour_gradientn(trans = "log", lim = c(5.75, 12.8),
    colours = c("red","yellow","green","lightblue","darkblue","blueviolet"), breaks = 3:6*2)
  p2 <- p2 + theme(legend.position = c(0.85, 0.2), legend.background = element_rect(fill = "white", colour = NA))
  p2 <- p2 + labs(x = "\nRelative AUC", y = "Density\n", colour = "  Number of\nObservations")
  p2 <- p2 + coord_cartesian(xlim = c(0.9, 1.1))
  p2 <- p2 + facet_wrap(~dataf, nrow = 2, scale = "free_y")
  centre_legend(p2)

  png("Figure 2c.png", width = 23.2, height = 11.2, units = "cm", res = 300)
  centre_legend(p2)
  dev.off()
