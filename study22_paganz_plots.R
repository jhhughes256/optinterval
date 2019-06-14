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
  run.string <- "broadnewdt2-AR3024"
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
    ref.par.nexp <- ceiling(dim(get(data.names[i])[["par"]])[1]/2)
    ref.par.m <- apply(get(data.names[i])[["par"]], 2, function(x) {
      x[1:ceiling(length(x)/2)]
    })
    ref.par.c <- apply(get(data.names[i])[["par"]], 2, function(x) {
      x[-(1:ceiling(length(x)/2))]
    })
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
        ref.m1 = rep(ref.par.m[1,], 25),
        ref.m2 = rep(ref.par.m[2,], 25),
        ref.m3 = if (dim(ref.par.m)[1] > 2) {
            rep(ref.par.m[3,], 25)
          } else {
            NA
          },
        ref.m4 = if (dim(ref.par.m)[1] > 3) {
            rep(ref.par.m[4,], 25)
          } else {
            NA
          },
        ref.c1 = if (is.vector(ref.par.c)) {
            rep(ref.par.c, 25)
          } else {
            rep(ref.par.c[1,], 25)
          },
        ref.c2 = if (is.vector(ref.par.c)) {
            NA
          } else if (dim(ref.par.c)[1] > 1) {
            rep(ref.par.c[2,], 25)
          } else {
            NA
          },
        ref.c3 = if (is.vector(ref.par.c)) {
            NA
          } else if (dim(ref.par.c)[1] > 2) {
            rep(ref.par.c[3,], 25)
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
        time.7 = c(rep(get(data.names[i])[["tbas"]][7], niter),
          get(data.names[i])[["t000"]][7,], get(data.names[i])[["t001"]][7,],
          get(data.names[i])[["t002"]][7,], get(data.names[i])[["t010"]][7,],
          get(data.names[i])[["t001"]][7,], get(data.names[i])[["t012"]][7,],
          get(data.names[i])[["t020"]][7,], get(data.names[i])[["t021"]][7,],
          get(data.names[i])[["t022"]][7,], get(data.names[i])[["t030"]][7,],
          get(data.names[i])[["t031"]][7,], get(data.names[i])[["t032"]][7,],
          get(data.names[i])[["t100"]][7,], get(data.names[i])[["t101"]][7,],
          get(data.names[i])[["t102"]][7,], get(data.names[i])[["t110"]][7,],
          get(data.names[i])[["t111"]][7,], get(data.names[i])[["t112"]][7,],
          get(data.names[i])[["t120"]][7,], get(data.names[i])[["t121"]][7,],
          get(data.names[i])[["t122"]][7,], get(data.names[i])[["t130"]][7,],
          get(data.names[i])[["t131"]][7,], get(data.names[i])[["t132"]][7,]),
        time.8 = c(rep(get(data.names[i])[["tbas"]][8], niter),
          get(data.names[i])[["t000"]][8,], get(data.names[i])[["t001"]][8,],
          get(data.names[i])[["t002"]][8,], get(data.names[i])[["t010"]][8,],
          get(data.names[i])[["t001"]][8,], get(data.names[i])[["t012"]][8,],
          get(data.names[i])[["t020"]][8,], get(data.names[i])[["t021"]][8,],
          get(data.names[i])[["t022"]][8,], get(data.names[i])[["t030"]][8,],
          get(data.names[i])[["t031"]][8,], get(data.names[i])[["t032"]][8,],
          get(data.names[i])[["t100"]][8,], get(data.names[i])[["t101"]][8,],
          get(data.names[i])[["t102"]][8,], get(data.names[i])[["t110"]][8,],
          get(data.names[i])[["t111"]][8,], get(data.names[i])[["t112"]][8,],
          get(data.names[i])[["t120"]][8,], get(data.names[i])[["t121"]][8,],
          get(data.names[i])[["t122"]][8,], get(data.names[i])[["t130"]][8,],
          get(data.names[i])[["t131"]][8,], get(data.names[i])[["t132"]][8,]),
        time.9 = c(rep(get(data.names[i])[["tbas"]][9], niter),
          get(data.names[i])[["t000"]][9,], get(data.names[i])[["t001"]][9,],
          get(data.names[i])[["t002"]][9,], get(data.names[i])[["t010"]][9,],
          get(data.names[i])[["t001"]][9,], get(data.names[i])[["t012"]][9,],
          get(data.names[i])[["t020"]][9,], get(data.names[i])[["t021"]][9,],
          get(data.names[i])[["t022"]][9,], get(data.names[i])[["t030"]][9,],
          get(data.names[i])[["t031"]][9,], get(data.names[i])[["t032"]][9,],
          get(data.names[i])[["t100"]][9,], get(data.names[i])[["t101"]][9,],
          get(data.names[i])[["t102"]][9,], get(data.names[i])[["t110"]][9,],
          get(data.names[i])[["t111"]][9,], get(data.names[i])[["t112"]][9,],
          get(data.names[i])[["t120"]][9,], get(data.names[i])[["t121"]][9,],
          get(data.names[i])[["t122"]][9,], get(data.names[i])[["t130"]][9,],
          get(data.names[i])[["t131"]][9,], get(data.names[i])[["t132"]][9,]),
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
        x$bioq <- ifelse(x$prop < 0.9 | x$prop > 1.1, T, F)
        x$bioq[x$metric == "tmax"] <- ifelse(x$prop[x$metric == "tmax"] < 0.8 | x$prop[x$metric == "tmax"] > 1.2, T, F)
        x
      })
      d.out <- rbind(d.out, d.mid)
    }
    res <- rbind(res, r.out)
    plotdata <- rbind(plotdata, d.out)
  }
# -----------------------------------------------------------------------------
# Plots
  setwd("E:/Hughes/Git/splines/fn_diag")
# True vs. Pred
  plotdata$bioqf <- factor(plotdata$bioq)
  lineofid <- 10^seq(
    from = log10(min(plotdata$test[plotdata$test >= 0], na.rm = T)+0.0001),
    to = log10(max(plotdata$test[plotdata$test != Inf], na.rm = T)),
    length.out = 1000
  )
  bioqhiline <- data.frame(ref = lineofid, test = lineofid*1.25)
  bioqloline <- data.frame(ref = lineofid, test = lineofid*0.8)
  bioqhiline2 <- data.frame(ref = lineofid, test = lineofid*1.1)
  bioqloline2 <- data.frame(ref = lineofid, test = lineofid*0.9)
  plotdata2 <- plotdata[plotdata$type %in% c("000", "030"),]
  plotdata2$typef <- factor(plotdata2$type)
  levels(plotdata2$typef) <- c("Base", "Adjunct")
  statdata <- ddply(plotdata2, .(data, metric, typef), function(x) {
    data.frame(
      percbioq = paste0("bioq = ", signif(sum(x$bioq, na.rm = T)/length(x$bioq)*100, 3), "%"),
      meanratio = paste("mean =", round(mean(x$prop), 2))
    )
  })

  truevpred.plot <- function(data, metric, n.col, scale) {
    plotsub <- plotdata2[plotdata2$data == data & plotdata2$metric == metric, ]
    statsub <- statdata[statdata$data == data & statdata$metric == metric, ]
    allval <- c(plotsub$ref)
    minval <- min(allval[allval != 0], na.rm = T)
    maxval <- max(allval[allval != Inf], na.rm = T)
    valrange <- log10(maxval) - log10(minval)
    p0 <- NULL
    p0 <- ggplot()
    p0 <- p0 + geom_point(aes(x = ref, y = test, colour = bioqf), data = plotsub, alpha = 0.2)
    if (metric == "tmax") {
      p0 <- p0 + geom_line(aes(x = ref, y = test), data = bioqhiline,
        colour = "green4", linetype = "dashed")
      p0 <- p0 + geom_line(aes(x = ref, y = test), data = bioqloline,
        colour = "green4", linetype = "dashed")
      p0 <- p0 + scale_x_log10("True tmax", lim = c(minval, maxval), breaks = scale)
      p0 <- p0 + scale_y_log10("Predicted tmax", lim = c(minval, maxval), breaks = scale)
    } else {
      p0 <- p0 + geom_line(aes(x = ref, y = test), data = bioqhiline2,
        colour = "green4", linetype = "dashed")
      p0 <- p0 + geom_line(aes(x = ref, y = test), data = bioqloline2,
        colour = "green4", linetype = "dashed")
      if (metric == "cmax") {
        p0 <- p0 + scale_x_log10("True Cmax", lim = c(minval, maxval), breaks = scale)
        p0 <- p0 + scale_y_log10("Predicted Cmax", lim = c(minval, maxval), breaks = scale)
      } else {
        p0 <- p0 + scale_x_log10("True AUC", lim = c(minval, maxval), breaks = scale, labels = scales::comma)
        p0 <- p0 + scale_y_log10("Predicted AUC", lim = c(minval, maxval), breaks = scale, labels = scales::comma)
      }
    }
    p0 <- p0 + geom_abline(intercept = 0, slope = 1)
    p0 <- p0 + scale_colour_manual(values = c("blue", "red"))
    p0 <- p0 + theme(legend.position = "none")
    p0 <- p0 + facet_wrap(~typef, ncol = n.col)
    p0
  }

  # truevpred.plot("d3a", "aucinf", 2)
  # ggsave("truevpred_paganz.png", width = 23.2, height = 11.2, units = "cm")
  #
  # theme_update(axis.text.x = element_text(angle = 335))
  #
  # p1 <- truevpred.plot("d2b", "aucinf", 2, c(10, 100, 1000, 10000))
  # p2 <- truevpred.plot("d3b", "aucinf", 2, c(10, 1000, 100000))
  # p3 <- truevpred.plot("d3a", "aucinf", 2, c(3, 300, 30000))
  # p4 <- truevpred.plot("d3a", "cmax", 2, c(1, 10, 100))
  # p5 <- truevpred.plot("d3a", "tmax", 2, c(1, 3, 10))
  #
  # plot_grid(p3, p4, p5, p2, ncol = 3, label_fontface = "plain")
  # ggsave("truevpred_paganz_cow_adj.png", width = 24.2, height = 12.2, units = "cm")
  #
  # theme_bw2 <- theme_set(theme_bw(base_size = 14))
  # theme_update(plot.title = element_text(hjust = 0.5))

  plotdata2$dataf <- factor(plotdata2$data)
  levels(plotdata2$dataf) <- c("2-comp bolus", "3-comp bolus", "1-comp w/ abs", "2-comp w/ abs", "3-comp w/ abs")
  subdata <- ddply(plotdata2, .(data, dataf, typef, metric), function(x) {
    data.frame(
      prop = x$prop,
      CI90lo = quantile(x$prop, probs = 0.05, na.rm = T, names = F),
      CI75lo = quantile(x$prop, probs = 0.25, na.rm = T, names = F),
      CI50 = median(x$prop),
      CI75hi = quantile(x$prop, probs = 0.75, na.rm = T, names = F),
      CI90hi = quantile(x$prop, probs = 0.95, na.rm = T, names = F)
    )
  })
  plotaucinf <- subdata[subdata$metric == "aucinf", ]
  plotcmax <- subdata[subdata$metric == "cmax" & subdata$data %in% c("d1a", "d2a", "d3a"), ]
  plottmax <- subdata[subdata$metric == "tmax" & subdata$data %in% c("d1a", "d2a", "d3a"), ]

  p7 <- NULL
  p7 <- ggplot(plotaucinf, aes(x = typef, y = prop))
  p7 <- p7 + geom_hline(yintercept = c(0.9, 1.1), color = "green4", linetype = "dashed")
  p7 <- p7 + geom_hline(yintercept = 1, color = "red", linetype = "dashed")
  # p7 <- p7 + geom_boxplot(aes(ymin = CI90lo, ymax = CI90hi))
  p7 <- p7 + geom_boxplot(stat = "identity",
    aes(ymin = CI90lo, lower = CI75lo, middle = CI50, upper = CI75hi, ymax = CI90hi)
  )
  p7 <- p7 + xlab("\nMethod")
  p7 <- p7 + scale_y_continuous("Method/Reference AUC Ratio\n", breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2))
  p7 <- p7 + facet_wrap(~dataf, nrow = 1)
  p7 + coord_cartesian(ylim = c(0.55, 1.15))


  ggsave("boxplot_paganz_adj.png", width = 23.2, height = 11.2, units = "cm")
