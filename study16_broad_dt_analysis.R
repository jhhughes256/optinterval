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
  run.string <- "dt15-RMSE"
  d2b <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d2b-broad-", run.string, ".rds"), sep = "/"))
  d3b <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d3b-broad-", run.string, ".rds"), sep = "/"))
  d1a <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d1a-broad-", run.string, ".rds"), sep = "/"))
  d2a <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d2a-broad-", run.string, ".rds"), sep = "/"))
  d3a <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/d3a-broad-", run.string, ".rds"), sep = "/"))

# -----------------------------------------------------------------------------
# Data structure
  data.names <- c("d2b", "d3b", "d1a", "d2a", "d3a")
  slot.names <- c("auc", "cmax", "tmax")
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
    if (i == 1) {
      fit.par.mse <- list(NULL)
      fit.par.cse <- list(NULL)
      fit.par.tse <- list(NULL)
      for (k in 1:niter) {
        fit.par.mse[[k]] <- rep(NA, 4)
        fit.par.cse[[k]] <- rep(NA, 4)
        fit.par.tse[[k]] <- NA
      }
    } else {
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
      fit.par.tse <- lapply(get(data.names[i])[["interv"]], function(x) {
        vc_mat <- try(solve(x$hessian))
        if(class(vc_mat) != "try-error") {
          se <- sqrt(diag(vc_mat))
          se_percent <- abs(se/x$par*100)
          max(se_percent, na.rm = T)
        } else {
          NA
        }
      })
    }
    for (j in 1:3) {
      d.in <- get(data.names[i])[[slot.names[j]]]
      d.melt <- data.frame(
        id = rep(1:niter, 3),
        data = data.names[i],
        metric = slot.names[j],
        ref = rep(d.in$true, 3),
        test = c(d.in$basic,
          d.in$optint,
          d.in$optintwCmax),
        type = c(rep("bas", niter), rep("opt", niter), rep("optc", niter)),
        ref.m1 = ref.par.m[1,],
        ref.m2 = ref.par.m[2,],
        ref.m3 = if (dim(ref.par.m)[1] > 2) {
            ref.par.m[3,]
          } else {
            NA
          },
        ref.m4 = if (dim(ref.par.m)[1] > 3) {
            ref.par.m[4,]
          } else {
            NA
          },
        ref.c1 = if (is.vector(ref.par.c)) {
            ref.par.c
          } else {
            ref.par.c[1,]
          },
        ref.c2 = if (is.vector(ref.par.c)) {
            NA
          } else if (dim(ref.par.c)[1] > 1) {
            ref.par.c[2,]
          } else {
            NA
          },
        ref.c3 = if (is.vector(ref.par.c)) {
            NA
          } else if (dim(ref.par.c)[1] > 2) {
            ref.par.c[3,]
          } else {
            NA
          },
        test.nexp = fit.par.nexp,
        test.m1 = unlist(lapply(fit.par.m, function(x) x[1])),
        test.m2 = unlist(lapply(fit.par.m, function(x) x[2])),
        test.m3 = unlist(lapply(fit.par.m, function(x) x[3])),
        test.m4 = unlist(lapply(fit.par.m, function(x) x[4])),
        test.c1 = unlist(lapply(fit.par.c, function(x) x[1])),
        test.c2 = unlist(lapply(fit.par.c, function(x) x[2])),
        test.c3 = unlist(lapply(fit.par.c, function(x) x[3])),
        test.m1se = unlist(lapply(fit.par.mse, function(x) x[1])),
        test.m2se = unlist(lapply(fit.par.mse, function(x) x[2])),
        test.m3se = unlist(lapply(fit.par.mse, function(x) x[3])),
        test.m4se = unlist(lapply(fit.par.mse, function(x) x[4])),
        test.c1se = unlist(lapply(fit.par.cse, function(x) x[1])),
        test.c2se = unlist(lapply(fit.par.cse, function(x) x[2])),
        test.c3se = unlist(lapply(fit.par.cse, function(x) x[3])),
        t2.1 = get(data.names[i])[["t2"]][1,],
        t2.2 = get(data.names[i])[["t2"]][2,],
        t2.3 = get(data.names[i])[["t2"]][3,],
        t2.4 = get(data.names[i])[["t2"]][4,],
        t2.5 = get(data.names[i])[["t2"]][5,],
        t2.6 = get(data.names[i])[["t2"]][6,],
        t2.7 = get(data.names[i])[["t2"]][7,],
        t2.8 = get(data.names[i])[["t2"]][8,],
        t2.9 = get(data.names[i])[["t2"]][9,],
        t3.1 = get(data.names[i])[["t3"]][1,],
        t3.2 = get(data.names[i])[["t3"]][2,],
        t3.3 = get(data.names[i])[["t3"]][3,],
        t3.4 = get(data.names[i])[["t3"]][4,],
        t3.5 = get(data.names[i])[["t3"]][5,],
        t3.6 = get(data.names[i])[["t3"]][6,],
        t3.7 = get(data.names[i])[["t3"]][7,],
        t3.8 = get(data.names[i])[["t3"]][8,],
        t3.9 = get(data.names[i])[["t3"]][9,],
        tse = unlist(fit.par.tse)
      )
      d.melt$prop <- with(d.melt, test/ref)
      r.out <- rbind(r.out,
        ddply(d.melt, .(type), function(x) {
          c(data = data.names[i], metric = slot.names[j], sumfuncBOX(x$prop)
          )
        })
      )
      d.mid <- ddply(d.melt, .(type), function(x) {
        q75 <- as.numeric(r.out[r.out$type == unique(x$type), ]$q75)[j]
        q25 <- as.numeric(r.out[r.out$type == unique(x$type), ]$q25)[j]
        lower <- q25 - (q75 - q25)
        upper <- q75 + (q75 - q25)
        x$outlier <- ifelse(x$prop < lower | x$prop > upper, T, F)
        x$bioq <- ifelse(x$prop < 0.8 | x$prop > 1.25, F, T)
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
  write.csv(broad.out, "broad-output.csv")

# Plot
  plotdata <- plotdata[plotdata$type != "optc", ]
  plotdata$typef <- as.factor(plotdata$type)
  levels(plotdata$typef) <- c("Basic", "Optint", "Optint w/ Cmax")
  plotauc <- plotdata[plotdata$metric == "auc", ]
  plotcmax <- plotdata[plotdata$metric == "cmax" & plotdata$data %in% c("d1a", "d2a", "d3a"), ]
  plottmax <- plotdata[plotdata$metric == "tmax" & plotdata$data %in% c("d1a", "d2a", "d3a"), ]

  p1 <- ggplot(plotauc, aes(x = typef, y = prop))
  p1 <- p1 + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
  p1 <- p1 + geom_boxplot()
  p1 <- p1 + ggtitle("Broad Study AUC")
  p1 <- p1 + xlab("\nMethod")
  p1 <- p1 + ylab("Method/Reference Ratio\n")
  p1 <- p1 + facet_wrap(~data, nrow = 1)
  p1

  p2 <- ggplot(plotcmax, aes(x = typef, y = prop))
  p2 <- p2 + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
  p2 <- p2 + geom_boxplot()
  p2 <- p2 + ggtitle("Broad Study Cmax")
  p2 <- p2 + xlab("\nMethod")
  p2 <- p2 + ylab("Method/Reference Ratio\n")
  p2 <- p2 + facet_wrap(~data, nrow = 1)
  p2

  p3 <- ggplot(plottmax, aes(x = typef, y = prop))
  p3 <- p3 + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
  p3 <- p3 + geom_boxplot()
  p3 <- p3 + ggtitle("Broad Study Tmax")
  p3 <- p3 + xlab("\nMethod")
  p3 <- p3 + ylab("Method/Reference Ratio\n")
  p3 <- p3 + facet_wrap(~data, nrow = 1)
  p3
