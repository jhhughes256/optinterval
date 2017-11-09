# A script designed to determine algorithm improvement
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
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Source functions
  source(paste(git.dir, reponame, "sumstat_functions.R", sep = "/"))

# Source files of interest
  sumexp.string <- "d3a"
  ref.string <- "verbose"
  test.string <- "dt"
  ref <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/", sumexp.string, "-broad-", ref.string, ".rds"), sep = "/"))
  test <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/", sumexp.string, "-broad-", test.string, ".rds"), sep = "/"))

# -----------------------------------------------------------------------------
# Data structure
  data.names <- c("ref", "test")
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
      if (i == 1) d.in$optint <- d.in$optintwCmax
      if (i == 2) d.in$optintwCmax <- d.in$optint
      d.melt <- data.frame(
        id = rep(1:niter, 2),
        data = data.names[i],
        metric = slot.names[j],
        ref = rep(d.in$true, 2),
        test = c(d.in$basic,
          d.in$optint),
        type = c(rep("bas", niter), rep("opt", niter)),
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
#
  refdata <- plotdata[plotdata$data == "ref",]
  testdata <- plotdata[plotdata$data == "test",]
  refdata.opt <- refdata[refdata$type == "opt",]
  testdata.opt <- testdata[testdata$type == "opt",]
  refdata.optauc <- refdata.opt[refdata.opt$metric == "auc",]
  testdata.optauc <- testdata.opt[testdata.opt$metric == "auc",]

  length(which(refdata.optauc$outlier))
  length(which(testdata.optauc$outlier))

  optauc.outliers <- which(refdata.optauc$outlier)
  subref.optauc <- refdata.optauc[optauc.outliers,]
  subtest.optauc <- testdata.optauc[optauc.outliers,]

  sumfuncCV(abs(1-subref.optauc$prop)-abs(1-subtest.optauc$prop))
  sumfuncBOX(abs(1-subref.optauc$prop)-abs(1-subtest.optauc$prop))

  refdata.optc <- refdata[refdata$type == "optc",]
  testdata.optc <- testdata[testdata$type == "optc",]
  refdata.optcauc <- refdata.optc[refdata.optc$metric == "auc",]
  testdata.optcauc <- testdata.optc[testdata.optc$metric == "auc",]

  length(which(refdata.optcauc$outlier))
  length(which(testdata.optcauc$outlier))

  optcauc.outliers <- which(refdata.optcauc$outlier)
  subref.optcauc <- refdata.optcauc[optcauc.outliers,]
  subtest.optcauc <- testdata.optcauc[optcauc.outliers,]

  sumfuncCV(abs(1-subref.optcauc$prop)-abs(1-subtest.optcauc$prop))
  sumfuncBOX(abs(1-subref.optcauc$prop)-abs(1-subtest.optcauc$prop))

  optctmax.outliers <- which(refdata.optctmax$outlier)
  subref.optctmax <- refdata.optctmax[optctmax.outliers,]
  subtest.optctmax <- testdata.optctmax[optctmax.outliers,]

  sumfuncCV(abs(1-subref.optctmax$prop)-abs(1-subtest.optctmax$prop))
  sumfuncBOX(abs(1-subref.optctmax$prop)-abs(1-subtest.optctmax$prop))

  box.plot.fn <- function(metric, data, x, zoom, layout = NULL) {
    subplot <- plotdata[plotdata$metric == metric & plotdata$data == data, ]
    subplot$type <- as.factor(subplot$type)
    levels(subplot$type) <- c("Basic", "Optint", "Optint w/ Cmax")

    plotobj <- ggplot(subplot, aes(x = type, y = prop))
    plotobj <- plotobj + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
    plotobj <- plotobj + geom_boxplot()
    # plotobj <- plotobj + ggtitle(paste("Metric:", metric))
    plotobj <- plotobj + ggtitle(paste(metric, data))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + ylab("Method/Reference Ratio\n")
    if (zoom) {
      # ylim.box <- boxplot.stats(subplot$prop)$stats[c(1, 5)]
      ylim.box <- c(0.2, 1.5)
      plotobj <- plotobj + coord_cartesian(ylim = ylim.box)
    }
    if (x) plotobj <- plotobj + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 4)
    if (is.null(layout)) print(plotobj)
    else print(plotobj, vp = layout)
  }

  forest.plot.fn <- function(metric, data, layout = NULL) {
    subplot <- res[res$metric == metric & res$data == data, ]
    subplot$type <- factor(subplot$type, levels = rev(subplot$type))
    levels(subplot$type) <- c("Optint w/ Cmax", "Optint", "Basic")
    subplot$median <- as.numeric(subplot$median)
    subplot$q025 <- as.numeric(subplot$q025)
    subplot$q25 <- as.numeric(subplot$q25)
    subplot$q75 <- as.numeric(subplot$q75)
    subplot$q975 <- as.numeric(subplot$q975)

    plotobj <- ggplot(subplot, aes(x = type, y = median))
    plotobj <- plotobj + geom_linerange(aes(ymin = q025, ymax = q975), color = "red")
    plotobj <- plotobj + geom_linerange(aes(ymin = q25, ymax = q75), size = 1.2)
    plotobj <- plotobj + geom_pointrange(aes(ymin = median, ymax = median), size = 0.8)
    plotobj <- plotobj + geom_hline(yintercept = 1, lty = 2, color = "green4")
    # plotobj <- plotobj + ggtitle(paste("Metric:", metric))
    plotobj <- plotobj + ggtitle(paste(metric, "j ==", unique(subplot$j)))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + scale_y_continuous("Method/Reference Ratio\n",
      lim = c(0.2, 4)) # auc - lim = c(0.35, 2.4)), cmax - c(0.2, 1.05)), tmax - c(0.2, 4))
    plotobj <- plotobj + coord_flip()
    if (is.null(layout)) print(plotobj)
    else print(plotobj, vp = layout)
  }

# -----------------------------------------------------------------------------
# Plot data
  # png("broad_new_boxplot_auc.png", width = 720, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))

  box.plot.fn("auc", "ref", F, T, vplayout(1,1))
  box.plot.fn("auc", "test", F, T, vplayout(1,2))

  box.plot.fn("cmax", "ref", F, T, vplayout(1,1))
  box.plot.fn("cmax", "test", F, T, vplayout(1,2))

  box.plot.fn("tmax", "ref", F, F, vplayout(1,1))
  box.plot.fn("tmax", "test", F, F, vplayout(1,2))

  dev.off()

  png("broad_new_forestplot_tmax.png", width = 720, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))

  forest.plot.fn("tmax", "ref", vplayout(1,1))
  forest.plot.fn("tmax", "test", vplayout(2,1))

  dev.off()
