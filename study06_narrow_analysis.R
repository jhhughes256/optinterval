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
  #library(GA)
  library(ggplot2)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  library(plyr)
  library(grid)

# Source scripts and r objects to set up environment
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "sumstat_functions.R", sep = "/"))

# Source data
  d1a <- readRDS(paste(git.dir, reponame, "d1a-narrow.rds", sep = "/"))
  d1a.user <- readRDS(paste(git.dir, reponame, "d1a-narrow-user.rds", sep = "/"))

  user.auc <- unique(with(d1a.user$auc, user/true))
  user.cmax <- unique(with(d1a.user$cmax, user/true))
  user.tmax <- unique(with(d1a.user$tmax, user/true))
  c(user.auc, user.cmax, user.tmax)

# -----------------------------------------------------------------------------
# Data structure
# > names(d1a)
# [1] "auc"  "cmax" "tmax"

  slot.names <- c("auc", "cmax", "tmax")
  res <- data.frame(NULL)
  plotdata <- data.frame(NULL)
  for (i in 1:3) {  # 1:5) {
    d.in <- d1a[[slot.names[i]]]
    d.melt <- data.frame(
      metric = slot.names[i],
      ref = rep(d.in$true, 3),
      test = c(d.in$basic,
        d.in$optint,
        d.in$optintwCmax),
      type = c(rep("bas", 1000), rep("opt", 1000), rep("optc", 1000)))
    d.melt$prop <- with(d.melt, test/ref)
    res <- rbind(res, ddply(d.melt, .(type), function(x) {
      c(metric = slot.names[i], sumfuncBOX(x$prop))
    }))
    plotdata <- rbind(plotdata, d.melt)
  }

# -----------------------------------------------------------------------------
  box.plot.fn <- function(metric, user, x, zoom, layout = NULL) {
    subplot <- plotdata[plotdata$metric == metric & plotdata$type != "bas", ]
    subplot$type <- as.factor(subplot$type)
    levels(subplot$type) <- c("Basic", "Optint", "Optint w/ Cmax")

    plotobj <- ggplot(subplot, aes(x = type, y = prop))
    plotobj <- plotobj + geom_boxplot()
    plotobj <- plotobj + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
    plotobj <- plotobj + geom_hline(yintercept = user, color = "red", linetype = "dashed")
    plotobj <- plotobj + ggtitle(paste("Metric:", metric))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + ylab("Method/Reference Ratio\n")
    if (zoom) {
      lim.check <- c(boxplot.stats(subplot$prop)$stats[c(1, 5)], user)
      ylim.box <- c(min(lim.check), max(lim.check))
      plotobj <- plotobj + coord_cartesian(ylim = ylim.box)
    }
    if (x) {
      plotobj <- plotobj + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 4)
    }
    if (is.null(layout)) print(plotobj)
    else print(plotobj, vp = layout)
  }

  forest.plot.fn <- function(metric, user, zoom, layout = NULL) {
    subplot <- res[res$metric == metric & res$type != "bas", ]
    subplot$type <- factor(subplot$type, levels = rev(subplot$type))
    levels(subplot$type) <- c("Optint w/ Cmax", "Optint", "Basic")
    subplot$median <- as.numeric(subplot$median)
    subplot$q025 <- as.numeric(subplot$q025)
    subplot$q25 <- as.numeric(subplot$q25)
    subplot$q75 <- as.numeric(subplot$q75)
    subplot$q975 <- as.numeric(subplot$q975)

    plotobj <- ggplot(subplot, aes(x = type, y = median))
    plotobj <- plotobj + geom_linerange(aes(ymin = q025, ymax = q975), color = "red")
    plotobj <- plotobj + geom_pointrange(aes(ymin = q25, ymax = q75))
    plotobj <- plotobj + geom_hline(yintercept = c(0.8, 1.25), lty = 2, color = "green4")
    plotobj <- plotobj + geom_hline(yintercept = user, color = "red", linetype = "dashed")
    plotobj <- plotobj + ggtitle(paste("Metric:", metric))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + scale_y_continuous("Method/Reference Ratio\n")
    if (zoom) {
      lim.check <- c(boxplot.stats(subplot$prop)$stats[c(1, 5)], user)
      ylim.box <- c(min(lim.check), max(lim.check))
      plotobj <- plotobj + coord_cartesian(ylim = ylim.box)
    }
    plotobj <- plotobj + coord_flip()
    if (is.null(layout)) print(plotobj)
    else print(plotobj, vp = layout)
  }

# -----------------------------------------------------------------------------

  box.plot.fn("auc", user.auc, F, F)
  box.plot.fn("cmax", user.cmax, F, T)
  box.plot.fn("tmax", user.tmax, F, T)

  forest.plot.fn("auc", user.auc, F)
  # ggsave("narrow_forest_auc.png")
  forest.plot.fn("cmax", user.cmax, F)
  # ggsave("narrow_forest_cmax.png")
  forest.plot.fn("tmax", user.tmax, F)
  # ggsave("narrow_forest_tmax.png")
