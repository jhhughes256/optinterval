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

# Source scripts and r objects to set up environment
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "sumstat_functions.R", sep = "/"))

# Source data
  d1a <- readRDS(paste(git.dir, reponame, "d1a-narrow.rds", sep = "/"))

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
  res.plot.fn <- function(metric, user, plot) {
    subplot <- plotdata[plotdata$metric == metric & plotdata$type != "bas", ]
    subplot$type <- as.factor(subplot$type)

    plotobj <- ggplot(subplot, aes(x = type, y = prop))
    plotobj <- plotobj + geom_boxplot()
    plotobj <- plotobj + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
    plotobj <- plotobj + geom_hline(yintercept = user, color = "red", linetype = "dashed")

    if (plot == 1) {
      lim.check <- c(boxplot.stats(subplot$prop)$stats[c(1, 5)], user)
      ylim.box <- c(min(lim.check), max(lim.check))
      plotobj <- plotobj + coord_cartesian(ylim = ylim.box)
      plotobj
    } else if (plot == 2) {
      plotobj <- plotobj + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 4)
      plotobj
    } else {
      plotobj
    }
  }

  res.plot.fn("auc", 1.016367, 0)
  res.plot.fn("auc", 1.016367, 1)
  res.plot.fn("auc", 1.016367, 2)
  res.plot.fn("cmax", 0.9905159, 0)
  res.plot.fn("cmax", 0.9905159, 1)
  res.plot.fn("cmax", 0.9905159, 2)
  res.plot.fn("tmax", 0.8571429, 0)
  res.plot.fn("tmax", 0.8571429, 1)
  res.plot.fn("tmax", 0.8571429, 2)
