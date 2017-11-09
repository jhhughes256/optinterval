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
  d1a1 <- readRDS(paste(git.dir, reponame, "d1a-broad-50-250.rds", sep = "/"))
  d1a2 <- readRDS(paste(git.dir, reponame, "d1a-broad-50-500.rds", sep = "/"))
  d1a3 <- readRDS(paste(git.dir, reponame, "d1a-broad-100-250.rds", sep = "/"))
  d1aj5 <- readRDS(paste(git.dir, reponame, "d1a-test-j5.rds", sep = "/"))
  d1aj6 <- readRDS(paste(git.dir, reponame, "d1a-test-j6.rds", sep = "/"))
  d1aj7 <- readRDS(paste(git.dir, reponame, "d1a-test-j7.rds", sep = "/"))

# -----------------------------------------------------------------------------
# Changes in maxiter & popSize
  # data.names <- c("d1a1", "d1a2", "d1a3")
  # iter.vals <- c(50, 50, 100)
  #

# Testing new repeated-ga algorithm
  data.names <- c("d1a1", "d1a2", "d1a3", "d1aj5", "d1aj6", "d1aj7")
  j.vals <- c(1, 1, 1, 5, 6, 7)
  iter.vals <- c(50, 50, 100, 50, 50, 50)
  popn.vals <- c(250, 500, 250, 250, 250, 250)

# Collate data
  slot.names <- c("auc", "cmax", "tmax")
  res <- data.frame(NULL)
  plotdata <- data.frame(NULL)
  for (i in 1:length(data.names)) {
    r.out <- data.frame(NULL)
    d.out <- data.frame(NULL)
    for (j in 1:length(slot.names)) {
      d.in <- get(data.names[i])[[slot.names[j]]]
      d.melt <- data.frame(
        data = data.names[i],
        iter = iter.vals[i],
        popn = popn.vals[i],
        j = j.vals[i],
        metric = slot.names[j],
        ref = rep(d.in$true, 3),
        test = c(d.in$basic,
          d.in$optint,
          d.in$optintwCmax),
        type = c(rep("bas", 1000), rep("opt", 1000), rep("optc", 1000)))
      d.melt$prop <- with(d.melt, test/ref)
      d.out <- rbind(d.out, d.melt)
      r.out <- rbind(r.out,
        ddply(d.melt, .(type), function(x) {
          c(data = data.names[i], iter = iter.vals[i], popn = popn.vals[i],
            j = j.vals[i], metric = slot.names[j], sumfuncBOX(x$prop)
          )
        })
      )
    }
    res <- rbind(res, r.out)
    plotdata <- rbind(plotdata, d.out)
  }

# -----------------------------------------------------------------------------
# Plot functions
  box.plot.fn <- function(metric, data, x, zoom, layout = NULL) {
    subplot <- plotdata[plotdata$metric == metric & plotdata$data == data, ]
    subplot$type <- as.factor(subplot$type)
    levels(subplot$type) <- c("Basic", "Optint", "Optint w/ Cmax")

    plotobj <- ggplot(subplot, aes(x = type, y = prop))
    plotobj <- plotobj + geom_hline(yintercept = 1, color = "green4", linetype = "dashed")
    plotobj <- plotobj + geom_boxplot()
    # plotobj <- plotobj + ggtitle(paste("Metric:", metric))
    plotobj <- plotobj + ggtitle(paste(metric, "j ==", unique(subplot$j)))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + ylab("Method/Reference Ratio\n")
    if (zoom) {
      # ylim.box <- boxplot.stats(subplot$prop)$stats[c(1, 5)]
      ylim.box <- c(0.2, 3)
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
  pushViewport(viewport(layout = grid.layout(1,3)))

  # box.plot.fn("tmax", "d1a1", F, F, vplayout(1,1))
  # box.plot.fn("tmax", "d1a2", F, F, vplayout(1,2))
  # box.plot.fn("tmax", "d1a3", F, F, vplayout(1,3))
  box.plot.fn("auc", "d1a1", F, T, vplayout(1,1))
  box.plot.fn("auc", "d1a2", F, T, vplayout(1,2))
  box.plot.fn("auc", "d1a3", F, T, vplayout(1,3))

  dev.off()

  png("broad_new_forestplot_tmax.png", width = 720, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))

  # forest.plot.fn("tmax", "d1a1", vplayout(1,1))
  # forest.plot.fn("tmax", "d1a2", vplayout(2,1))
  # forest.plot.fn("tmax", "d1a3", vplayout(3,1))
  forest.plot.fn("tmax", "d1a1", vplayout(1,1))
  forest.plot.fn("tmax", "d1a2", vplayout(2,1))
  forest.plot.fn("tmax", "d1a3", vplayout(3,1))

  dev.off()
