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
  d2b <- readRDS(paste(git.dir, reponame, "d2b-broad.rds", sep = "/"))
  d3b <- readRDS(paste(git.dir, reponame, "d3b-broad.rds", sep = "/"))
  d1a <- readRDS(paste(git.dir, reponame, "d1a-broad.rds", sep = "/"))
  d2a <- readRDS(paste(git.dir, reponame, "d2a-broad.rds", sep = "/"))
  d3a <- readRDS(paste(git.dir, reponame, "d3a-broad.rds", sep = "/"))

# -----------------------------------------------------------------------------
# Data structure
# d2b[[1]], d3b[[2]], d1a[[3]], d2a[[4]], d3a[[5]]
# names(d[[#]]) -> "data" "result"
# names(d[[#]]$result) -> "auc" "cmax" "tmax"

  data.names <- c("d2b", "d3b", "d1a", "d2a", "d3a")
  slot.names <- c("auc", "cmax", "tmax")

  res <- data.frame(NULL)
  plotdata <- data.frame(NULL)
  for (i in 1:5) {
    r.out <- data.frame(NULL)
    d.out <- data.frame(NULL)
    for (j in 1:3) {
      d.in <- get(data.names[i])[[i]]$result[slot.names[j]]
      d.melt <- data.frame(
        data = data.names[i],
        metric = slot.names[j],
        ref = rep(d.in[[1]]$true, 3),
        test = c(d.in[[1]]$basic,
          d.in[[1]]$optint,
          d.in[[1]]$optintwCmax),
        type = c(rep("bas", 1000), rep("opt", 1000), rep("optc", 1000)))
      d.melt$prop <- with(d.melt, test/ref)
      d.out <- rbind(d.out, d.melt)
      r.out <- rbind(r.out, ddply(d.melt, .(type), function(x) {
        c(data = data.names[i], metric = slot.names[j], sumfuncBOX(x$prop))
      }))
    }
    res <- rbind(res, r.out)
    plotdata <- rbind(plotdata, d.out)
  }

# -----------------------------------------------------------------------------

  box.plot.fn <- function(metric, data, x, zoom, layout = NULL) {
    subplot <- plotdata[plotdata$metric == metric & plotdata$data == data, ]
    subplot$type <- as.factor(subplot$type)
    levels(subplot$type) <- c("Basic", "Optint", "Optint w/ Cmax")

    plotobj <- ggplot(subplot, aes(x = type, y = prop))
    plotobj <- plotobj + geom_hline(yintercept = c(0.8, 1.25), color = "green4", linetype = "dashed")
    plotobj <- plotobj + geom_boxplot()
    plotobj <- plotobj + ggtitle(paste("Data:", data, "- Metric:", metric))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + ylab("Method/Reference Ratio\n")
    if (zoom) {
      ylim.box <- boxplot.stats(subplot$prop)$stats[c(1, 5)]
      # ylim.box <- c(0.6, 1.7) # "broad_boxplot_auc_zoom.png"
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
    plotobj <- plotobj + geom_pointrange(aes(ymin = q25, ymax = q75))
    plotobj <- plotobj + geom_hline(yintercept = 1, lty = 2, color = "green4")
    plotobj <- plotobj + ggtitle(paste("Data:", data, "- Metric:", metric))
    plotobj <- plotobj + xlab("\nMethod")
    plotobj <- plotobj + scale_y_continuous("Method/Reference Ratio\n")
    plotobj <- plotobj + coord_flip()
    if (is.null(layout)) print(plotobj)
    else print(plotobj, vp = layout)
  }

# -----------------------------------------------------------------------------

  # png("broad_boxplot_auc_zoom.png", width = 1200, height = 480)
  # png("broad_boxplot_auc_nozoom.png", width = 1200, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,5)))

  box.plot.fn("auc", "d1a", F, T, vplayout(1,1))
  box.plot.fn("auc", "d2a", F, T, vplayout(1,2))
  box.plot.fn("auc", "d3a", F, T, vplayout(1,3))
  box.plot.fn("auc", "d2b", F, T, vplayout(1,4))
  box.plot.fn("auc", "d3b", F, T, vplayout(1,5))

  dev.off()

  # png("broad_boxplot_cmax_nozoom.png", width = 720, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("cmax", "d1a", F, F, vplayout(1,1))
  box.plot.fn("cmax", "d2a", F, F, vplayout(1,2))
  box.plot.fn("cmax", "d3a", F, F, vplayout(1,3))

  dev.off()

  # png("broad_boxplot_tmax_nozoom.png", width = 720, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("tmax", "d1a", F, F, vplayout(1,1))
  box.plot.fn("tmax", "d2a", F, F, vplayout(1,2))
  box.plot.fn("tmax", "d3a", F, F, vplayout(1,3))

  dev.off()

# -----------------------------------------------------------------------------

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("auc", "d1a", F, T, vplayout(1,1))
  box.plot.fn("cmax", "d1a", F, T, vplayout(1,2))
  box.plot.fn("tmax", "d1a", F, T, vplayout(1,3))

  dev.off()

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))

  forest.plot.fn("auc", "d1a", vplayout(1,1))
  forest.plot.fn("cmax", "d1a", vplayout(2,1))
  forest.plot.fn("tmax", "d1a", vplayout(3,1))

  dev.off()

# -----------------------------------------------------------------------------

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("auc", "d2a", F, F, vplayout(1,1))
  box.plot.fn("cmax", "d2a", F, F, vplayout(1,2))
  box.plot.fn("tmax", "d2a", F, F, vplayout(1,3))

  dev.off()

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))

  forest.plot.fn("auc", "d2a", vplayout(1,1))
  forest.plot.fn("cmax", "d2a", vplayout(2,1))
  forest.plot.fn("tmax", "d2a", vplayout(3,1))

  dev.off()

# -----------------------------------------------------------------------------

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("auc", "d3a", F, T, vplayout(1,1))
  box.plot.fn("cmax", "d3a", F, T, vplayout(1,2))
  box.plot.fn("tmax", "d3a", F, T, vplayout(1,3))

  dev.off()

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))

  forest.plot.fn("auc", "d3a", vplayout(1,1))
  forest.plot.fn("cmax", "d3a", vplayout(2,1))
  forest.plot.fn("tmax", "d3a", vplayout(3,1))

  dev.off()

# -----------------------------------------------------------------------------

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("auc", "d2b", F, T, vplayout(1,1))
  box.plot.fn("cmax", "d2b", F, T, vplayout(1,2))
  box.plot.fn("tmax", "d2b", F, T, vplayout(1,3))

  dev.off()

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))

  forest.plot.fn("auc", "d2b", vplayout(1,1))
  forest.plot.fn("cmax", "d2b", vplayout(2,1))
  forest.plot.fn("tmax", "d2b", vplayout(3,1))

  dev.off()

# -----------------------------------------------------------------------------

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))

  box.plot.fn("auc", "d3b", F, T, vplayout(1,1))
  box.plot.fn("cmax", "d3b", F, T, vplayout(1,2))
  box.plot.fn("tmax", "d3b", F, T, vplayout(1,3))

  dev.off()

  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1)))

  forest.plot.fn("auc", "d2b", vplayout(1,1))
  forest.plot.fn("cmax", "d2b", vplayout(2,1))
  forest.plot.fn("tmax", "d2b", vplayout(3,1))

  dev.off()
