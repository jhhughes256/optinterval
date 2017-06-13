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
  d2b <- readRDS(paste(git.dir, reponame, "d2b-broad.rds", sep = "/"))
  d3b <- readRDS(paste(git.dir, reponame, "d3b-broad.rds", sep = "/"))
  d1a <- readRDS(paste(git.dir, reponame, "d1a-broad.rds", sep = "/"))
  #d2a <- readRDS(paste(git.dir, reponame, "d2a-broad.rds", sep = "/"))
  #d3a <- readRDS(paste(git.dir, reponame, "d3a-broad.rds3, sep = "/"))

# -----------------------------------------------------------------------------
# Data structure
# d2b[[1]], d3b[[2]], d1a[[3]], d2a[[4]], d3a[[5]]
# names(d[[#]]) -> "data" "result"
# names(d[[#]]$result) -> "auc" "cmax" "tmax"

  data.names <- c("d2b", "d3b", "d1a")  #, "d2a", "d3a")
  slot.names <- c("auc", "cmax", "tmax")

  res <- data.frame(NULL)
  plotdata <- data.frame(NULL)
  for (i in 1:3) {  # 1:5) {
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

  subplot <- plotdata[plotdata$metric == "auc", ]
  subplot$data <- as.factor(subplot$data)
  subplot$type <- as.factor(subplot$type)

  plotobj <- ggplot(subplot, aes(x = type, y = prop))
  plotobj <- plotobj + geom_boxplot()
  plotobj <- plotobj + facet_wrap(~data)
  plotobj

  plotobj + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 4)
