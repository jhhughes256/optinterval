# Create a function that is able to parse a matrix containing alot of data
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    graphics.off()
    if (getwd() == "C:/Users/Jim Hughes/Documents") {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "optinterval"
    } else if (getwd() == "C:/Users/hugjh001/Documents") {
      git.dir <- getwd()
      reponame <- "optinterval"
    } else if (getwd() == "C:/Users/hugjh001/Desktop") {
      git.dir <- "E:/Hughes/Git"
      reponanme <- "splines"
    }
  }

# Load packages
  library(GA)
  #library(ggplot2)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

# Source scripts and r objects to set up environment
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))

# Sourc
  d2b <- readRDS(paste(git.dir, reponame, "d2b-broad.rds", sep = "/"))
  d3b <- readRDS(paste(git.dir, reponame, "d3b-broad.rds", sep = "/"))
  d1a <- readRDS(paste(git.dir, reponame, "d1a-broad.rds", sep = "/"))
  #d2a <- readRDS(paste(git.dir, reponame, "d2a-broad.rds", sep = "/"))
  #d3a <- readRDS(paste(git.dir, reponame, "d3a-broad.rds3, sep = "/"))

# -----------------------------------------------------------------------------
