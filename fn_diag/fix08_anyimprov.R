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
  sumexp.string <- "d2a"
  ref.string <- "verbose"
  test.string <- "popsize"
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
        t3.9 = get(data.names[i])[["t3"]][9,]
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
  refdata.opttmax <- refdata.opt[refdata.opt$metric == "tmax",]
  testdata.opttmax <- testdata.opt[testdata.opt$metric == "tmax",]

  length(which(refdata.optauc$outlier))
  length(which(testdata.optauc$outlier))
  length(which(refdata.opttmax$outlier))
  length(which(testdata.opttmax$outlier))

  optauc.outliers <- which(refdata.optauc$outlier)
  subref.optauc <- refdata.optauc[optauc.outliers,]
  subtest.optauc <- testdata.optauc[optauc.outliers,]

  sumfuncCV(abs(1-subref.optauc$prop)-abs(1-subtest.optauc$prop))
  sumfuncBOX(abs(1-subref.optauc$prop)-abs(1-subtest.optauc$prop))

  opttmax.outliers <- which(refdata.opttmax$outlier)
  subref.opttmax <- refdata.opttmax[opttmax.outliers,]
  subtest.opttmax <- testdata.opttmax[opttmax.outliers,]

  sumfuncCV(abs(1-subref.opttmax$prop)-abs(1-subtest.opttmax$prop))
  sumfuncBOX(abs(1-subref.opttmax$prop)-abs(1-subtest.opttmax$prop))

  refdata.optc <- refdata[refdata$type == "optc",]
  testdata.optc <- testdata[testdata$type == "optc",]
  refdata.optcauc <- refdata.optc[refdata.optc$metric == "auc",]
  testdata.optcauc <- testdata.optc[testdata.optc$metric == "auc",]
  refdata.optctmax <- refdata.optc[refdata.optc$metric == "tmax",]
  testdata.optctmax <- testdata.optc[testdata.optc$metric == "tmax",]

  length(which(refdata.optcauc$outlier))
  length(which(testdata.optcauc$outlier))
  length(which(refdata.optctmax$outlier))
  length(which(testdata.optctmax$outlier))

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
