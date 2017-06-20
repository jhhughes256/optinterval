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
  d1a <- readRDS(paste(git.dir, reponame, "d1a-broad-verbose.rds", sep = "/"))

# -----------------------------------------------------------------------------
# Data structure
  data.names <- c("d1a")
  slot.names <- c("auc", "cmax", "tmax")

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
        id = rep(1:100, 3),
        data = data.names[i],
        metric = slot.names[j],
        ref = rep(d.in$true, 3),
        test = c(d.in$basic,
          d.in$optint,
          d.in$optintwCmax),
        type = c(rep("bas", 100), rep("opt", 100), rep("optc", 100)),
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
# Now to determine if the outliers are different in some way!
  d.outlier <- plotdata[plotdata$outlier & !plotdata$bioq, ]
  time.samp <- seq(0, 24, by = 0.1)
  plot.rdata <- function(ref, test, t, n, interv, log = F) {
    plotdata <- data.frame(
      id = rep(1:n, each = length(t)),
      time = rep(t, times = n),
      dv = as.vector(ref),
      pred = as.vector(test)
    )
    xlim <- c(t[1], t[length(t)])
    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
    plotobj <- plotobj + geom_line(aes(x = time, y = dv), colour = "red")
    plotobj <- plotobj + geom_line(aes(x = time, y = pred), colour = "blue", alpha = 0.5)
    plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
    if (!log) plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
    else plotobj <- plotobj + scale_y_log10("Concentration (mg/mL)\n")
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
    plotobj <- plotobj + facet_wrap(~id, ncol = round(sqrt(n)), scales = "free")
    return(plotobj)
  }
# Set up for 2 exponential Absorption
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d.out2 <- d.outlier[with(d.outlier, test.nexp == 2 & metric == "auc" & type == "opt"), ]
  # View(d.out2[c("id", paste0("ref.", c("m1", "m2", "c")), paste0("test.", c("m1", "m2", "c1")), "prop")])
  m.out2.ref <- t(as.matrix(ddply(d.out2, .(id), function(x) {
    with(x, c(ref.m1, ref.m2, ref.c))
  })[, -1]))
  m.out2.test <- t(as.matrix(ddply(d.out2, .(id), function(x) {
    with(x, c(test.m1, test.m2, test.c1))
  })[, -1]))
  d1a1.ref <- apply(m.out2.ref, 2, function(p, x) pred.d1a(x, p), x = time.samp)
  d1a1.test <- apply(m.out2.test, 2, function(p, x) pred.d1a(x, p), x = time.samp)
  plot.rdata(d1a1.ref, d1a1.test, time.samp, dim(m.out2.test)[2], -1, log = F)

# Set up for 3 exponential absorption
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  d.out3 <- d.outlier[with(d.outlier, test.nexp == 3 & metric == "auc" & type == "opt"), ]
  # View(d.out3[c("id", paste0("ref.", c("m1", "m2", "c")), paste0("test.", c("m1", "m2", "c1")), "prop")])
  m.out3.ref <- t(as.matrix(ddply(d.out3, .(id), function(x) {
    with(x, c(ref.m1, ref.m2, ref.c))
  })[, -1]))
  m.out3.test <- t(as.matrix(ddply(d.out3, .(id), function(x) {
    with(x, c(test.m1, test.m2, test.m3, test.c1, test.c2))
  })[, -1]))
  d1a2.ref <- apply(m.out3.ref, 2, function(p, x) pred.d1a(x, p), x = time.samp)
  d1a2.test <- apply(m.out3.test, 2, function(p, x) pred.d2a(x, p), x = time.samp)
  plot.rdata(d1a2.ref, d1a2.test, time.samp, dim(m.out3.test)[2], -1, log = F)

# set up for 4 exponential absorption
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }
  d.out4 <- d.outlier[with(d.outlier, test.nexp == 4 & metric == "auc" & type == "opt"), ]
  # View(d.out4[c("id", paste0("ref.", c("m1", "m2", "c")), paste0("test.", c("m1", "m2", "c1")), "prop")])
  m.out4.ref <- t(as.matrix(ddply(d.out4, .(id), function(x) {
    with(x, c(ref.m1, ref.m2, ref.c))
  })[, -1]))
  m.out4.test <- t(as.matrix(ddply(d.out4, .(id), function(x) {
    with(x, c(test.m1, test.m2, test.c1))
  })[, -1]))
  d1a3.ref <- apply(m.out4.ref, 2, function(p, x) pred.d1a(x, p), x = time.samp)
  d1a3.test <- apply(m.out4.test, 2, function(p, x) pred.d3a(x, p), x = time.samp)
  plot.rdata(d1a3.ref, d1a3.test, time.samp, dim(m.out4.test)[2], -1, log = F)
