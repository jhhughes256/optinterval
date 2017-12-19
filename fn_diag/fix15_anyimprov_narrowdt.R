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
  library(grid)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Source functions
  source(paste(git.dir, reponame, "sumstat_functions.R", sep = "/"))

# Source files of interest
  sumexp.string <- "d2a"
  ref.string <- "dt05"
  test.string <- "dt15"
  ref <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/", sumexp.string, "-narrow-", ref.string, ".rds"), sep = "/"))
  test <- readRDS(paste(git.dir, reponame,
    paste0("fn_diag/", sumexp.string, "-narrow-", test.string, ".rds"), sep = "/"))

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
    ref.par.m <- get(data.names[i])[["par"]][1:ceiling(length(get(data.names[i])[["par"]])/2)]
    ref.par.c <- get(data.names[i])[["par"]][-(1:ceiling(length(get(data.names[i])[["par"]])/2))]
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
    }
    for (j in 1:3) {
      d.in <- get(data.names[i])[[slot.names[j]]]
      d.melt <- data.frame(
        id = rep(1:niter, 3),
        data = data.names[i],
        metric = slot.names[j],
        ref = rep(d.in$true, 3),
        test = c(d.in$user,
          d.in$optint,
          d.in$optintwCmax),
        type = c(rep("bas", niter), rep("opt", niter), rep("optc", niter)),
        ref.m1 = ref.par.m[1],
        ref.m2 = ref.par.m[2],
        ref.m3 = if (length(ref.par.m) > 2) {
            ref.par.m[3]
          } else {
            NA
          },
        ref.m4 = if (length(ref.par.m) > 3) {
            ref.par.m[4]
          } else {
            NA
          },
        ref.c1 = ref.par.c[1],
        ref.c2 = if (length(ref.par.c) > 1) {
            ref.par.c[2]
          } else {
            NA
          },
        ref.c3 = if (length(ref.par.c) > 2) {
            ref.par.c[3]
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
        # test.m1se = unlist(lapply(fit.par.mse, function(x) x[1])),
        # test.m2se = unlist(lapply(fit.par.mse, function(x) x[2])),
        # test.m3se = unlist(lapply(fit.par.mse, function(x) x[3])),
        # test.m4se = unlist(lapply(fit.par.mse, function(x) x[4])),
        # test.c1se = unlist(lapply(fit.par.cse, function(x) x[1])),
        # test.c2se = unlist(lapply(fit.par.cse, function(x) x[2])),
        # test.c3se = unlist(lapply(fit.par.cse, function(x) x[3])),
        t2.1 = get(data.names[i])[["t2"]][2,],
        t2.2 = get(data.names[i])[["t2"]][3,],
        t2.3 = get(data.names[i])[["t2"]][4,],
        t2.4 = get(data.names[i])[["t2"]][5,],
        t2.5 = get(data.names[i])[["t2"]][6,],
        t2.6 = get(data.names[i])[["t2"]][7,],
        t2.7 = get(data.names[i])[["t2"]][8,],
        t2.8 = get(data.names[i])[["t2"]][9,],
        t2.9 = get(data.names[i])[["t2"]][10,],
        t3.1 = get(data.names[i])[["t3"]][1,],
        t3.2 = get(data.names[i])[["t3"]][2,],
        t3.3 = get(data.names[i])[["t3"]][3,],
        t3.4 = get(data.names[i])[["t3"]][4,],
        t3.5 = get(data.names[i])[["t3"]][5,],
        t3.6 = get(data.names[i])[["t3"]][6,],
        t3.7 = get(data.names[i])[["t3"]][7,],
        t3.8 = get(data.names[i])[["t3"]][8,],
        t3.9 = get(data.names[i])[["t3"]][9,]#,
        # tse = unlist(fit.par.tse)
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
# Plots

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
      ylim.box <- boxplot.stats(subplot$prop)$stats[c(1, 5)]
      # ylim.box <- c(0.2, 1.5)
      plotobj <- plotobj + coord_cartesian(ylim = ylim.box)
    }
    if (x) plotobj <- plotobj + stat_summary(fun.y = mean, geom = "point", shape = 4, size = 4)
    if (is.null(layout)) print(plotobj)
    else print(plotobj, vp = layout)
  }

# -----------------------------------------------------------------------------
# Plot data
  # png("broad_new_boxplot_auc.png", width = 360, height = 480)
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))

  box.plot.fn("auc", "ref", F, F, vplayout(1,1))
  box.plot.fn("auc", "test", F, F, vplayout(1,2))

  box.plot.fn("cmax", "ref", F, F, vplayout(1,1))
  box.plot.fn("cmax", "test", F, F, vplayout(1,2))

  box.plot.fn("tmax", "ref", F, F, vplayout(1,1))
  box.plot.fn("tmax", "test", F, F, vplayout(1,2))

  dev.off()

# -----------------------------------------------------------------------------
  pred.sumexp <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    m <- x[1:n]
    ord <- order(m, decreasing = T)
    p <- c(m[ord], x[(n+1):l])
    for (i in 1:n) {
      if (i == 1) y <- p[i]^d*exp(p[i]*t + p[n+i])
      else if (i != n | a == 0) y <- y + p[i]^d*exp(p[i]*t + p[n+i])
      else if (a == 1) y <- y - p[i]^d*exp(p[i]*t)*sum(exp(p[(n+1):(2*n-1)]))
    }
    return(y)
  }


  subdata <- plotdata[plotdata$type == "opt" & plotdata$metric == "auc",]
  timedata <- data.frame(mean = NULL, CI90lo = NULL, CI90hi = NULL)

  subpar <- as.numeric(na.omit(as.numeric(subdata[1, 7:13])))
  subpartimes <- seq(0, 24, by = 0.1)
  subpardata <- data.frame(
    time = subpartimes,
    dv = pred.sumexp(subpar, subpartimes)
  )

  timedata <- ddply(subdata, .(data), function(x) {
    meantime <- c(0, mean(x$t2.2), mean(x$t2.3), mean(x$t2.4), mean(x$t2.5), mean(x$t2.6), mean(x$t2.7), mean(x$t2.8), 24)
    ci90lotime <- c(
      quantile(x$t2.2, probs = 0.05, na.rm = T, names = F),
      quantile(x$t2.3, probs = 0.05, na.rm = T, names = F),
      quantile(x$t2.4, probs = 0.05, na.rm = T, names = F),
      quantile(x$t2.5, probs = 0.05, na.rm = T, names = F),
      quantile(x$t2.6, probs = 0.05, na.rm = T, names = F),
      quantile(x$t2.7, probs = 0.05, na.rm = T, names = F),
      quantile(x$t2.8, probs = 0.05, na.rm = T, names = F)
    )
    ci90hitime <- c(
      quantile(x$t2.2, probs = 0.95, na.rm = T, names = F),
      quantile(x$t2.3, probs = 0.95, na.rm = T, names = F),
      quantile(x$t2.4, probs = 0.95, na.rm = T, names = F),
      quantile(x$t2.5, probs = 0.95, na.rm = T, names = F),
      quantile(x$t2.6, probs = 0.95, na.rm = T, names = F),
      quantile(x$t2.7, probs = 0.95, na.rm = T, names = F),
      quantile(x$t2.8, probs = 0.95, na.rm = T, names = F)
    )
    data.frame(
      time = meantime,
      dv = pred.sumexp(subpar, meantime),
      ci90lo = c(NA, ci90lotime, NA),
      ci90hi = c(NA, ci90hitime, NA)
    )
  })

  p <- NULL
  p <- ggplot()
  p <- p + geom_line(aes(x = time, y = dv), data = subpardata, colour = "blue", size = 0.8)
  p <- p + geom_line(aes(x = time, y = dv), data = timedata, colour = "green4", linetype = "dashed", size = 0.8)
  p <- p + geom_point(aes(x = time, y = dv), data = timedata, size = 1.5, colour = "red")
  p <- p + geom_errorbarh(aes(y = dv, xmin = ci90lo, xmax = ci90hi, x = time), data = timedata, colour = "black", size = 1, height = 0.6)
  p <- p + facet_wrap(~data, ncol = 1)
  p
