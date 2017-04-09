# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

# Load libraries
  library(splines)
  library(ggplot2)
  library(plyr)
  library(stringr)

# Set working directory and load functions
  #setwd("C:/Users/hugjh001/Documents/optinterval")
  setwd("C:/Users/Jim Hughes/Documents/GitRepos/optinterval")
  source("spline_functions_ns.R")
  source("pk_data.R")

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

#
  which.knots <- function(data, max.knots = NULL) {
  # Setup environment
    if (is.null(max.knots)) max.knots <- 3
    x <- data[,1]
    y <- data[,2]
    err <- NULL
    knots <- list(NULL)
  # Determine initial estimate using tmax (xmax)
    xmax <- x[which(y == max(y))]
  # Run optim for each number of knots
    for (i in 1:max.knots) {
      par <- rep(xmax, i)
    # Determine upper and lower boundaries
      #low.lim <- c(min(x), rep(1, i - 1))
      #upp.lim <- rep(max(x), i)
    # Optimise
      opt.res <- optim(
        par,  # Initial parameter estimates
        mlespline.knotfit,  # Fitting function
        method = "L-BFGS-B", lower = min(x) + 0.001, upper = max(x),
        data = data  # Function arguments
      )
    # Collate results
      err[i] <- opt.res$value
      knots[[i]] <- knot.func(opt.res$par)
    }
  # Organise output
    list(
      error = err,
      knots = knots
    )
  }

# -----------------------------------------------------------------------------
# Set the possible sampling times
# These sample.times will be used by indexing 1:(number of obs)
# This way we have the important sampling points always and with each increase
# in number of obs, the observations that are added are relevant
# These times are ideally put in ascending order before use
  sample.times <- c(
    0, 0.5, 1, 2, 4, 8, 24, 0.25, 12, 3, 6,
    18, 10, 1.5, 0.75, 5, 7, 14, 16, 20, 22
  )

# Sample from the simulated data
# pkdata1[pkdata1$x %in% sample.times[1:7], ]
# Then add error on using an error model as seen in mle script


  spline.to.plot <- function(data, spline) {
    plot.df <- ldply(seq_len(length(spline)), function(n) {
      knots <- spline[[n]]$knots
      mle.vec <- spline[[n]]$error
      ldply(knots, function(K) {
        param <- list(
          X = data$time,
          Y = data$conc,
          K = K,
          n = n
        )
        K.len <- length(K)
        K.vec <- c(K, rep(NULL, length(knots) - K.len))
        mle <- signif(mle.vec[[K.len]], 3)
        out <- data.frame(
          time = seq(from = 0, to = 24, by = 0.25),
          conc = predict.spline.func(param, times),
          degree = n,
          nknots = length(K),
          spline.err = mle
        )
        err <- conc - out$conc
        sigma <- sd(err)
        loglik <- dnorm(conc, out$conc, sigma, log = T)
        out$error <- signif(-1*sum(loglik), 3)
        for (i in 1:length(knots)) {
          out[[paste0("knot", i)]] <- K.vec[i]
        }
        out
      })
    })  # ldply
    plot.df$knots <- unlist(
      dlply(plot.df, .(degree, nknots), function(x) {
        K <- unique(x$nknots)
        df.knots <- NULL
        for (i in 1:K) {
          df.knots <- paste(df.knots, round(x[[paste0("knot", i)]], 1))
        }
        df.knots
      })
    )
    plot.df$facet <- factor(
      paste0("Degree: ", plot.df$degree,
        "   Knots: ", plot.df$nknots, "\n",
        plot.df$knots
      )
    )
    plot.df
  }  # spline.to.plot

# -----------------------------------------------------------------------------
# Naturally a 6 degree spline fits the data the best but we aren't interested in
# the sse for spline against the observed, but infact the spline against the
# truth!
  plotdata <- spline.to.plot(onecomp, spline.data)
  truedata <- data.frame(
    times = rep(times, times = length(levels(plotdata$facet))),
    conc = rep(conc, times = length(levels(plotdata$facet))),
    facet = factor(rep(levels(plotdata$facet), each = length(times)))
  )

  plotobj <- NULL
  plotobj <- ggplot()
  plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = onecomp)
  plotobj <- plotobj + geom_line(aes(x = time, y = conc), data = plotdata, size = 1, colour = "red")
  plotobj <- plotobj + geom_line(aes(x = times, y = conc), data = truedata, size = 0.5, linetype = 2)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.6, label = error), data = plotdata)
  plotobj <- plotobj + geom_text(aes(x = 20, y = 0.3, label = spline.err), data = plotdata, colour = "red")
  plotobj <- plotobj + scale_y_log10("Concentration (mg/mL) \n", lim = c(0.01, 0.8))
  plotobj <- plotobj + scale_x_continuous("\n Time (hours)", lim = c(0, 24))
  plotobj <- plotobj + facet_wrap(~facet, ncol = 4)
  plotobj
