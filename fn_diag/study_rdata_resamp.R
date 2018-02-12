# Random sum of exponential data for sourcing
# -----------------------------------------------------------------------------
# The datasets in order are:
#   d2b - two compartment kinetics given iv
#   d3b - three compartment kinetics given iv
#   d1a - one compartment kinetics given orally
#   d2a - two compartment kinetics given orally
#   d3a - three compartment kinetics given orally
# -----------------------------------------------------------------------------
# Uncomment below code to view random data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# library(ggplot2)
# theme_bw2 <- theme_set(theme_bw(base_size = 14))
# theme_update(plot.title = element_text(hjust = 0.5))
# plot.rdata <- function(data, t, n, interv, log = F) {
#   plotdata <- data.frame(
#     id = rep(1:n, each = length(t)),
#     time = rep(t, times = n),
#     dv = as.vector(data)
#   )
#   xlim <- c(t[1], t[length(t)])
#   plotobj <- NULL
#   plotobj <- ggplot(data = plotdata)
#   plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
#   plotobj <- plotobj + geom_line(aes(x = time, y = dv), colour = "red")
#   plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
#   if (!log) plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
#   else plotobj <- plotobj + scale_y_log10("Concentration (mg/mL)\n")
#   plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
#   plotobj <- plotobj + facet_wrap(~id, ncol = round(sqrt(n)), scales = "free")
#   return(plotobj)
# }
# -----------------------------------------------------------------------------
  if (!exists("niter")) niter <- 16
  if (!exists("sdev")) sdev <- 4
  time.samp <- seq(0, 24, by = 0.1)
# Two Compartment Kinetics
  pred.d2b <- function(x, p) {
    exp(p[1]*x + p[3]) + exp(p[2]*x + p[4])
  }
  d2b.p <- matrix(nrow = 4, ncol = niter)
  d2b.sim <- 1:niter
  repeat {
    nsim <- length(d2b.sim)
    d2b.p[1, d2b.sim] <- runif(nsim, -1, -0.1)
    d2b.p[2, d2b.sim] <- runif(nsim, d2b.p[1, d2b.sim]*0.8, d2b.p[1, d2b.sim]*0.05)
    d2b.p[3, d2b.sim] <- runif(nsim, 1, 6)
    d2b.p[4, d2b.sim] <- runif(nsim, d2b.p[3, d2b.sim] - 1, d2b.p[3, d2b.sim] - 0.2)
    d2b.all <- apply(d2b.p, 2, function(p, x) pred.d2b(x, p), x = time.samp)
    d2b.tmax <- apply(d2b.all, 2, function(x) time.samp[which(x == max(x))])
    d2b.sim <- which(d2b.tmax > 12)
    if (length(d2b.sim) == 0) break
    else print(paste("d2b repeat", length(d2b.sim)))
  }
  d2b.t <- c(0, 0.5, 1, 2, 4, 8, 12, 16, 24)
  #plot.rdata(d2b, time.samp, niter, -1, log = F)

# Three Compartment Kinetics
  pred.d3b <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) + exp(p[3]*x + p[6])
  }
  d3b.p <- matrix(nrow = 6, ncol = niter)
  d3b.sim <- 1:niter
  repeat {
    nsim <- length(d3b.sim)
    d3b.p[1, d3b.sim] <- runif(nsim, -1, -0.1)
    d3b.p[2, d3b.sim] <- runif(nsim, d3b.p[1, d3b.sim]*0.8, d3b.p[1, d3b.sim]*0.05)
    d3b.p[3, d3b.sim] <- runif(nsim, d3b.p[2, d3b.sim]*0.8, d3b.p[2, d3b.sim]*0.05)
    d3b.p[4, d3b.sim] <- runif(nsim, 1, 6)
    d3b.p[5, d3b.sim] <- runif(nsim, d3b.p[4, d3b.sim] - 1, d3b.p[4, d3b.sim] - 0.2)
    d3b.p[6, d3b.sim] <- runif(nsim, d3b.p[5, d3b.sim] - 1, d3b.p[5, d3b.sim] - 0.2)
    d3b.all <- apply(d3b.p, 2, function(p, x) pred.d3b(x, p), x = time.samp)
    d3b.tmax <- apply(d3b.all, 2, function(x) time.samp[which(x == max(x))])
    d3b.sim <- which(d3b.tmax > 12)
    if (length(d3b.sim) == 0) break
    else print(paste("d3b repeat", length(d3b.sim)))
  }
  d3b.t <- c(0, 0.5, 1, 2, 4, 8, 12, 16, 24)

  #plot.rdata(d3b, time.samp, niter, -1, log = F)

# One Compartment Kinetics w/ Absorption
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a.p <- matrix(nrow = 3, ncol = niter)
  d1a.sim <- 1:niter
  repeat {
    nsim <- length(d1a.sim)
    d1a.p[1, d1a.sim] <- runif(nsim, -1, -0.1)
    d1a.p[2, d1a.sim] <- runif(nsim, d1a.p[1, d1a.sim]*2, d1a.p[1, d1a.sim]*1.1)
    d1a.p[3, d1a.sim] <- runif(nsim, 1, 6)
    d1a.all <- apply(d1a.p, 2, function(p, x) pred.d1a(x, p), x = time.samp)
    d1a.tmax <- apply(d1a.all, 2, function(x) time.samp[which(x == max(x))])
    d1a.sim <- which(d1a.tmax > 12)
    if (length(d1a.sim) == 0) break
    else print(paste("d1a repeat", length(d1a.sim)))
  }
  d1a.t <- c(0, 0.5, 1, 3, 5, 7, 10, 16, 24)
  #plot.rdata(d1a, time.samp, niter, -1, log = F)

# Two Compartment Kinetics w/ Absorption
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  d2a.p <- matrix(nrow = 5, ncol = niter)
  d2a.sim <- 1:niter
  repeat {
    nsim <- length(d2a.sim)
    d2a.p[1, d2a.sim] <- runif(nsim, -1, -0.1)
    d2a.p[2, d2a.sim] <- runif(nsim, d2a.p[1, d2a.sim]*0.8, d2a.p[1, d2a.sim]*0.05)
    d2a.p[3, d2a.sim] <- runif(nsim, d2a.p[1, d2a.sim]*2, d2a.p[1, d2a.sim]*1.1)
    d2a.p[4, d2a.sim] <- runif(nsim, 1, 6)
    d2a.p[5, d2a.sim] <- runif(nsim, d2a.p[4, d2a.sim] - 1, d2a.p[4, d2a.sim] - 0.2)
    d2a.all <- apply(d2a.p, 2, function(p, x) pred.d2a(x, p), x = time.samp)
    d2a.tmax <- apply(d2a.all, 2, function(x) time.samp[which(x == max(x))])
    d2a.sim <- which(d2a.tmax > 12)
    if (length(d2a.sim) == 0) break
    else print(paste("d2a repeat", length(d2a.sim)))
  }
  d2a.t <- c(0, 0.5, 1, 3, 5, 8, 12, 16, 24)
  #plot.rdata(d2a, time.samp, niter, -1, log = F)

# Three Compartment Kinetics w/ Absorption
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }
  # pred.d3a.dv <- function(x, p) {
  #   p[1]^3*exp(p[1]*x + p[5]) + p[2]^3*exp(p[2]*x + p[6]) + p[3]^3*exp(p[3]*x + p[7]) - p[4]^3*exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  # }
  d3a.p <- matrix(nrow = 7, ncol = niter)
  d3a.sim <- 1:niter
  repeat {
    nsim <- length(d3a.sim)
    d3a.p[1, d3a.sim] <- runif(nsim, -1, -0.1)
    d3a.p[2, d3a.sim] <- runif(nsim, d3a.p[1, d3a.sim]*0.8, d3a.p[1, d3a.sim]*0.05)
    d3a.p[3, d3a.sim] <- runif(nsim, d3a.p[2, d3a.sim]*0.8, d3a.p[2, d3a.sim]*0.05)
    d3a.p[4, d3a.sim] <- runif(nsim, d3a.p[1, d3a.sim]*2, d3a.p[1, d3a.sim]*1.1)
    d3a.p[5, d3a.sim] <- runif(nsim, 1, 6)
    d3a.p[6, d3a.sim] <- runif(nsim, d3a.p[5, d3a.sim] - 1, d3a.p[5, d3a.sim] - 0.2)
    d3a.p[7, d3a.sim] <- runif(nsim, d3a.p[6, d3a.sim] - 1, d3a.p[6, d3a.sim] - 0.2)
    d3a.all <- apply(d3a.p, 2, function(p, x) pred.d3a(x, p), x = time.samp)
    d3a.tmax <- apply(d3a.all, 2, function(x) time.samp[which(x == max(x))])
    d3a.sim <- which(d3a.tmax > 12)
    if (length(d3a.sim) == 0) break
    else print(paste("d3a repeat", length(d3a.sim)))
  }
  #plot.rdata(d3a.all, time.samp, niter, -1, log = F)
  d3a.t <- c(0, 0.5, 1.5, 3, 5, 8, 12, 16, 24)

### Error models
# Two Compartment Kinetics
  if (sdev == 1) err.sig <- 0.05
  if (sdev == 2) err.sig <- 0.15
  if (sdev == 3) err.sig <- 0.05
  if (sdev == 4) err.sig <- 0.1
  if (sdev == 5) err.sig <- 0.15
  d2b.sub <- d2b.all[which(time.samp %in% d2b.t),]
  eps1 <- matrix(nrow = length(d2b.t), ncol = niter)
  eps2 <- matrix(nrow = length(d2b.t), ncol = niter)
  d2b.err <- 1:niter
  repeat {
    nsim <- length(d2b.err)
    eps1[, d2b.err] <- matrix(
      rnorm(n = length(d2b.t)*nsim, mean = 0, sd = err.sig),
      nrow = length(d2b.t), ncol = nsim
    )
    if (sdev <= 2) {
      d2b <- d2b.sub*(1 + eps1)
    } else {
      if (class(d2b.sub[, d2b.err]) == "matrix") {
        eps2[, d2b.err] <- apply(d2b.sub[, d2b.err], 2, function(x) {
          rnorm(n = length(d2b.t), mean = 0, sd = tail(x, 1)*err.sig)
        })
      } else {
        eps2[, d2b.err] <- rnorm(
          n = length(d2b.t), mean = 0, sd = tail(d2b.sub[, d2b.err], 1)*err.sig
        )
      }
      d2b <- d2b.sub*(1 + eps1) + eps2
    }
    d2b.err.tmax <- apply(d2b, 2, function(x) d2b.t[which(x == max(x))])
    d2b.err <- which(d2b.err.tmax > 12)
    if (length(d2b.err) == 0) break
    else print(paste("d2b err repeat", length(d2b.err)))
  }

# Three Compartment Kinetics
  d3b.sub <- d3b.all[which(time.samp %in% d3b.t),]
  eps1 <- matrix(nrow = length(d3b.t), ncol = niter)
  eps2 <- matrix(nrow = length(d3b.t), ncol = niter)
  d3b.err <- 1:niter
  repeat {
    nsim <- length(d3b.err)
    eps1[, d3b.err] <- matrix(
      rnorm(n = length(d3b.t)*nsim, mean = 0, sd = err.sig),
      nrow = length(d3b.t), ncol = nsim
    )
    if (sdev <= 2) {
      d3b <- d3b.sub*(1 + eps1)
    } else {
      if (class(d3b.sub[, d3b.err]) == "matrix") {
        eps2[, d3b.err] <- apply(d3b.sub[, d3b.err], 2, function(x) {
          rnorm(n = length(d3b.t), mean = 0, sd = tail(x, 1)*err.sig)
        })
      } else {
        eps2[, d3b.err] <- rnorm(
          n = length(d3b.t), mean = 0, sd = tail(d3b.sub[, d3b.err], 1)*err.sig
        )
      }
      d3b <- d3b.sub*(1 + eps1) + eps2
    }
    d3b.err.tmax <- apply(d3b, 2, function(x) d3b.t[which(x == max(x))])
    d3b.err <- which(d3b.err.tmax > 12)
    if (length(d3b.err) == 0) break
    else print(paste("d3b err repeat", length(d3b.err)))
  }

# One Compartment Kinetics w/ Absorption
  d1a.sub <- d1a.all[which(time.samp %in% d1a.t),]
  eps1 <- matrix(nrow = length(d1a.t), ncol = niter)
  eps2 <- matrix(nrow = length(d1a.t), ncol = niter)
  d1a.err <- 1:niter
  repeat {
    nsim <- length(d1a.err)
    eps1[, d1a.err] <- matrix(
      rnorm(n = length(d1a.t)*nsim, mean = 0, sd = err.sig),
      nrow = length(d1a.t), ncol = nsim
    )
    if (sdev <= 2) {
      d1a <- d1a.sub*(1 + eps1)
    } else {
      if (class(d1a.sub[, d1a.err]) == "matrix") {
        eps2[, d1a.err] <- apply(d1a.sub[, d1a.err], 2, function(x) {
          rnorm(n = length(d1a.t), mean = 0, sd = tail(x, 1)*err.sig)
        })
      } else {
        eps2[, d1a.err] <- rnorm(
          n = length(d1a.t), mean = 0, sd = tail(d1a.sub[, d1a.err], 1)*err.sig
        )
      }
      d1a <- d1a.sub*(1 + eps1) + eps2
    }
    d1a.err.tmax <- apply(d1a, 2, function(x) d1a.t[which(x == max(x))])
    d1a.err <- which(d1a.err.tmax > 12)
    if (length(d1a.err) == 0) break
    else print(paste("d1a err repeat", length(d1a.err)))
  }

# Two Compartment Kinetics w/ Absorption
  d2a.sub <- d2a.all[which(time.samp %in% d2a.t),]
  eps1 <- matrix(nrow = length(d2a.t), ncol = niter)
  eps2 <- matrix(nrow = length(d2a.t), ncol = niter)
  d2a.err <- 1:niter
  repeat {
    nsim <- length(d2a.err)
    eps1[, d2a.err] <- matrix(
      rnorm(n = length(d2a.t)*nsim, mean = 0, sd = err.sig),
      nrow = length(d2a.t), ncol = nsim
    )
    if (sdev <= 2) {
      d2a <- d2a.sub*(1 + eps1)
    } else {
      if (class(d2a.sub[, d2a.err]) == "matrix") {
        eps2[, d2a.err] <- apply(d2a.sub[, d2a.err], 2, function(x) {
          rnorm(n = length(d2a.t), mean = 0, sd = tail(x, 1)*err.sig)
        })
      } else {
        eps2[, d2a.err] <- rnorm(
          n = length(d2a.t), mean = 0, sd = tail(d2a.sub[, d2a.err], 1)*err.sig
        )
      }
      d2a <- d2a.sub*(1 + eps1) + eps2
    }
    d2a.err.tmax <- apply(d2a, 2, function(x) d2a.t[which(x == max(x))])
    d2a.err <- which(d2a.err.tmax > 12)
    if (length(d2a.err) == 0) break
    else print(paste("d2a err repeat", length(d2a.err)))
  }

# Three Compartment Kinetics w/ Absorption
  d3a.sub <- d3a.all[which(time.samp %in% d3a.t),]
  eps1 <- matrix(nrow = length(d3a.t), ncol = niter)
  eps2 <- matrix(nrow = length(d3a.t), ncol = niter)
  d3a.err <- 1:niter
  repeat {
    nsim <- length(d3a.err)
    eps1[, d3a.err] <- matrix(
      rnorm(n = length(d3a.t)*nsim, mean = 0, sd = err.sig),
      nrow = length(d3a.t), ncol = nsim
    )
    if (sdev <= 2) {
      d3a <- d3a.sub*(1 + eps1)
    } else {
      if (class(d3a.sub[, d3a.err]) == "matrix") {
        eps2[, d3a.err] <- apply(d3a.sub[, d3a.err], 2, function(x) {
          rnorm(n = length(d3a.t), mean = 0, sd = tail(x, 1)*err.sig)
        })
      } else {
        eps2[, d3a.err] <- rnorm(
          n = length(d3a.t), mean = 0, sd = tail(d3a.sub[, d3a.err], 1)*err.sig
        )
      }
      d3a <- d3a.sub*(1 + eps1) + eps2
    }
    d3a.err.tmax <- apply(d3a, 2, function(x) d3a.t[which(x == max(x))])
    d3a.err <- which(d3a.err.tmax > 12)
    if (length(d3a.err) == 0) break
    else print(paste("d3a err repeat", length(d3a.err)))
  }
