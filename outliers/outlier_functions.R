# A source script for the analysis of outliers
# Load libraries
  library(ggplot2)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
# Set global parameters
  time.samp <- seq(0, 24, 0.1)
  basic.samp <- seq(0, 24, length.out = 9)

# Exponential functions
  pred.d1b <- function(x, p) {
    exp(p[1]*x + p[2])
  }
  pred.d2b <- function(x, p) {
    exp(p[1]*x + p[3]) + exp(p[2]*x + p[4])
  }
  pred.d3b <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) + exp(p[3]*x + p[6])
  }
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[5]) + exp(p[2]*x + p[6]) + exp(p[3]*x + p[7]) - exp(p[4]*x + log(sum(exp(p[5]), exp(p[6]), exp(p[7]))))
  }

# Plotting functions
  plot.sumexp.facet <- function(data, t, n, interv, facet.scale = "fixed", log = F) {
    plotdata <- data.frame(
      id = rep(1:n, each = length(t)),
      time = rep(t, times = n),
      dv = as.vector(data)
    )
    xlim <- c(t[1], t[length(t)])
    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
    plotobj <- plotobj + geom_line(aes(x = time, y = dv), colour = "red")
    plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
    if (!log) plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
    else plotobj <- plotobj + scale_y_log10("Concentration (mg/mL)\n")
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
    plotobj <- plotobj + facet_wrap(~id, ncol = round(sqrt(n)), scales = facet.scale)
    return(plotobj)
  }

  plot.sumexp.overlay <- function(data, t, n, interv, log = F) {
    plotdata <- data.frame(
      id = rep(1:n, each = length(t)),
      time = rep(t, times = n),
      dv = as.vector(data)
    )
    xlim <- c(t[1], t[length(t)])
    plotobj <- NULL
    plotobj <- ggplot(data = plotdata)
    plotobj <- plotobj + ggtitle("Random Concentration Time Curves")
    plotobj <- plotobj + geom_point(aes(x = time, y = dv, colour = id))
    plotobj <- plotobj + geom_vline(xintercept = interv, colour = "green4", linetype = "dashed")
    if (!log) plotobj <- plotobj + scale_y_continuous("Concentration (mg/mL)\n")
    else plotobj <- plotobj + scale_y_log10("Concentration (mg/mL)\n")
    plotobj <- plotobj + scale_x_continuous("\nTime after dose (hrs)", lim = xlim)
    return(plotobj)
  }
