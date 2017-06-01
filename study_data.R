# Sum of exponential data for sourcing
# -----------------------------------------------------------------------------
# The datasets in order are:
#   onedata - one compartment kinetics
#   twodata - two compartment kinetics
#   thrdata - three compartment kinetics
#   onedata.abs - one compartment kinetics w/ absorption
#   twodata.abs - two compartment kinetics w/ absorption
#   thrdata.abs - three compartment kinetics w/ absorption
# -----------------------------------------------------------------------------
  time.samp <- seq(0, 48, by = 0.1)
# One Compartment Kinetics
  pred.d1b <- function(x, p) {
    exp(p[1]*x + p[2])
  }
  d1b.p <- c(-0.5, 5)
  d1b <- data.frame(
    time = time.samp,
    conc = pred.d1b(time.samp, d1b.p)
  )
  #with(d1b, plot(time, log(conc)))
# Two Compartment Kinetics
  pred.d2b <- function(x, p) {
    exp(p[1]*x + p[3]) + exp(p[2]*x + p[4])
  }
  d2b.p <- c(-0.5, -0.05, 6, 5)
  d2b <- data.frame(
    time = time.samp,
    conc = pred.d2b(time.samp, d2b.p)
  )
  #with(d2b, plot(time, log(conc)))
# Three Compartment Kinetics
  pred.d3b <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) + exp(p[3]*x + p[6])
  }
  d3b.p <- c(-0.5, -0.05, -0.15, 6, 3.5, 4.7)
  thrdata <- data.frame(
    time = time.samp,
    conc = pred.d3b(time.samp, d3b.p)
  )
  #with(d3b, plot(time, log(conc)))
# One Compartment Kinetics w/ Absorption
  pred.d1a <- function(x, p) {
    exp(p[1]*x + p[3]) - exp(p[2]*x + p[3])
  }
  d1a <- data.frame(
    time = time.samp,
    conc = pred.d1a(time.samp, d1a.p)
  )
  onedata.abs$conc <- with(onedata.abs, exp(line2) - exp(line1))
  #with(onedata.abs, plot(time, log(conc)))
# Two Compartment Kinetics w/ Absorption
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + p[6])
  }
  d2a <- data.frame(
    time = time.samp,
    conc = pred.d2a(time.samp, d2a.p)
  )
  twodata.abs$conc <- with(twodata.abs, exp(line2) + exp(line3) - exp(line1))
  #with(twodata.abs, plot(time, log(conc)))
# Three Compartment Kinetics w/ Absorption
  pred.d3a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + p[6])
  }
  d3a <- data.frame(
    time = time.samp,
    conc = pred.d3a(time.samp, d3a.p)
  )
  thrdata.abs$conc <- with(thrdata.abs, exp(line2) + exp(line3) + exp(line4) - exp(line1))
  #with(thrdata.abs, plot(time[-1], log(conc)[-1]))
