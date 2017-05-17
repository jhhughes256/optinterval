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
  time.samp <- seq(0, 48, by = 2)
# One Compartment Kinetics
  onedata <- data.frame(
    time = time.samp,
    line1 = -0.5*time.samp + 5
  )
  onedata$sumexp <- exp(onedata$line1)
# Two Compartment Kinetics
  twodata <- data.frame(
    time = time.samp,
    line1 = -0.5*time.samp + 6,
    line2 = -0.05*time.samp + 5
  )
  twodata$sumexp <- with(twodata, exp(line1) + exp(line2))
# Three Compartment Kinetics
# One Compartment Kinetics w/ Absorption
  onedata.abs <- data.frame(
    time = time.samp,
    line1 = -0.2*time.samp + 4,
    line2 = -0.1*time.samp + 4
  )
  onedata.abs$sumexp <- with(onedata.abs, exp(line2) - exp(line1))
# Two Compartment Kinetics w/ Absorption
  twodata.abs <- data.frame(
    time = time.samp,
    line1 = -0.4*time.samp + 4,
    line2 = -0.2*time.samp + log(exp(4)*0.8),
    line3 = -0.02*time.samp + log(exp(4)*0.2)
  )
  twodata.abs$sumexp <- with(twodata.abs, exp(line2) + exp(line3) - exp(line1))
# Three Compartment Kinetics w/ Absorption
