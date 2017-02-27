# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)
  library(splines)

# Set knots
  K <- c(5, 10)

# Set up one-compartment kinetic model
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  ERR <- 0.3	# Standard deviation of error
  dose <- 50	# mg
  times <- seq(from = 0,to = 24,by = 0.25)	# Time sequence for simulating concentrations
  sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)	# Sampling times
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))*(1+rnorm(n = length(times),mean = 0,sd = ERR))
  sample.conc <- conc[times %in% sample.times]
  onecomp <- data.frame(
    time = sample.times,
    conc = sample.conc
  )

# Breakdown of bs() function
# Arguments
  x <- onecomp$time
  df <- NULL
  knots <- c(K)
  degree <- 1
  intercept <- FALSE
  Boundary.knots <- range(x)

  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1)
    stop("'degree' must be integer >= 1")
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))

  splineDesign(Aknots, x, ord)
