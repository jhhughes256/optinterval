# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

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

# Function
# Determine order and check to make sure degree isn't less than 1
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1)
      stop("'degree' must be integer >= 1")
# Setting up environment
  nx <- names(x)
  x <- as.vector(x)
# Check if x to see if there are any na's and remove them
  nax <- is.na(x)
  if (nas <- any(nax))
      x <- x[!nax]
# Check to see if x is able to be converted to boundary knots
# Then sort the knots and check what points are outside the knots
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  }
  else FALSE
# Check if degrees of freedom were specified (instead of knots)
# Knot interested in this part
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept)
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
        ord - (1L - intercept)), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots +
        2L)[-c(1L, nIknots + 2L)]
      quantile(x[!outside], knots)
    }
  }
# Write all knots to an object and sort in ascending order
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
# If there are any points outside the boundary.knots
  if (any(outside)) {
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree -
      1L))
    e <- 1/4
    if (any(ol)) {
      k.pivot <- (1 - e) * Boundary.knots[1L] + e * Aknots[ord +
        1]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree,
        "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
        derivs)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      k.pivot <- (1 - e) * Boundary.knots[2L] + e * Aknots[length(Aknots) -
        ord]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree,
        "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord,
        derivs)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splineDesign(Aknots, x[inside],
        ord)
  }
  else basis <- splineDesign(Aknots, x, ord)  #THIS IS THE FUNCTION WE NEED TO HONE IN ON
  if (!intercept) 
      basis <- basis[, -1L, drop = FALSE]
  n.col <- ncol(basis)
  if (nas) {
      nmat <- matrix(NA, length(nax), n.col)
      nmat[!nax, ] <- basis
      basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots,
      Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
