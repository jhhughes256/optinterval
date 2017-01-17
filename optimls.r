# Fit a 1-compartment model to some data using a least squares objective function
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	set.seed(123)

	simulate.1comp.abs <- function(x) {
		k10 <- x$CL/x$V2
		A <- x$KA/(x$V2*(x$KA - k10))
		x$AMT*A*(exp(-k10*x$TIME) - exp(-x$KA*x$TIME))
	}

# ------------------------------------------------------------------------------
# Create some "data" that will be fitted
# The parameters defined here can be our "answers" for the fitting process
	obs.param <- list(
		AMT = 50,  # mg
	  CL = 10,  # Clearance, 10 L/h
		V2 = 50,  # Volume of distribution, 50 L
		KA = 0.5,  # Absorption rate constant, h^-1
		TIME = seq(from = 0,to = 24,by = 0.25)  # Time sequence for simulating concentrations
	)
	ERR <- 0.3
	conc <- simulate.1comp.abs(obs.param)*(1 + rnorm(n = length(obs.param$TIME),mean = 0,sd = ERR))
#	sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)  # Sampling times
	sample.times <- c(0:24)
	sample.conc <- conc[obs.param$TIME %in% sample.times]
	plot(sample.conc ~ sample.times, col = "red")
	testdf <- data.frame(time = sample.times, Cp = sample.conc)

	#spline <- smooth.spline(testdf$time, testdf$Cp, df = 5)
  spline <- splinefun(testdf$time, testdf$Cp, "natural")
	library(splines)
	spline <- interpSpline(testdf$time, testdf$Cp)

	ggplot(testdf, aes(x = time, y = Cp)) + stat_function(fun = spline) + geom_point() + scale_x_continuous(limits = c(0,24))

	summary(fm1 <- lm(Cp ~ ns(time, df = 4), data = testdf))
	predict(fm1, time = sample.times)

# ------------------------------------------------------------------------------
# Create a function for fitting a 1-compartment model to our "observed data"
	one.comp.function <- function(par, Cobs, dose, time) {
		param <- list(
			AMT = dose,
			CL = par[1],
			V2 = par[2],
			KA = par[3],
			TIME = time
		)
		Chat <- simulate.1comp.abs(param)
		err <- Cobs-Chat
		squ.err <- err^2
		sse <- sum(squ.err)  # Sum of squared errors to be minimised
		sse
	}

# Initial estimates
	init.par <- c(15,30,0.2)
  #par <- init.par

# Minimise the sum of squared errors (error between fitted concs and observations)
# Using the optim function in R
	result <- optim(
		init.par,  # Initial parameter estimates
		one.comp.function,  # Fitting function
		Cobs = sample.conc, dose = obs.param$AMT, time = sample.times  # Function arguments
	)
	result  # Print the result

# ------------------------------------------------------------------------------
# Plot the model fit on top of the observed data
	sim.param <- list(
		AMT = obs.param$AMT,
		CL = result$par[1],
		V2 = result$par[2],
		KA = result$par[3],
		TIME = obs.param$TIME
	)
	conc.sim <- simulate.1comp.abs(sim.param)
	points(conc.sim ~ sim.param$TIME,type = "l",col = "blue")


	nobs <- length(sample.times)-1
	onecomp.fn <- function(x, AMT, CL, V2, KA) AMT*KA/(V2*(KA - CL/V2))*(exp(-CL/V2*x) - exp(-KA*x))

	ddfn <- function(x, .par, .fun) {
		.ddfun <- D(D(body(.fun), "x"), "x")
		with(.par, eval(.ddfun))
	}
	ffn <- function(x, .par, .fun) {
		with(.par, integrate(.fun, x[1], x[2], AMT = AMT, CL = CL, V2 = V2, KA = KA))
	}

	trap.par <- sim.param[1:4]
	trap.err <- ldply(seq_len(nobs), function(i) {
		t <- sample.times[c(i,i+1)]
		abs.err <- -(t[2] - t[1])^3*ddfn(t[1], trap.par, onecomp.fn)/12
		c(a = t[1],
			b = t[2],
			abs = abs(abs.err),
			rel = abs(abs.err/ffn(c(t[1],t[2]), trap.par, onecomp.fn)$value),
			auc = ffn(c(t[1],t[2]), trap.par, onecomp.fn)$value
		)
	})

	plot(trap.err$b, trap.err$abs)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Fit a 2-compartment model to some data using a least squares objective function
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	set.seed(123)

	simulate.2comp.abs <- function(x) {
		k10 <- x$CL/x$V2
		k12 <- x$Q/x$V2
		k21 <- x$Q/x$V3
		apb <- k10+k12+k21            # alpha + beta
		amb <- k10*k21                # alpha * beta
		alpha <- ((apb)+sqrt((apb)^2-4*amb))/2
		beta <- ((apb)-sqrt((apb)^2-4*amb))/2
		A <- x$KA*(k21-alpha)/(x$V2*(x$KA-alpha)*(beta-alpha))
		B <- x$KA*(k21-beta)/(x$V2*(x$KA-beta)*(alpha-beta))
		x$AMT*(A*exp(-alpha*x$TIME)+B*exp(-beta*x$TIME)-(A+B)*exp(-x$KA*x$TIME))
	}

# ------------------------------------------------------------------------------
# Create some "data" that will be fitted
# The parameters defined here can be our "answers" for the fitting process
	obs.param <- list(
		AMT = 50,  # mg
		CL = 10,  # Clearance, 10 L/h
		V2 = 50,  # Volume of distribution, 50 L
		Q = 8,  #
		V3 = 100,  #
		KA = 0.5,  # Absorption rate constant, h^-1
		TIME = seq(from = 0,to = 24,by = 0.25)  # Time sequence for simulating concentrations
	)
	ERR <- 0.3  # Standard deviation of error
	conc <- simulate.2comp.abs(obs.param)*(1 + rnorm(n = length(obs.param$TIME),mean = 0,sd = ERR))
	sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)  # Sampling times
	sample.conc <- conc[obs.param$TIME %in% sample.times]
	plot(sample.conc ~ sample.times, col = "red")

# ------------------------------------------------------------------------------
# Create a function for fitting a 2-compartment model to our "observed data"
	two.comp.function <- function(par, Cobs, dose, time) {
		param <- list(
			AMT = dose,
			CL = par[1],
			V2 = par[2],
			Q = par[3],
			V3 = par[4],
			KA = par[5],
			TIME = time
		)
		Chat <- simulate.2comp.abs(param)
		err <- Cobs-Chat
		squ.err <- err^2
		sse <- sum(squ.err)  # Sum of squared errors to be minimised
		sse
	}

# Initial estimates
	init.par <- c(15,30,4,80,0.2)
  #par <- init.par

# Minimise the sum of squared errors (error between fitted concs and observations)
# Using the optim function in R
	result <- optim(
		init.par,  # Initial parameter estimates
		two.comp.function,  # Fitting function
		Cobs = sample.conc, dose = obs.param$AMT, time = sample.times  # Function arguments
	)
	result  # Print the result

# ------------------------------------------------------------------------------
# Plot the model fit on top of the observed data
	sim.param <- list(
		AMT = obs.param$AMT,
		CL = result$par[1],
		V2 = result$par[2],
		Q = result$par[3],
		V3 = result$par[4],
		KA = result$par[5],
		TIME = obs.param$TIME
	)
	conc.sim <- simulate.2comp.abs(sim.param)
	points(conc.sim ~ sim.param$TIME,type = "l",col = "blue")
