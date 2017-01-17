# Fit a 1-compartment model to some data using a least squares objective function
# ------------------------------------------------------------------------------
# Remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	set.seed(123)

# ------------------------------------------------------------------------------
# Create some "data" that will be fitted
# The parameters defined here can be our "answers" for the fitting process
	CL <- 10	# Clearance, 10 L/h
	V <- 50	# Volume of distribution, 50 L
	KA <- 0.5	# Absorption rate constant, h^-1
	ERR <- 0.3	# Standard deviation of error
	dose <- 50	# mg
	times <- seq(from = 0,to = 24,by = 0.25)	# Time sequence for simulating concentrations
	sample.times <- c(0,0.25,0.5,1,2,4,8,12,24)	# Sampling times
	conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))*(1+rnorm(n = length(times),mean = 0,sd = ERR))
	sample.conc <- conc[times %in% sample.times]
	plot(sample.conc ~ sample.times, col = "red")

# ------------------------------------------------------------------------------
# Create a function for fitting a 1-compartment model to our "observed data"
	one.comp.function <- function(par, Cobs, dose, time) {
		CLfit <- par[1]
		Vfit <- par[2]
		KAfit <- par[3]
		Chat <- dose*KAfit/(Vfit*(KAfit-(CLfit/Vfit)))*(exp(-CLfit/Vfit*time)-exp(-KAfit*time))
		err <- Cobs-Chat
		squ.err <- err^2
		sse <- sum(squ.err)	# Sum of squared errors to be minimised
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
		Cobs = sample.conc, dose = dose, time = sample.times  # Function arguments
	)
	result	# Print the result

# ------------------------------------------------------------------------------
# Plot the model fit on top of the observed data
	CLfit <- result$par[1]	# Extract the fitted parameter estimates
	Vfit <- result$par[2]
	KAfit <- result$par[3]
	conc.sim <- dose*KAfit/(Vfit*(KAfit-(CLfit/Vfit)))*(exp(-CLfit/Vfit*times)-exp(-KAfit*times))
	points(conc.sim ~ times,type = "l",col = "blue")
