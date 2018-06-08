# Remove all current objects in the workspace
  rm(list = ls(all = TRUE))
  set.seed(123)

# Then again why use splines?
# PK Model
  CL <- 10	# Clearance, 10 L/h
  V <- 50	# Volume of concribution, 50 L
  KA <- 0.5	# Absorption rate constant, h^-1
  dose <- 50	# mg
  times <- seq(from = 0,to = 24,by = 0.5)	# Time sequence for simulating concentrations
  conc <- dose*KA/(V*(KA-(CL/V)))*(exp(-CL/V*times)-exp(-KA*times))
  plot(conc ~ times)
  
# Calculate slopes using basic difference in concentration / difference in times
  out2$m2 <- diff(conc)/diff(times)

# Calculate intercepts by tracing the slope back to time = 0
  out2$c2 <- head(conc, -1) - head(times, -1)*out$m2
  
# Turn this into a function
  interpLinear <- function(x, y) {
    m <- diff(conc)/diff(times)
    data.frame(
      K = paste(head(x, -1), tail(x, -1), sep = "-" ),
      m = m,
      c = head(conc, -1) - head(times, -1)*m
    )
  }
  
# Plot and call it a day
  out3 <- interpLinear(times, conc)
  plot(conc ~ times)
  for(i in 1:4) {
    lines(out3$m[i]*times + out3$c[i] ~ times)
  }
  