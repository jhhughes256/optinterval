# Fixing predsumexp as per Issue #8
# -----------------------------------------------------------------------------
# First set up the functions
# pred.sumexp with for loop, no changes
  fn.for.old <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    for (i in 1:n) {
      if (i == 1) y <- x[i]^d*exp(x[i]*t + x[n+i])
      else if (i != n | a == 0) y <- y + x[i]^d*exp(x[i]*t + x[n+i])
      else if (a == 1) y <- y - x[i]^d*exp(x[i]*t)*sum(exp(x[(n+1):(2*n-1)]))
    }
    return(y)
  }

# pred.sumexp vectorised, no changes
  fn.vec.old <- function(p, x, d = 0) {
    l <- length(p)
    a <- ifelse(l %% 2 == 0, F, T)
    n <- ceiling(l/2)
    m <- p[1:(n-a)]
    b <- p[(n+1):l]
    y <- colSums(exp(b)*outer(m, x, function(m, x) m^d*exp(m*x))) -
      a*p[n]^d*exp(p[n]*x)*sum(exp(b))
    return(y)
  }

# pred.sumexp with for loop, changed for issue #8
  fn.for.new <- function(x, t, d = 0) {
    l <- length(x)
    a <- ifelse(l %% 2 == 0, 0, 1)
    n <- ceiling(l/2)
    m <- x[1:n]  # new
    ord <- order(m, decreasing = T)  # new
    p <- c(m[ord], x[(n+1):l])  # new
    for (i in 1:n) {
      if (i == 1) y <- p[i]^d*exp(p[i]*t + p[n+i])
      else if (i != n | a == 0) y <- y + p[i]^d*exp(p[i]*t + p[n+i])
      else if (a == 1) y <- y - p[i]^d*exp(p[i]*t)*sum(exp(p[(n+1):(2*n-1)]))
    }
    return(y)
  }

# pred.sumexp vectorised, changed for issue #8
  fn.vec.new <- function(p, t, d = 0) {
    l <- length(p)
    a <- ifelse(l %% 2 == 0, F, T)
    n <- ceiling(l/2)
    m <- p[1:n]
    ord <- order(m, decreasing = T)
    ke <- m[ord][1:(n-a)]
    ka <- m[ord][n]
    int <- p[(n+1):l]
    y <- colSums(exp(int)*outer(ke, t, function(ke, t) ke^d*exp(ke*t))) -
      a*ka^d*exp(ka*t)*sum(exp(int))
    return(y)
  }

# -----------------------------------------------------------------------------
# Test to ensure that the functions are working as intended
  time.samp <- seq(0, 24, by = 0.1)
# Test environment 1.1 - Regular One Compartment Absorption:
  d1a1.p <- c(-0.1, -0.4, 4)
  d1a1 <- data.frame(
    fo = fn.for.old(d1a1.p, time.samp),
    vo = fn.vec.old(d1a1.p, time.samp),
    fn = fn.for.new(d1a1.p, time.samp),
    vn = fn.vec.new(d1a1.p, time.samp)
  )
  with(d1a1, c(all.equal(fo, vo), all.equal(fo, fn), all.equal(fo, vn)))

# Test environment 1.2 - Regular Two Compartment Absorption:
  d2a1.p <- c(-0.01, -0.2, -0.4, log(exp(4)*0.2), log(exp(4)*0.8))
  d2a1 <- data.frame(
    fo = fn.for.old(d2a1.p, time.samp),
    vo = fn.vec.old(d2a1.p, time.samp),
    fn = fn.for.new(d2a1.p, time.samp),
    vn = fn.vec.new(d2a1.p, time.samp)
  )
  with(d2a1, c(all.equal(fo, vo), all.equal(fo, fn), all.equal(fo, vn)))

# Test environment 1.3 - Regular Three Compartment Absorption:
  d3a1.p <- c(-0.001, -0.08, -0.5, -0.8, log(exp(4)*0.15), log(exp(4)*0.25), log(exp(4)*0.6))
  d3a1 <- data.frame(
    fo = fn.for.old(d3a1.p, time.samp),
    vo = fn.vec.old(d3a1.p, time.samp),
    fn = fn.for.new(d3a1.p, time.samp),
    vn = fn.vec.new(d3a1.p, time.samp)
  )
  with(d3a1, c(all.equal(fo, vo), all.equal(fo, fn), all.equal(fo, vn)))

# Test environment 1.4 - Regular Two Compartment Bolus:
  d2b1.p <- c(-0.05, -0.5, 5, 6)
  d2b1 <- data.frame(
    fo = fn.for.old(d2b1.p, time.samp),
    vo = fn.vec.old(d2b1.p, time.samp),
    fn = fn.for.new(d2b1.p, time.samp),
    vn = fn.vec.new(d2b1.p, time.samp)
  )
  with(d2b1, c(all.equal(fo, vo), all.equal(fo, fn), all.equal(fo, vn)))

# Test environment 1.5 - Regular Three Compartment Bolus:
  d3b1.p <- c(-0.01, -0.1, -0.4, 4.1, 4.7, 6)
  d3b1 <- data.frame(
    fo = fn.for.old(d3b1.p, time.samp),
    vo = fn.vec.old(d3b1.p, time.samp),
    fn = fn.for.new(d3b1.p, time.samp),
    vn = fn.vec.new(d3b1.p, time.samp)
  )
  with(d3b1, c(all.equal(fo, vo), all.equal(fo, fn), all.equal(fo, vn)))

# -----------------------------------------------------------------------------
# Test to ensure that the functions solve the problems outlined in issue #8
# Test environment 2.1 - Regular One Compartment Absorption:
  d1a2.p <- c(-0.2, -0.1, 4)
  d1a2 <- data.frame(
    fo = fn.for.old(d1a2.p, time.samp),
    fn = fn.for.new(d1a2.p, time.samp),
    vn = fn.vec.new(d1a2.p, time.samp)
  )
  with(d1a2, c(all.equal(fo, fn), all.equal(fn, vn)))

# Test environment 2.2 - Regular Two Compartment Absorption:
  d2a2.p <- c(-0.8, -0.01, -0.4, log(exp(4)*0.8), log(exp(4)*0.2))
  d2a2 <- data.frame(
    fo = fn.for.old(d2a2.p, time.samp),
    fn = fn.for.new(d2a2.p, time.samp),
    vn = fn.vec.new(d2a2.p, time.samp)
  )
  with(d2a2, c(all.equal(fo, fn), all.equal(fn, vn)))

# Test environment 2.3 - Regular Three Compartment Absorption:
  d3a2.p <- c(-0.001, -0.8, -0.08, -0.5, log(exp(4)*0.15), log(exp(4)*0.25), log(exp(4)*0.6))
  d3a2 <- data.frame(
    fo = fn.for.old(d3a2.p, time.samp),
    fn = fn.for.new(d3a2.p, time.samp),
    vn = fn.vec.new(d3a2.p, time.samp)
  )
  with(d3a2, c(all.equal(fo, fn), all.equal(fn, vn)))

# -----------------------------------------------------------------------------
# Benchmarking functions
  benchmark(fn.for.new(d1a1.p, time.samp), replications = 10000)
  benchmark(fn.vec.new(d1a1.p, time.samp), replications = 10000)

  benchmark(fn.for.new(d3a1.p, time.samp), replications = 10000)
  benchmark(fn.vec.new(d3a1.p, time.samp), replications = 10000)

  benchmark(fn.for.new(d2b1.p, time.samp), replications = 10000)
  benchmark(fn.vec.new(d2b1.p, time.samp), replications = 10000)

# It seems vectorising the function doesn't speed it up, worth trying though! :)
