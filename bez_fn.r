# Spline Primer
#As found from: https://cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf

# Function notation
	rm(list=ls(all=TRUE))
  source("F:/2016 PhD/Shiny/Gadgets/Function_Visualiser.r")
  plot(c(0,1),c(0,1),ann=F,bty="n",type="n",xaxt="n",yaxt="n")
  text(x = 0.15, y = 0.95, expression(B(x) == sum(b[i], i == 0, n) ~~ bgroup("(",atop(n,i),")") ~~ (1-x)^{n-i} * x^i))
  text(x = 0.17, y = 0.85, expression(B(x) == b[0]*(1-x) + b[1]*x ~~ ~~ ~~ linear))
  text(x = 0.27, y = 0.77, expression(B(x) == b[0]*(1-x)^2 + b[1]*(1-x)*x + b[2]*x^2 ~~ ~~ ~~ quadratic))
  text(x = 0.32, y = 0.69, expression(B(x) == b[0]*(1-x)^3 + b[1]*(1-x)^2*x + b[2]*(1-x)*x^2 + b[3]*x^3 ~~ ~~ ~~ cubic))

  # Bezier curve, n = 2  (quadratic)
  Bez.n1 <- function(x, b) b[1]*(1-x) + b[2]*x

# Bezier curve, n = 2  (quadratic)
  Bez.n2 <- function(x, b) b[1]*(1-x)^2 + 2*b[2]*(1-x)*x + b[3]*x^2
  funcVis(xscale = 1, func = "1*(1-x)^2 + -2*(1-x)*x + 2*x^2")

# Bezier curve, n = 3  (cubic)
  Bez.n3 <- function(x, b) b[1]*(1-x)^3 + 3*b[2]*(1-x)^2*x + 3*b[3]*(1-x)*x^2 + b[4]*x^3
  funcVis(xscale = 1, func = "1*(1-x)^3 + -6*(1-x)^2*x + 6*(1-x)*x^2 + -1*x^3")

# Bezier cruve, n = 4  (quartic)
  Bez.n4 <- function(x, b) b[1]*(1-x)^3 + 4*b[2]*(1-x)^3*x + 6*b[3]*(1-x)^2*x^2 + 4*b[4]*(1-x)*x^3 + b[5]*x^4
  funcVis(xscale = 1, func = "1*(1-x)^3 + -4*(1-x)^3*x + 12*(1-x)^2*x^2 + -18*(1-x)*x^3 + 4*x^4")

# B(x) = sum(b[i]*B(x)[i,n], i == 0, n)

n <- 5
b <- rnorm(n)
par(mar = rep(0,4))
