# Script for data creation of 6 drugs for testing
# -----------------------------------------------------------------------------
# Set up functions and global objects
  times <- seq(from = 0, to = 24, by = 0.25)

  onecomp.func <- function(par, x) {
    A <- par$KA/(par$V*(par$KA-(par$CL/par$V)))
    amt <- par$dose*par$F1
    data.frame(
      x = x,
      y = amt*A*(exp(-par$CL/par$V*x)-exp(-par$KA*x))
    )
  }

  twocomp.func <- function(par, x) {
    k10 <- par$CL/par$V1
    k12 <- par$Q/par$V1
    k21 <- par$Q/par$V2
    apb <- k10 + k12 + k21            # alpha + beta
    amb <- k10*k21                # alpha * beta
    alpha <- (apb + sqrt(apb^2 - 4*amb))/2
    beta <- (apb - sqrt(apb^2 - 4*amb))/2
    A <- par$KA*(k21 - alpha)/(par$V1*(par$KA - alpha)*(beta - alpha))
    B <- par$KA*(k21 - beta)/(par$V1*(par$KA - beta)*(alpha - beta))
    amt <- par$dose*par$F1
    data.frame(
      x = x,
      y = amt*(A*exp(-alpha*x) + B*exp(-beta*x) - (A+B)*exp(-par$KA*x))
    )
  }

  threecomp.func <- function(par, x) {
    k10 <- par$CL/par$V1
    k12 <- par$Q2/par$V1
    k21 <- par$Q2/par$V2
    k13 <- par$Q3/par$V1
    k31 <- par$Q3/par$V3

    a0 <- k10*k21*k31
    a1 <- k10*k31 + k21*k31 + k21*k13 + k10*k21 + k31*k12
    a2 <- k10 + k12 + k13 + k21 + k31

    p <- a1 - a2^2/3
    q <- 2*a2^3/27 - a1*a2/3 + a0
    r1 <- sqrt(-p^3/27)
    r2 <- 2*r1^(1/3)

    phi <- acos(-q/2*r1)/3
    alpha <- -(cos(phi)*r2 - a2/3)
    beta <- -(cos(phi + 2*pi/3)*r2 - a2/3)
    gamma <- -(cos(phi + 4*pi/3)*r2 - a2/3)

    Anum <- par$KA*(k21 - alpha)*(k31 - alpha)
    Aden <- par$V1*(par$KA - alpha)*(alpha - beta)*(alpha - gamma)
    A <- Anum/Aden
    Bnum <- par$KA*(k21 - beta)*(k31 - beta)
    Bden <- par$V1*(par$KA - beta)*(beta - alpha)*(beta - gamma)
    B <- Bnum/Bden
    Cnum <- par$KA*(k21 - gamma)*(k31 - gamma)
    Cden <- par$V1*(par$KA - gamma)*(gamma - beta)*(gamma - alpha)
    C <- Cnum/Cden
    amt <- par$dose*par$F1
    data.frame(
      x = x,
      y = amt*(A*exp(-alpha*x) + B*exp(-beta*x) + C*exp(-gamma*x) - (A+B+C)*exp(-par$KA*x))
    )
  }

# -----------------------------------------------------------------------------
# Drug 1: One Compartment No Flip-Flop
  drug1 <- list(
    CL = 10,
    V = 50,
    KA = 0.5,
    F1 = 1,
    dose = 50
  )
  pkdata1 <- onecomp.func(drug1, times)
# -----------------------------------------------------------------------------
# Drug 2: One Compartment Flip-Flop
  drug2 <- list(
    CL = 8,
    V = 10,
    KA = 0.3,
    F1 = 1,
    dose = 50
  )
  pkdata2 <- onecomp.func(drug2, times)
# -----------------------------------------------------------------------------
# Drug 3: Two Compartment
  drug3 <- list(
    CL = 10,
    V1 = 30,
    Q = 15,
    V2 = 100,
    KA = 0.5,
    F1 = 1,
    dose = 50
  )
  pkdata3 <- twocomp.func(drug3, times)
# -----------------------------------------------------------------------------
# Drug 4: Two Compartment
  drug4 <- list(
    CL = 5,
    V1 = 30,
    Q = 20,
    V2 = 200,
    KA = 0.3,
    F1 = 1,
    dose = 50
  )
  pkdata4 <- twocomp.func(drug4, times)
# -----------------------------------------------------------------------------
# Drug 5: Three Compartment
  drug5 <- list(
    CL = 10,
    V1 = 30,
    Q2 = 20,
    V2 = 50,
    Q3 = 10,
    V3 = 200,
    KA = 0.5,
    F1 = 1,
    dose = 50
  )
  pkdata5 <- threecomp.func(drug5, times)
# -----------------------------------------------------------------------------
# Drug 6: Three Compartment
  drug6 <- list(
    CL = ,
    V1 = ,
    Q2 = ,
    V2 = ,
    Q3 = ,
    V3 = ,
    KA = ,
    F1 = 1,
    dose = 50
  )
  pkdata6 <- threecomp.func(drug6, times)
