# Incorporate sigma into model fitting
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
            "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      gir.dir <- paste0(getwd(), "/GitRepos")
      reponame <- "optinterval"
    } else if (getwd() == wd[2]) {
      git.dir <- getwd()
      reponame <- "optinterval"
    } else if (getwd() == wd[3] | getwd() == wd[4]) {
      git.dir <- "E:/Hughes/Git"
      reponame <- "splines"
    }
    rm("wd")
  }

# Load packages
  library(GA)

# Source scripts to set up environment
  set.seed(256256)
  # niter <- 1000
  source(paste(git.dir, reponame, "study_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "study_rdata.R", sep = "/"))

# -----------------------------------------------------------------------------
  best.sumexp.lrt <- function(opt) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    i <- 1
    for (j in 2:length(opt$par)) {
      degf <- length(opt$par[[j]]) - length(opt$par[[i]])
      x <- ofv[i] - ofv[j]
      p <- pchisq(x, degf, lower.tail = F)
      if (p < 0.01) {
        i <- i + 1
      }
    }
    return(sapply(opt, function(x) x[i]))
  }

  best.sumexp.aic <- function(opt) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    k <- unlist(lapply(opt$par, length))
    aic <- ofv + 2*k
    return(sapply(opt, function(x) x[which(aic == min(aic))]))
  }

  best.sumexp.bic <- function(opt, nobs) {
    values <- unlist(opt$value)
    ofv <- values[which(names(values) == "ofv")]
    k <- unlist(lapply(opt$par, length))
    bic <- ofv + log(nobs)*k
    return(sapply(opt, function(x) x[which(bic == min(bic))]))
  }

  mle.sumexp.sig <- function(par, x, y, errmod, ga = F) {
  # Determine log likelihood of given model parameters
  # par = sum of exponential parameters
  # x = independent variable (time)
  # y = observed dependent variable (drug concentration)
  # errmod = error model to be used c("add", "prop", "both")
  # ga = genetic algorithm status
    nerr <- ifelse(errmod == "both", 2, 1) # set number of error parameters
    fit.par <- head(par, -nerr)  # define sum of exponentail parameters
    err.par <- tail(par, nerr)  # define error paramters
    z <- ifelse(ga, 2, -2)  # adjust ofv for minimisation or maximisation
    yhat <- pred.sumexp(fit.par, x)  # sum of exponential model prediction
  # Define standard deviation of normal distribution
    if (errmod == "add") sd <- err.par
    else if (errmod == "prop") sd <- abs(yhat)*err.par
    else if (errmod == "both") {
      add <- err.par[1]
      prop <- abs(yhat)*err.par[2]
      sd <- sqrt(add^2 + prop^2)
    }
    else stop("Please enter valid error model; \"add\", \"prop\" or \"both\"")
  # Determine log likelihood
    loglik <- suppressWarnings(dnorm(y, yhat, sd, log = T))
    return(z*sum(loglik))
  }

  optim.sumexp.new <- function(data, oral = F, nexp = 3) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
  # Set up objects
    res <- list(par = NULL, sumexp = NULL, value = NULL,
    error = NULL, hessian = NULL, message = NULL)
    # opt.par <- list(NULL)
    # opt.val <- list(NULL)
    # opt.gra <- list(NULL)
    # opt.con <- list(NULL)
    # opt.mes <- list(NULL)
    # opt.hes <- list(NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    repeat {
      lm.mod <- lm(log(y[lm.sub]) ~ x[lm.sub])
      lmres <- unname(lm.mod$coefficients)
      if (is.na(lmres[2])) {
        lmres <- c(max(y), -0.00001)
        break
      }
      if (lmres[2] < 0) break
      else lm.sub <- tail(lm.sub, -1)
    }
  # Estimate candidate model parameters
    for (i in 1:nexp) {
      gares <- try(ga("real-valued",
        mle.sumexp, x = x, y = y, ga = T, sigma = 0.01,
        min = c(rep(lmres[2]*50, i + oral), rep(lmres[1]-2, i)),
        max = c(rep(lmres[2]/50, i + oral), rep(lmres[1]+2, i)),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      ))
      if (class(gares) == "try-error") browser()
      optres <- try(optim(
        gares@solution[1, ],
        mle.sumexp,
        method = "BFGS", hessian = T,
        x = x, y = y, sigma = 0.01
      ))
      if (class(optres) == "try-error") {
        optres <- list(
          par = gares@solution[1,],
          value = mle.sumexp(gares@solution[1,], x, y, 0.01),
          counts = c("function" = 501, gradient = NA),
          convergence = 99,
          message = "zero gradient",
          hessian = matrix(NA,
            ncol = length(gares@solution[1,]),
            nrow = length(gares@solution[1,])
          )
        )
      }
    # Create output
      par.ord <- order.sumexp(optres$par, i, oral)
      res$par[[i]] <- optres$par
      res$sumexp[[i]] <- par.ord
      res$value[[i]] <- c(ofv = optres$value, optres$counts)
      res$error[[i]] <- c(0.01, "fixed")
      res$hessian[[i]] <- optres$hessian
      res$message[[i]] <- c(convergence = optres$convergence,
        message = ifelse(is.null(optres$message), "NULL", optres$message)
      )
    }
    res
  }

  #   opt.par[[i]] <- par.ord
  #   opt.val[[i]] <- optres$value
  #   opt.gra[[i]] <- optres$counts
  #   opt.con[[i]] <- optres$convergence
  #   opt.mes[[i]] <- ifelse(is.null(optres$message),
  #     "NULL", optres$message)
  #   opt.hes[[i]] <- optres$hessian
  # }
  # res <- list(par = opt.par, value = opt.val, counts = opt.gra,
  #   convergence = opt.con, message = opt.mes, hessian = opt.hes)
  # res
  # }

  optim.sumexp.sig <- function(data, oral = F, nexp = 3) {
  # Determines best sum of exponential for given data
  # data = mean pooled data;
  #        first column independent variable;
  #        second column dependent variable
  # oral = whether the drug displays an absorption curve
  # nexp = maximum fitted number of exponentials
  # Prepare data (remove t == 0, remove negative concentrations)
    x <- data[which(data[, 2] > 0 & data[, 1] != 0), 1]
    y <- data[which(data[, 2] > 0 & data[, 1] != 0), 2]
  # Set up objects
    res <- list(par = NULL, sumexp = NULL, value = NULL,
      error = NULL, hessian = NULL, message = NULL)
    lm.sub <- which(y == max(y))[1]:length(y)
    repeat {
      lm.mod <- lm(log(y[lm.sub]) ~ x[lm.sub])
      lmres <- unname(lm.mod$coefficients)
      if (is.na(lmres[2])) {
        lmres <- c(max(y), -0.00001)
        break
      }
      if (lmres[2] < 0) break
      else lm.sub <- tail(lm.sub, -1)
    }
  # Estimate candidate model parameters
    lm.sd <- sd(residuals(lm.mod))
    lm.add <- matrix(c(0, max(y)*lm.sd), nrow = 2)
    lm.prop <- matrix(c(lm.sd/50, lm.sd*50), nrow = 2)
    cand.mod <- expand.grid(1:nexp, c("add", "prop", "both"))  # candidate models
    for (i in 1:nrow(cand.mod)) {
      mod <- cand.mod[i, ]
      if (mod[[2]] == "both") lm.err <- cbind(lm.add, lm.prop)
      else lm.err <- get(paste0("lm.", mod[[2]]))
      gares <- try(ga("real-valued",
        mle.sumexp.sig, x = x, y = y, ga = T, errmod = mod[[2]],
        min = c(rep(lmres[2]*50, mod[[1]] + oral), rep(lmres[1]-2, mod[[1]]), lm.err[1,]),
        max = c(rep(lmres[2]/50, mod[[1]] + oral), rep(lmres[1]+2, mod[[1]]), lm.err[2,]),
        selection = gareal_lrSelection,
        crossover = gareal_spCrossover,
        mutation = gareal_raMutation,
        maxiter = 50,
        popSize = 250,
        monitor = F
      ))
      if (class(gares) == "try-error") browser()
      optres <- try(optim(
        gares@solution[1, ],
        mle.sumexp.sig,
        method = "BFGS", hessian = T,
        x = x, y = y, errmod = mod[[2]]
      ))
      if (class(optres) == "try-error") {
        optres <- list(
          par = gares@solution[1,],
          value = mle.sumexp.sig(gares@solution[1,], x, y, errmod = mod[[2]]),
          counts = c("function" = 501, gradient = NA),
          convergence = 99,
          message = "zero gradient",
          hessian = matrix(NA,
            ncol = length(gares@solution[1,]),
            nrow = length(gares@solution[1,])
          )
        )
      }
    # Create output
      fit.par <- head(optres$par, -ncol(lm.err))
      err.par <- tail(optres$par, ncol(lm.err))
      par.ord <- order.sumexp(fit.par, mod[[1]], oral)
      res$par[[i]] <- optres$par
      res$sumexp[[i]] <- par.ord
      res$value[[i]] <- c(ofv = optres$value, optres$counts)
      res$error[[i]] <- c(signif(err.par, 5), type = paste(mod[[2]]))
      res$hessian[[i]] <- optres$hessian
      res$message[[i]] <- c(convergence = optres$convergence,
        message = ifelse(is.null(optres$message), "NULL", optres$message)
      )
    }
    res
  }


# -----------------------------------------------------------------------------
# Set up
  data <- d1a
  par <- d1a.p
  t1 <- d1a.t

  niter <- dim(data)[2]
  absorp <- ifelse((dim(par)[1] %% 2) != 0, T, F)
  err <- matrix(
    1 + rnorm(n = length(t1)*niter, mean = 0, sd = 0.05),
    nrow = length(t1), ncol = niter
  )
  subd <- data[which(time.samp %in% t1),]*err
  auc.tlast.true <- apply(par, 2, function(x) pred.tlast.lam(x))


# Original
  set.seed(256256)
  res.sumexp <- apply(subd, 2, function(x) {
    out <- optim.sumexp.hes(
      data.frame(time = t1, conc = x), oral = absorp
    )
    print("done")
    chisq.sumexp.aic(out)
  })
  fit.par.hes <- lapply(res.sumexp, FUN = function(x) {
    x$par
  })
  auc.tlast.hes <- sapply(fit.par, function(x) pred.tlast.lam(x))

# New Method - Original
  set.seed(256256)
  res.sumexp <- apply(subd, 2, function(x) {
    out <- optim.sumexp.new(
      data.frame(time = t1, conc = x), oral = absorp
    )
    print("done")
    best.sumexp.aic(out)
  })
  fit.par.new <- lapply(res.sumexp, FUN = function(x) {
    x$sumexp
  })
  auc.tlast.new <- sapply(fit.par, function(x) pred.tlast.lam(x))

# Estimating standard deviation
  set.seed(256256)
  i <- 1
  res.sumexp <- apply(subd, 2, function(x) {
    if (i == 14) browser()
    out <- optim.sumexp.sig(
      data.frame(time = t1, conc = x), oral = absorp
    )
    print("done")
    i <<- i + 1
    best.sumexp.aic(out)
  })
  # problem with the 14th dataset, returning error
  # argument "errmod" is missing, with no default
  fit.par <- lapply(res.sumexp, FUN = function(x) {
    x$sumexp
  })
  auc.tlast.sig <- sapply(fit.par, function(x) pred.tlast.lam(x))

# New Method - BIC
  set.seed(256256)
  res.sumexp <- apply(subd, 2, function(x) {
    out <- optim.sumexp.new(
      data.frame(time = t1, conc = x), oral = absorp
    )
    print("done")
    best.sumexp.bic(out, length(t1))
  })
  fit.par.new <- lapply(res.sumexp, FUN = function(x) {
    x$sumexp
  })
  auc.tlast.new.bic <- sapply(fit.par, function(x) pred.tlast.lam(x))

# New Method - LRT
  set.seed(256256)
  res.sumexp <- apply(subd, 2, function(x) {
    out <- optim.sumexp.new(
      data.frame(time = t1, conc = x), oral = absorp
    )
    print("done")
    best.sumexp.lrt(out)
  })
  fit.par.new <- lapply(res.sumexp, FUN = function(x) {
    x$sumexp
  })
  auc.tlast.new.lrt <- sapply(fit.par, function(x) pred.tlast.lam(x))
