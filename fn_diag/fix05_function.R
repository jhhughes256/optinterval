# Function containing the framework setup in the intial study design
# -----------------------------------------------------------------------------
# Set up directories
  if (!exists("git.dir")) {
    rm(list = ls(all = T))
    wd <- c("C:/Users/Jim Hughes/Documents", "C:/Users/hugjh001/Documents",
      "C:/Users/hugjh001/Desktop", "C:/windows/system32")

    graphics.off()
    if (getwd() == wd[1]) {
      git.dir <- paste0(getwd(), "/GitRepos")
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
  #library(ggplot2)
  #theme_bw2 <- theme_set(theme_bw(base_size = 14))
  #theme_update(plot.title = element_text(hjust = 0.5))

# Source scripts to set up environment
  source(paste(git.dir, reponame, "fn_diag/fix_functions.R", sep = "/"))
  source(paste(git.dir, reponame, "study_data.R", sep = "/"))

# -----------------------------------------------------------------------------
# Set basic parameters
  data.names <- paste0("d", as.vector(outer(1:3, c("a", "b"), paste0)))
  par.names <- paste(data.names, "p", sep = ".")
  fn.names <- paste("pred", data.names, sep = ".")
  t2.names <- paste(data.names, "t", sep = ".")

  study.fn <- function(data, par, fn, nobs, t2, sdev = 0.05, tlast = 24, logauc = F) {
    absorp <- ifelse((length(par) %% 2) != 0, T, F)
    t1 <- round(seq(0, tlast, length.out = nobs), 1)
    e1 <- 1 + rnorm(n = length(t1), mean = 0, sd = sdev)
    subd <- data.frame(
      time = t1,
      conc = with(data, conc[time %in% t1])*e1
    )
    fit.par <- chisq.sumexp(optim.sumexp(subd, oral = absorp))$par
    int.t3 <- optim.interv.rep(fit.par, 2)
    t3 <- c(0, round(int.t3$par, 1), tlast)

    t4.tmax <- tmax.sumexp(fit.par)
    int.t4 <- optim.interv.rep(fit.par, 2, tmax = t4.tmax)
    nr.t4 <- c(0, round(int.t4$par, 1), tmax.sumexp(fit.par), tlast)
    t4 <- nr.t4[order(nr.t4)]

    fit.fn <- paste0("pred.d", floor(length(fit.par)/2), ifelse(absorp, "a", "b"))
    auc <- c(
      true = integrate(fn, 0, tlast, p = par)$value,
      sumexp = integrate(get(fit.fn), 0, tlast, p = fit.par)$value,
      basic = auc.interv(t1, par, fn, log = logauc),
      user = auc.interv(t2, par, fn, log = logauc),
      optint = auc.interv(t3, par, fn, log = logauc),
      optintwCmax = auc.interv(t4, par, fn, log = logauc)
    )
    cmax <- c(
      true = pred.sumexp(par, tmax.sumexp(par)),
      basic = max(pred.sumexp(par, t1)),
      user = max(pred.sumexp(par, t2)),
      optint = max(pred.sumexp(par, t3)),
      optintwCmax = max(pred.sumexp(par, t4))
    )
    tmax <- c(
      true = tmax.sumexp(par),
      basic = t1[which(pred.sumexp(par, t1) == cmax["basic"])],
      user = t2[which(pred.sumexp(par, t2) == cmax["user"])],
      optint = t3[which(pred.sumexp(par, t3) == cmax["optint"])],
      optintwCmax = t4[which(pred.sumexp(par, t4) == cmax["optintwCmax"])]
    )
    return(list(auc = auc, cmax = cmax, tmax = tmax))
  }

  fin.res <- list(NULL)
  for (i in 1:6) {
    fin.res[[i]] <- list(
      data = data.names[i],
      result = study.fn(get(data.names[i]),
        par = get(par.names[i]), fn = get(fn.names[i]),
        t2 = get(t2.names[i]), nobs = 9
      )
    )
  }
