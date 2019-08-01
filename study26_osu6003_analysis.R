# Analysing effect of adding/removing data points
# -----------------------------------------------------------------------------
# Enhancing the algorithm paper by determining the loss of accuracy/precision
#   when reducing the number of sample times. Aiming to prove that less time
#   points are needed.
# Data is from OSU6003; papers below:

# Blum W, Klisovic RB, Becker H, et al. Dose Escalation of Lenalidomide in 
#   Relapsed or Refractory Acute Leukemias. J Clin Oncol. 2010;28(33):4919-4925.
# Maddocks K, Ruppert AS, Browning R, et al. A dose escalation feasibility study
#   of lenalidomide for treatment of symptomatic, relapsed chronic lymphocytic 
#  leukemia. Leuk Res. 2014;38(9):1025-1029.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load packages
  library(GA)  # genetic algorithm (ga)
  library(plyr)  # iterative functions (ddply)

# Load data visualisation package (set theme, define palette)
  library(ggplot2)  # plots (ggsave etc.)
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  
# Source algorithm functions
  source("fn_diag/fix_functions.R")

# Define additional functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Last observation carried forward
  # Finds an NA and carries forward the previous value
  locf <- function(x) {
    good <- !is.na(x)
    positions <- seq(length(x))
    good.positions <- good * positions
    last.good.position <- cummax(good.positions)
    last.good.position[last.good.position == 0] <- NA
    x[last.good.position]
  }
  
# Determine AUC (regular)
  auc.regular <- function(times, C, log = F) {
    auc <- c(NULL)
    for (i in 2:length(C)) {
      h <- times[i] - times[i-1]
      dC <- C[i-1] - C[i]
      if (log & dC > 0) auc[i-1] <- dC*h/log(C[i-1]/C[i])
      else auc[i-1] <- (C[i-1] + C[i])*h/2
    }
    return(sum(auc))
  }
  auc.regular.lam <- function(t, dv) {
    lambz <- try(-pred.lambdaz(dv, t)[2])
    if (class(lambz) == "try-error") return(0)
    else return(tail(dv, 1)/lambz)
  }
    
  
# 2 comp. w/ abs sumexp
  pred.d2a <- function(x, p) {
    exp(p[1]*x + p[4]) + exp(p[2]*x + p[5]) - exp(p[3]*x + log(sum(exp(p[4]), exp(p[5]))))
  }
  
# Load/Tidy/Prepare the Data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source lenalidomide data
  nmprep <- read.csv(stringsAsFactors = F, 
    "C:/Users/Jim Hughes/Documents/GitRepos/LEN_PK/Data/data.csv"
  )
# Subset the data to look at OSU6003, without MDVs (excluding time = 0)
# Additionally the weekly troughs taken after the first 48 hours
  clindata <- nmprep[
    nmprep$STUDY == 6003 & 
    (nmprep$MDV == 0 | nmprep$TIME == 0) &
    nmprep$TIME <= 48,
  ]
  
# Check the times in the dataset
  table(clindata$TIME)
  # times are not uniform, must be binned to make mean pooled data
# Define nominal sample times
  samp.times <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 24)
  
# Create breaks to use with cut function
  # First find the mean of each pair of values to choose a break half way 
  # between each sample time
  inner.breaks <- rowMeans(cbind(head(samp.times, -1), tail(samp.times, -1)))
  # Then create first and last break points by extending them before and after
  # the final times by the first and final inner break points
  cut.breaks <- c(head(samp.times, 1) - head(inner.breaks, 1), 
    inner.breaks, 
    tail(samp.times, 1) + tail(inner.breaks, 1))
  # and now use cut
  bin.groups <- cut(clindata$TIME, cut.breaks)
  
# Use these bins to define the nominal time after dose 
# Remembering to convert from factor to numeric
  levels(bin.groups) <- samp.times
  clindata$TADNOM <- as.numeric(paste(bin.groups))
  
# Check nominal times in the dataset
  table(clindata$TADNOM)
  
# Normalise DV for dose size
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Check the AMT column
  table(clindata$AMT)
  # needs to be filled using last observation carried forward
  clindata$AMT[clindata$AMT == 0] <- NA
  clindata$DOSE <- locf(clindata$AMT)
  
# Check the DOSE column
  table(clindata$DOSE)
  head(clindata[clindata$TADNOM == 2,])
  # appears to have worked

# Normalise for dose of 25mg
  clindata$DOSENORM <- with(clindata, 25*DV/DOSE)
  
# Create mean pooled data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Use ddply to mean the normalised DV for each of the bins
  meanpool <- ddply(clindata, .(TADNOM), function(df) {
    out <- data.frame(MEAN = mean(df$DOSENORM))
  })
  
# Visualisation of mean pooled data
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + geom_line(aes(x = TADNOM, y = MEAN), data = meanpool, 
    col = cbPalette$red, size = 0.8)
  p1 <- p1 + geom_point(aes(x = TIME, y = DOSENORM), data = clindata,
    col = cbPalette$blue, size = 1, shape = 1, alpha = 0.7)
  p1  # looks pretty decent
  
# Fit sum of exponentials to mean pooled data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimise sum of exponential equation, then choose one with best aic
  sumexp.res <- best.sumexp.aic(optim.sumexp.new(meanpool, oral = T, nexp = 2))
  opt.sumexp <- sumexp.res$sumexp
  
# Visualisation of sum of exponentials
  # Create dataset
  sumexpdata <- data.frame(
    TIME = seq(0, 24, 0.1), 
    PRED = pred.sumexp(opt.sumexp, seq(0, 24, 0.1))
  )
  # Plot dataset against mean pooled
  p2 <- p1 + geom_line(aes(x = TIME, y = PRED), sumexpdata,
    col = cbPalette$blue, size = 0.8, alpha = 0.8)
  p2  # looks good!
  
# Determine optimal sample times and compare AUC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimise time intervals for various numbers of times (6-13)
# Output is list due to variable lengths of time vectors
  interv.res <- llply(6:13, function(i) {
    res <- optim.interv.dtmax(opt.sumexp, seq(0, 24, length.out = i))
    out <- res$times
  })

# Determine AUC0-24 of original mean pooled data
  out <- list(
    auc024 = data.frame(
      data = with(meanpool, auc.regular(TADNOM, MEAN)),
      sumexp = integrate(pred.d2a, 0, 24, p = opt.sumexp)$value,
      t6 = auc.interv.sumexp(interv.res[[1]], opt.sumexp),
      t7 = auc.interv.sumexp(interv.res[[2]], opt.sumexp),
      t8 = auc.interv.sumexp(interv.res[[3]], opt.sumexp),
      t9 = auc.interv.sumexp(interv.res[[4]], opt.sumexp),
      t10 = auc.interv.sumexp(interv.res[[5]], opt.sumexp),
      t11 = auc.interv.sumexp(interv.res[[6]], opt.sumexp),
      t12 = auc.interv.sumexp(interv.res[[7]], opt.sumexp),
      t13 = auc.interv.sumexp(interv.res[[8]], opt.sumexp)
    ),
    aucinf = data.frame(
      data = with(meanpool, 
        auc.regular(TADNOM, MEAN) + auc.regular.lam(TADNOM, MEAN)
      ),
      sumexp = integrate(pred.d2a, 0, 168, p = opt.sumexp)$value,
      t6 = auc.interv.sumexp(interv.res[[1]], opt.sumexp) 
        + auc.interv.lam(opt.sumexp, interv.res[[1]]),
      t7 = auc.interv.sumexp(interv.res[[2]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[2]]),
      t8 = auc.interv.sumexp(interv.res[[3]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[3]]),
      t9 = auc.interv.sumexp(interv.res[[4]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[4]]),
      t10 = auc.interv.sumexp(interv.res[[5]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[5]]),
      t11 = auc.interv.sumexp(interv.res[[6]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[6]]),
      t12 = auc.interv.sumexp(interv.res[[7]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[7]]),
      t13 = auc.interv.sumexp(interv.res[[8]], opt.sumexp)
        + auc.interv.lam(opt.sumexp, interv.res[[8]])
    ),
    cmax = data.frame(
      data = with(meanpool, max(MEAN)),
      sumexp = pred.sumexp(opt.sumexp, tmax.sumexp(opt.sumexp)),
      t6 = max(pred.sumexp(opt.sumexp, interv.res[[1]])),
      t7 = max(pred.sumexp(opt.sumexp, interv.res[[2]])),
      t8 = max(pred.sumexp(opt.sumexp, interv.res[[3]])),
      t9 = max(pred.sumexp(opt.sumexp, interv.res[[4]])),
      t10 = max(pred.sumexp(opt.sumexp, interv.res[[5]])),
      t11 = max(pred.sumexp(opt.sumexp, interv.res[[6]])),
      t12 = max(pred.sumexp(opt.sumexp, interv.res[[7]])),
      t13 = max(pred.sumexp(opt.sumexp, interv.res[[8]]))
    ),
    tmax = data.frame(
      data = with(meanpool, TADNOM[which(MEAN == max(MEAN))]),
      sumexp = tmax.sumexp(opt.sumexp),
      t6 = interv.res[[1]][which(
        pred.sumexp(opt.sumexp, interv.res[[1]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[1]]))
      )],
      t7 = interv.res[[2]][which(
        pred.sumexp(opt.sumexp, interv.res[[2]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[2]]))
      )],
      t8 = interv.res[[3]][which(
        pred.sumexp(opt.sumexp, interv.res[[3]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[3]]))
      )],
      t9 = interv.res[[4]][which(
        pred.sumexp(opt.sumexp, interv.res[[4]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[4]]))
      )],
      t10 = interv.res[[5]][which(
        pred.sumexp(opt.sumexp, interv.res[[5]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[5]]))
      )],
      t11 = interv.res[[6]][which(
        pred.sumexp(opt.sumexp, interv.res[[6]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[6]]))
      )],
      t12 = interv.res[[7]][which(
        pred.sumexp(opt.sumexp, interv.res[[7]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[7]]))
      )],
      t13 = interv.res[[8]][which(
        pred.sumexp(opt.sumexp, interv.res[[8]]) == 
        max(pred.sumexp(opt.sumexp, interv.res[[8]]))
      )]
    )
  )
  
  llply(out, function(x) {
    rbind(x, x/x$sumexp)
  })