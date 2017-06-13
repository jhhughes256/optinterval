# DATA SUMMARY FUNCTIONS
# -----------------------------------------------------------------------------
# Convenience function for summarising a dataframe
  headtail <- function(x) {
    print(dim(x))
    print(head(x))
    print(tail(x))
  }

# 90% confidence interval functions
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)

# Define a function for geometric mean
  geomean <- function(x, na.rm = F) {
    if (na.rm == T) x <- x[is.na(x) == F]
    exp(mean(log(x)))
  }
  # Note x cannot be negative, zero

# Median and 90% tolerance intervals
  sumfunc90 <- function(x) {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  quantile(x, probs = 0.05, na.rm = T, names = F)  # 90%CI
    stat3 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    stat4 <-  length(na.omit(x))
    result <- c("median" = stat1, "low90" = stat2, "hi90" = stat3, "n" = stat4)
    result
  }

# Median and 95% tolerance intervals
  sumfunc95 <- function(x) {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  quantile(x, probs = 0.025, na.rm = T, names = F)  # 95%CI
    stat3 <-  quantile(x, probs = 0.975, na.rm = T, names = F)
    stat4 <-  length(na.omit(x))
    result <- c("median" = stat1, "low95" = stat2, "hi95" = stat3, "n" = stat4)
    result
  }

# Mean, sd and CV
  sumfuncCV <- function(x) {
    stat1 <-  mean(x, na.rm = T)
    stat2 <-  sd(x, na.rm = T)
    stat3 <-  stat2/stat1*100
    stat4 <-  length(na.omit(x))
    result <- c("mean" = stat1, "sd" = stat2, "cv" = stat3, "n" = stat4)
    result
  }

# Geomean, mean, sd, CV, min, max, n - for Millenium MLN 8237
  sumfuncMLN <- function(x) {
    stat1 <-  geomean(x, na.rm = T)
    stat2 <-  mean(x, na.rm = T)
    stat3 <-  sd(x, na.rm = T)
    stat4 <-  stat3/stat2*100
    stat5 <-  min(x, na.rm = T)
    stat6 <-  quantile(x, probs = 0.05, na.rm = T, names = F)  #90%CI
    stat7 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    stat8 <-  max(x, na.rm = T)
    stat9 <-  length(na.omit(x))
    result <- c("gmean" = stat1, "mean" = stat2, "sd" = stat3, "cv" = stat4,
      "min" = stat5, "lo90" = stat6, "hi90" = stat7, "max" = stat8, "n" = stat9)
    result
  }

# Geomean, mean, sd, CV, min, max, n - for CBIO
  sumfuncCBIO <- function(x) {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  mean(x, na.rm = T)
    stat3 <-  sd(x, na.rm = T)
    stat4 <-  stat3/stat2*100
    stat5 <-  min(x, na.rm = T)
    stat6 <-  quantile(x, probs = 0.05, na.rm = T, names = F)  #90%CI
    stat7 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    stat8 <-  max(x, na.rm = T)
    stat9 <-  length(na.omit(x))
    result <- c("median" = stat1, "mean" = stat2, "sd" = stat3, "cv" = stat4,
      "min" = stat5, "lo90" = stat6, "hi90" = stat7, "max" = stat8, "n" = stat9)
    result
  }

# Median, mean, sd, CV, 95%CI - for bootstrap parameter summary
  sumfuncBOOT <- function(x) {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  mean(x, na.rm = T)
    stat3 <-  sd(x, na.rm = T)
    stat4 <-  stat3/stat2*100
    stat5 <-  quantile(x, probs = 0.025, na.rm = T, names = F)  #95%CI
    stat6 <-  quantile(x, probs = 0.975, na.rm = T, names = F)
    result <- c("median" = stat1, "mean" = stat2, "sd" = stat3,
      "cv" = stat4, "lo95" = stat5, "hi95" = stat6)
    result
  }

# Summarize distribution by percentiles
  sumfuncPercentile <- function(x) {
    stat1 <-  quantile(x, probs = 0.05, na.rm = T, names = F)
    stat2 <-  quantile(x, probs = 0.10, na.rm = T, names = F)
    stat3 <-  quantile(x, probs = 0.25, na.rm = T, names = F)
    stat4 <-  quantile(x, probs = 0.50, na.rm = T, names = F)
    stat5 <-  quantile(x, probs = 0.75, na.rm = T, names = F)
    stat6 <-  quantile(x, probs = 0.90, na.rm = T, names = F)
    stat7 <-  quantile(x, probs = 0.95, na.rm = T, names = F)
    result <- c("05perct" = stat1, "10perct" = stat2,
      "25perct" = stat3, "50perct" = stat4, "75perct" = stat5,
      "90perct" = stat6, "95perct" = stat7)
    result
  }

#M ean, sd, min and max & n
  sumfuncRange <- function(x) {
    stat1 <-  mean(x, na.rm = T)
    stat2 <-  sd(x, na.rm = T)
    stat3 <-  min(x, na.rm = T)
    stat4 <-  max(x, na.rm = T)
    stat5 <-  length(na.omit(x))
    result <- c("mean" = stat1,"sd" = stat2, "min" = stat3,
      "max" = stat4, "n" = stat5)
    result
  }

# Median etc for boxplot
  sumfuncBOX <- function(x) {
    stat1 <-  median(x, na.rm = T)
    stat2 <-  quantile(x, probs = 0.025, na.rm = T, names = F)
    stat3 <-  quantile(x, probs = 0.25, na.rm = T, names = F)
    stat4 <-  quantile(x, probs = 0.75, na.rm = T, names = F)
    stat5 <-  quantile(x, probs = 0.975, na.rm = T, names = F)
    result <- c("median" = stat1, "q025" = stat2, "q25" = stat3,
      "q75" = stat4, "q975" = stat5)
    result
  }

# Define a function for geometric mean and 90% CI of the sem
  geomeansemCI <- function(x, na.rm = F) {
    logx <- log(x)
    logmean <- mean(logx)
    n <- length(x)
    logsem <- sd(logx)/sqrt(n)
    # Critical value of the t-distribution for two one-sided p=0.05
    critt <- qt(.95, df = (n-1))
    loglo95 <- logmean - critt*logsem
    loghi95 <- logmean + critt*logsem
    gmean <- exp(logmean)
    glo95 <- exp(loglo95)
    ghi95 <- exp(loghi95)
    result <- c("gmean" = gmean, "glo95" = glo95, "ghi95" = ghi95,
      "crit.t" = critt)
    result
  }
  # Note x cannot be negative, zero

# Computes the time of Cmax
  tmax <- function(dv, time) {
    cmax <- max(dv)
    tindex <- which(dv == cmax)
    tmax <- time[tindex]
    head(tmax, n = 1)  #as there can be 2 or more equal Cmax's, choose the first
  }
