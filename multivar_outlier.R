## contains functions for multivariate outlier detection/removal
library(mvoutlier) # used for chi-square plots
library(robustbase)

###########################################
# Multivariate outlier detection using robust Mahalanobis distance
# calculated using minimum covariance determinent (MCD).
# Outliers are identified as points that cause the distribution
# of robust distances to deviate from a chi-squared distribution
# as follows: a linear model was fit to the ordered robust distances
# against the quantiles of a chi-squared distribution. Subsequently,
# outliers were identified by iteratively removing the point
# corresponding to the largest robust distance until the coefficient
# of the determination for the linear model passes a threshold.
# Inputs:
#        dat: a 2D matrix of size nxp, where n is the number of
#            observations and p is the number of variables.
#        adj.r.square.thresh: threshold on adjusted r-square on the
#            linear model.
#        plt: boolean value indicating whether to plot qq-plots
#        title: title of the qq-plots
# Output:
#        a boolean vector of length n, where TRUE indicates outlier
##########################################
outlier.mcd.qq_chisq <- function(dat, adj.r.square.thresh = 0.9, plt = F, plot_title = ''){
  mcd <- covMcd(dat, nsamp = "deterministic")
  cm <- mcd$center
  S <- mcd$cov
  # Take square root to get the distance
  d <- sqrt(apply(dat, 1, function(x) t(x-cm) %*% solve(S) %*% (x-cm)))
  
  outlier <- logical(length(d))
  
  x <- sort(d^2)
  y <- qchisq(ppoints(length(d)),df = ncol(dat))
  model <- lm(y ~ x)
  m <- summary(model)$coefficients['x','Estimate']
  r <- summary(model)$adj.r.squared
  if (plt){
    qqplot(x,y)
    abline(model$coefficients)
    title(paste(plot_title,'\n','slope=',m,'\n','adj.r.squared=',r))
  }
  while ((r<adj.r.square.thresh) & (sum(!is.na(d))>2)){
    outlier[which.max(d)] <- T
    d[which.max(d)] <- NA
    x <- sort(d^2)
    
    y <- qchisq(ppoints(sum(! is.na(d))),df = ncol(dat))
    model <- lm(y ~ x)
    m <- summary(model)$coefficients['x','Estimate']
    r <- summary(model)$adj.r.squared
    if (plt){
      qqplot(x,y)
      abline(model$coefficients)
      title(paste(plot_title,'\n','slope=',m,'\n','adj.r.squared=',r))
    }
  }
  return(outlier)
}


###########################################
# Multivariate outlier detection using robust Mahalanobis distance
# calculated using minimum covariance determinent (MCD).
# Points above a certain critical quantile of the chi-squared distribution
# are identified as outlier.
# Inputs:
#        dat: a 2D matrix of size nxp, where n is the number of
#            observations and p is the number of variables.
#        chisquare.crit.val: chi-squared critical value
#        plt: boolean value indicating whether to plot qq-plots
#        title: title of the qq-plots
# Output:
#        a boolean vector of length n, where TRUE indicates outlier
##########################################
rm.mcd_outlier <- function(dat, chisquare.crit.val = 0.975, plt = F, plot_title = ''){
  # multivariate outlier detection based on minimum covariance determinent (MCD),
  # i.e., robust Mahalanobis distance.
  # dat is a 2D matrix of size nxp, where n is the number of observations
  # and p is the number of variables
  # output is a boolean vector of length n, where T indicates outlier
  mcd <- covMcd(dat)
  cm <- mcd$center
  S <- mcd$cov
  # Take square root to get the distance
  d_rd <- sqrt(apply(dat, 1, function(x) t(x-cm) %*% solve(S) %*% (x-cm)))
  #d_rd <- sqrt(mcd$mah)
  # calculate critical value
  crit = sqrt(qchisq(chisquare.crit.val,df = ncol(dat)))
  if (plt) {
    plot(d_rd)
    abline(a=crit, b = 0)
    title(plot_title)
  }
  return(d_rd > crit)
}


###########################################
# Multivariate outlier detection using Mahalanobis distance.
# Points above a certain critical quantile of the chi-squared distribution
# are identified as outlier.
# Inputs:
#        dat: a 2D matrix of size nxp, where n is the number of
#            observations and p is the number of variables.
#        chisquare.crit.val: chi-squared critical value
#        plt: boolean value indicating whether to plot qq-plots
#        title: title of the qq-plots
# Output:
#        a boolean vector of length n, where TRUE indicates outlier
##########################################
rm.md_outlier <- function(dat, chisquare.crit.val = 0.975, plt = F, plot_title = ''){
  # multivariate outlier detection based on Mahalanobis distance
  # dat is a 2D matrix of size nxp, where n is the number of observations
  # and p is the number of variables
  # output is a boolean vector of length n, where T indicates outlier
  cm <- colMeans(dat)
  S <- cov(dat)
  # Take square root to get the distance
  d_md <- sqrt(apply(dat, 1, function(x) t(x-cm) %*% solve(S) %*% (x-cm)))
  # calculate critical value
  crit = sqrt(qchisq(chisquare.crit.val,df = ncol(dat)))
  if (plt) {
    plot(d_md)
    abline(a=crit, b = 0)
    title(plot_title)
  }
  return(d_md > crit)
}
