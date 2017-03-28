#################
# This script forms the foundation for building a multivariate dynamic regression 
# model for timeseries data by building a regression model using ARIMA errors. The use
# of AR errors corrects the regression model for autocorrelation in the data.
#
# Data needs to be a ts() object for the Arima functions. 
#################

#################
# Load libs
library(forecast)
library(tseries)

#################
# Define functions

num_diffs <- function(df, a, t){
  apply(df, 2, function(x){ndiffs(x, alpha = a, test = t)})
}

num_sdiffs <- function(df, f){
  apply(df, 2, function(x){nsdiffs(ts(x, frequency = f))})
  }

apply_same_diffs <- function(df, f, s, a = 0.05, t = c("adf"), d=FALSE, D=FALSE){
  # This function applies seasonal and standard differencing using passed parameters or as calculated by nsdiffs() and 
  # ndiffs(). If parameters passed, all series in df will be differenced the same. Otherwise each series will be 
  # differenced according to the maximum level of differencing as returned by nsdiffs()/ndiffs()
  #
  # All series will have the same level of differencing. 
  #
  # df is ts dataframe
  # f is frequency of data (monthly = 12, quarterly = 4, etc.)
  # s is start arg for ts() (ex, c(1997,1) for monthly data starting Jan 1997)
  # a is alpha level for ndiffs() 
  # t is test for ndiffs()
  # d is level of normal differencing to be applied
  # D is level of seasonal differencing to be applied
  
  if(d & D > 0) {
    diff_data <- apply(df, 2, function(x) {diff(x, lag = f, differences = D)})
  } else if (d) {
    diff_data <- apply(df, 2, function(X) {diff(x, differences = d)})
  } else if (D > 0){
    sdiff_data <- apply(df, 2, function(x) {diff(x, lag = f, differences = D)})
    diff_data <- apply(sdiff_data, 2, function(X) {diff(x, differences = d)})    
  } else { 
    
    # First determine max order of Seasonal Differencing and difference individual series if needed. 
    max_sdiff <- max(num_sdiffs(df, f))

    # Seasonal differencing. NOTE: If done this will remove the first year from data.
    # Need to include freq (f) and cast as timeseries in order for nsdiffs to work. 
    if(max_sdiff > 0) {
      sdiff_data <- ts(as.data.frame(apply(df, 2, function(x) 
        {diff(x, lag = f, differences = max_sdiff)})), 
        frequency = f, start = s)
    } else {
      sdiff_data <- df 
    }
    
    # Second, determine if any additional differencing must be done.  
    max_diff <- max(num_diffs(sdiff_data, a, t))

    # Apply differencing and subset to ensure data of the same length. NOTE: This assumes that 
    # max_diff <= 2. If it is larger the function won't align the data correctly. 
    if(max_diff > 0){
      diff_data <- ts(as.data.frame(apply(sdiff_data, 2, function(x) 
        {diff(x, differences = max_diff)})), 
        frequency = f, start = s)
    }
  }
  
  cat(c("D = ", max_sdiff)) # Print out max_sdiff for reference 
  if(max_sdiff > 0){cat("nNOTE: The returned series is shorter than the original series because of seasonal differencing.")}
  cat(c("\nd = ", max_diff)) # Print out max_diff for reference
  
  return(diff_data)
}

apply_diffs <- function(df, a, t, f, d=FALSE, D=FALSE){
  # This function applies seasonal and standard differencing using passed parameters or as 
  # calculated by ndiffs() and nsdiffs(). If parameters passed, all series in df will be differenced
  # the same. Otherwise each series will be differenced according to the ndiffs()/nsdiffs() results.
  # ** This may result in each series having a different level of differencing **
  
  if(d & D > 0) {
    diff_data <- apply(df, 2, function(x) {diff(x, lag = f, differences = D)})
  } else if (d) {
    diff_data <- apply(df, 2, function(X) {diff(x, differences = d)})
  } else if (D > 0){
    sdiff_data <- apply(df, 2, function(x) {diff(x, lag = f, differences = D)})
    diff_data <- apply(sdiff_data, 2, function(X) {diff(x, differences = d)})    
  } else { 
    
    # First determine max order of Seasonal Differencing and difference individual series if needed. 
    max_sdiff <- max(num_sdiffs(df, f))
    
    # Seasonal differencing. NOTE: If done this will remove the first year from data.
    # Need to include freq (f) and cast as timeseries in order for nsdiffs to work. 
    sdiff_data <- ts(as.data.frame(apply(df, 2, function(x) 
    {if (nsdiffs(ts(x, frequency = f)) <= 0) 
      tail(x, -max_sdiff*f)
      else if (nsdiffs(ts(x, frequency = f)) == max_sdiff) 
        diff(x, lag = f, differences = nsdiffs(ts(x, frequency = f)))
      else if (nsdiffs(ts(x, frequency = f)) == (max_sdiff - 1))
        diff(x, lag=f, differences = nsdiffs(ts(x, frequency = f)))})), 
    frequency = f, start = c(1997,1))
    ### HOW TO PROGRAMMATICALLY DETERMINE START FOR ABOVE??
    
    # Second, determine if any additional differencing must be done.  
    max_diff <- max(num_diffs(sdiff_data, a, t))
    
    # Apply differencing and subset to ensure data of the same length. NOTE: This assumes that 
    # max_diff <= 2. If it is larger the function won't align the data correctly. 
    diff_data <- ts(as.data.frame(apply(sdiff_data, 2, function(x) 
    {if (ndiffs(x, alpha = a, test=t) == max_diff) 
      diff(x, differences = ndiffs(x, alpha = a, test = t)) 
      else if (ndiffs(x, alpha = a, test = t) == (max_diff - 1))
        tail(diff(x, differences = ndiffs(x, alpha = a, test = t)), -max_diff+1)
      else tail(x, -max_diff)})), 
    frequency = f, start = c(1997,2))
    ### HOW TO PROGRAMMATICALLY DETERMINE START FOR ABOVE??
  }
  
  cat("D = ", max_sdiff) # Print out max_sdiff for reference  
  cat("d = ", max_diff) # Print out max_diff for reference
  
  return(diff_data)
}

cointegration<-function(x, y){
  vals<-data.frame(x,y)
  beta<-coef(lm(vals[,2]~vals[,1]+0,data=vals))[1]
  (adf.test(vals[,2]-beta*vals[,1], alternative="stationary", k=0))$p.value
}

coint_matrix <- function(df){
  # This function takes a dataframe as input and outputs a matrix comparing pairwise cointegration of 
  # variables. The value in the matrix is the p-value from the adf test in the cointegration() function. 
  #
  # df is a dataframe holding the series to be compaired. All pairs of variables will be examined. 
  
  m <- matrix(NA, ncol=NCOL(df), nrow = NCOL(df))
  colnames(m) <- colnames(df)
  rownames(m) <- colnames(df)
  
  for(n in 1:(NCOL(df)-1)){
    for(i in (n+1):NCOL(df)){
      m[n,i] <- cointegration(df[,n], df[,i])
    }
  }
  
  return(m)
}

getrmse <- function(x, y, h, ...){
  train.end <- time(x)[length(x)-h]
  test.start <- time(x)[length(x)-h+1]
  train.x <- window(x, end=train.end)
  test.x <- window(x, start=test.start)
  train.y <- window(y, end=train.end)
  test.y <- window(y, start=test.start)
  
  fit <- Arima(train.x, xreg = train.y, ...)
  fc <- forecast(fit, xreg = test.y, h=h)
  plot(fc)
  lines(x, col="red")
  #print(fc$mean)
  return(accuracy(fc, test.x)[2, "RMSE"])
}

revBoxCox <- function(x, lambda){
  if(lambda==0){
    out <- sapply(x, function(x1){exp(x1)})
  } else {
    out <- sapply(x, function(x1){(abs(lambda * x1 + 1))^(1/lambda)})
  }
  return(out)
}

revDiff <- function(x, y, d, D=0, f=0){
  # x is original data set (need full set for seasonal diffed data)
  # y is differenced data set
  # d is normal differencing factor
  # D is seasonal differencing factor
  # f is seasonal frequency
  
  if(d > 0){
    while(d > 0){
      y <- c(x[1], x[1] + cumsum(y))
      d = d -1
    }
  }
  
  if(D){
    if (D > 1){
      stop("Function not implemented for D > 1")
    } else {
      for(i in 1:NROW(y)){
        y[i] <- y[i] + x[i]
      }
      y <- c(x[1:f], y)
    }
  }
  
  return(y)
  
}

lm_pvalue <- function (modelobject) {
  # Returns p-value for a lm() model.
  # From : http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
  
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

build_lm_arima <- function(x, y){
  # This method calculates several statistics for each variable in x against the dependent variable y. 
  # It performs simple linear regression to calculate p-value, f statistic, r squared, and adjusted r squared;
  # utilizes auto.arima to get AICc, and uses the getrmse() function to get RMSE.
  # 
  # This function one regressor for each test and does not account for any interaction among regressors.
  #
  # x is a dataframe of potential regressors
  # y is the dependent variable

  
  m <- matrix(NA, ncol = 6, nrow = NCOL(x))
  rownames(m) <- colnames(x)
  colnames(m) <- c("p-value", "fstatistic", "r.squared", "adj.r.squared", "AICc", "RMSE")
  
  for(i in 1:NCOL(x)){
    reg_model <- lm(y ~ x[,i])
    fit <- auto.arima(y, xreg = x[,i])
    rmse <- getrmse(x[,i], y, h=24, order=arimaorder(fit)[1:3], seasonal = arimaorder(fit)[4:6])
    
    m[i,1] <- lm_pvalue(reg_model)
    m[i,2] <- summary(reg_model)$fstatistic[1]
    m[i,3] <- summary(reg_model)$r.squared
    m[i,4] <- summary(reg_model)$adj.r.squared
    m[i,5] <- fit$aicc
    m[i,6] <- rmse
  }
  
  return(m)
}

build_ccf_plots <- function(x, y, filepath){
  # This function takes the dependent variable and a list of possible regressors to build a ccf plot. 
  # The plots will be saved as a pdf at the filepath location specified.
  #
  # x is a dataframe of potential regressors
  # y is the dependent variable
  # filepath is the location where the pdf should be saved (including file name)
  
  pdf(file=filepath)
  par(mfrow=c(3,1))
  
  for(i in 1:NCOL(x)){
    if(!is.null(colnames(x)[i])){
      if(!is.null(colnames(y))){
        ccf(x[,i], y, type="correlation", main=c(colnames(x)[i], " & ", colnames(y)))
      } else {
        ccf(x[,i], y, type="correlation", main=c(colnames(x)[i], " & Y"))
      } 
    } else {
      ccf(x[,i], y, type="correlation", main = "X & Y")
    }
  }
  
  #dev.off()
  
}

test_lags <- function(x, y, lags){
  # This method takes a dependent variable and a single regressor and builds an Arima model for each 
  # lag specified in the lags argument. Lags apply to the regressor should be passed like c(-1, -3, -5)
  #
  # y is the dependent variable
  # x is the regressor
  # lags is a list of lags to be tested
  
  m <- matrix(NA, nrow=length(lags), ncol = 5)
  colnames(m) <- c("AICc", "AIC", "BIC", "Sigma2", "Log Likelihood")
  rownames(m) <- lags
  
  for(i in 1:length(lags)){
    xx <- lag(x, lags[i])
    fit <- auto.arima(y[(abs(lags[i])+1):length(y)], xreg=xx[1:(length(xx)+lags[i])], d = 0)
    m[i, 1] <- fit$aicc
    m[i, 2] <- fit$aic
    m[i, 3] <- fit$bic
    m[i, 4] <- fit$sigma2
    m[i, 5] <- fit$loglik
  } 
  
  return(m)
}

find_max_ccf <- function(a, b, f){
  d <- ccf(ts(a, frequency = f), ts(b, frequency = f), plot=FALSE, lag.max = length(a)-5)
  cor = d$acf[,,1]
  abscor = abs(d$acf[,,1])
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  absres = data.frame(abscor,lag)
  absres_max = res[which.max(absres$abscor),]
  return(absres_max$cor)
}

#################
# Load data : Assumes that all series are the same frequency and length.
# Data should have already been cleaned, missing values filled, etc. 

raw_df <- read.csv()
df_ts <- ts(raw_df, frequency = 12, start = c(1996,1))

#################
# Transform using BoxCox 
#QUESTION: Apply to all series or only those that need it? How to keep track of all lambdas?

  # I'm going to skip this for the time being. 
  # lambda <- BoxCox.lambda(x)
  # BoxCox(x, lambda)
  # InvBoxCox(x, lambda) # Reverses BoxCox

#################
# Ensure data is stationary. Difference to same level if needed (using apply_same_diffs())
# Print out diff levels for reference and assign to variable (need it for revDiff() later)

n_diff <- max(num_diffs(df_ts, 0.05, c("adf"))); n_diff
n_sdiff <- max(num_sdiffs(df_ts, f = 12)); n_sdiff

df_diff <- apply_same_diffs(df_ts, f = 12, s = c(1996,1))

x <- df_diff[] # regressors
y <- df_diff[] # dependent variable

#################
# Determine which regressors have correlation and should be included in the final model.
# Three tests will be used: 
#     - a general lm() and auto.arima() to determine p-value and R2 (these will be the 
#       same regardless of the lag used), AICc, and RMSE 
#     - a cointegration test to determine which pairs of variables are cointegrated
#     - a ccf plot to determine which lags to test

# Generate lm(), auto.arima(), and getrmse() for each and output relavant stats to a csv
lm_arima_matrix <- build_lm_arima(x, y)
write.csv(lm_arima_matrix, file = "C:/Users/Zach.Bremmer/Desktop/lm_values.csv")

# Test for cointegraion among regressors using raw data and output results to csv
m <- coint_matrix(df_ts) 
write.csv(m, file = "C:/users/Zach.Bremmer/Desktop/cointegration.csv")

# Build ccf plots and output to file (dev.off() doesn't work in the function)
build_ccf_plots(x, y, "C:/Users/Zach.Bremmer/Desktop/ccf_plots.pdf"); dev.off();

### TESTED UP TO THIS POINT, WORKS WELL 

#################
#
# Using results of lm_arima_matrix, ccf and coint_matrix(), user needs to manually select the regressors 
# to use in final model.
#
# The variable selection process depends on the final goal for the model. If accuracy is the goal, values for R2
# and their p-values should be the primary criteria. If parsimony is the key goal, AICc should be primary criteria.
# In most cases they will be the same variables, but this isn't always the case. RMSE is also included to give an 
# indication of forecast errors using each individual regressor. 
#
# Once the initial regressors are chosen, determine if any are cointegrated (p < 0.05). If cointegrated, both cannot be
# included in the model. Instead choose the one that has the highest R2/lowest p-value OR the one that will be easier 
# to forecast for future periods. 
#
# Once final regressors are choosen, examine the ccf plots to determine which lags to test to include in 
# the final model. 
#
# Note that this selection process does not account for any interactions between regressors. 
#
#################


#################
# Test lags of predictor variables as indicated on ccf() results. This will need to be performed for each regressor.
# For each regressor and lag combo an Arima model will be built and summary stats returned. The lag with the lowest 
# AICc score is the one that should be used. 

m1 <- test_lags(x1, y, c(-1, -3, -5))
m2 <- test_lags(x2, y, c(-3, -4, -5))

#################
# Using the results above, build out an Arima model with the proper regressors and lags. We can also do cross validation 
# by training on a subset of the data and testing on the remaining values. An RMSE score will be returned by the 
# function and a plot of the forecast vs. actual value of y will be created.

regressors <- df_diff[,c("", "")] # need to adjust for lags too
y_adj <- y[] # adjust y for lags

fit <- auto.arima(y_adj, xreg = regressors, d = 0)
summary(fit)

#################
# Check residuals to ensure no autocorrelation. 
tsdisplay(arima.errors(fit), main = "ARIMA Errors")
Box.test(residuals(fit), fitdf=5, lag=10, type="Ljung") # There is correlation if p < 0.05
hist(fit$residuals)
acf(fit$residuals)

#################
# Cross validate with test data 

getrmse(x, y, h=24, order=c(3,0,3)) # Get order from summary(fit) above

#################
# Once we have a good model, we need to get the forecast and untransform those values.
# This includes reversing any differencing that was done and undoing BoxCox if done. 

# Forecast values can be retreived from the forecast(fit) model, and then transformed 
# using InvBoxCox() or revDiff()


fcast <- forecast(fit, h = 12, xreg = )
fcast$method # Gives ARIMA model used for forecast
fcast$mean # Will give point prediction for each forecast 1:h
fcast$upper # Gives upper CI level for each prediction - one col for 80% CI and one for 95% CI
fcast$lower # Gives lower CI level for each prediction - one col for 80% CI and one for 95% CI
plot(fcast$residuals) # May be useful for something. acf()?


# What can I do to get a reliable/intelligible coefficient for each regressor? At the very least I need
# to be able to quantify the strength and direction of the relationship. If I pull the ARIMA coeff and 
# transform (revDiff()) it does that make sense? I think it does 
# See http://stats.stackexchange.com/questions/23881/reproducing-arima-model-outside-r
