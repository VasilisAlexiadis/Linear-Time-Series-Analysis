# Load required packages
library(forecast)
library(fpp2)
library(signal)
library(stats)
install.packages("devtools")
library(tseries)
library(ggplot2)
install.packages("itsmr")
library(itsmr)
#################################################################################################################
# Load the time series dataset
x<-data(qgas)
x
summary(qgas )
str(qgas)
# Check for missing values in the time series
missing_values <- sum(is.na(qgas))
missing_values
# Plot the original time series
plot(qgas , main = "Original Time Series")

#################################################################################################################

# here we can see the periodicity
acf(qgas,20)
# From the acf plot we see that there is a 4 time-unit period

#################################################################################################################
# Convert the "qgas" dataset into a time series with a proper time index
 qgas <- ts(qgas, start = 1956, frequency = 4) 

# Plot the converted time series to check consistency
plot(qgas , main = "Converted Time Series")

acf(qgas,40)

#################################################################################################################

# Slice the time series to start from 1970
qgas <- window(qgas, start = c(1970))

# Plot the sliced time series
plot(qgas, main = "Sliced Time Series: qgas")


#################################################################################################################

# Taking the natural logaritm
qgas_log<-log(qgas)
# Plot the original time series
 plot(qgas_log , main = "Logarithmic Time Series")
# We notice that there is distortion and induces non linearity in trend with variance not much improved.
 
 #################################################################################################################

# here we can see the periodicity
acf(qgas)
# From the acf plot we see that there is a 4 time-unit period



 #################################################################################################################
# DECOMPOSITION
 #################################################################################################################
# First decomposition method
 #################################################################################################################

# We use the function decompose
decomp <- decompose(qgas)
# Plot the decomposition
plot(decomp)
# Extract the seasonal component from the decomposition
seasonal_component <- decomp$seasonal
# Plot the seasonal component
plot(seasonal_component, main="Seasonal Component")
# Extract the trend component from the decomposition
trend_component<- decomp$trend
# Plot the seasonal component
plot(trend_component, main="Trend Component")
# Extract the random component from the decomposition
random_component<- decomp$random
# Plot the random component
plot(random_component, main="Random Component")
plotc(random_component)
# seasonality adjasted
qgas_adj<-qgas-seasonal_component
plot(qgas_adj, main="Seasonality Adjusted")

diff_qgas_adj<-diff(qgas_adj,)
plot(diff_qgas_adj)

# Check for missing values
missing_values <- sum(is.na(random_component))
# If there are missing values, you can either remove them or impute them
# Removing missing values
random_component <- random_component[!is.na(random_component)]
adf.test(random_component)
# Now you can apply the acf() function to the updated random_component time series
pacf_values <- pacf(random_component)
print(pacf_values)
ggtsdisplay(random_component)

adf.test(random_component, k=5)

#################################################################################################################
# Second decomposition method
#################################################################################################################


# Extract the seasonal component 
seasonal_component1<-season(qgas,4) 
# Plot the seasonal component 
plot(seasonal_component1, main='Seasonal Component1')
plotc(seasonal_component1)
# ts without seasonality
qgas1<-qgas-seasonal_component1
# Plot ts without seasonality
plot(qgas1, main='ts without seasonality')
# Extract the trend component 
trend_component1<-trend(qgas1,8) # I used argument firs2 because the trend is nearly linear
# Plot trend component
plot(trend_component1, main='Trend component1')
# ts without seasonality and trend
qgas2<-qgas1-trend_component1
# Plot ts without seasonality and trend = random_component
plot(qgas2, main='ts without seasonality and trend')
plotc(qgas2)
# Deseasonilize and detrend = random_component1 = residuals
random_component1<- qgas-seasonal_component1-trend_component1
# Plot the random component1
plot(random_component1, main="Random Component1")
# plotc(random_component1)

test(random_component1)

pacf(random_component1)

adf.test(random_component1)
ggtsdisplay(random_component1)

#################################################################################################################
# Third decomposition method / Differencing
#################################################################################################################

# Differencing of degree 4 since our seasonality is of time units 4
diff_qgas<-diff(qgas, 4)
# Plot differenced ts
plot(diff_qgas, main="Differenced ts")
test(diff_qgas)

adf.test(diff_qgas)
ggtsdisplay(diff_qgas)


#################################################################################################################
# SEARCHING BEST ARMA MODEL
#################################################################################################################

model1 <- NULL
best_aicc <- Inf

# Iterate over different lag orders
for (p in 1:10) {
  for (q in 0:15) {
    # Estimate AR model
    model <- arma(random_component1, p = p, q = q) 
    
    # Retrieve the AIC value
    aicc <- model$aicc 
    
    # Check if the current model has a lower AIC than the previous best model
    if (aicc < best_aicc) {
      model1 <- model
      best_aicc <- aicc
    }
  }
}


# Print the best ARMA model
print(model1)


#################################################################################################################
# AUTOFIT
#################################################################################################################

model2<-autofit(random_component1, p = 0:10, q = 0:15)
model2


#################################################################################################################
# Burg's Algorithm
#################################################################################################################
model3 <- NULL
best_aicc <- Inf

# Iterate over different lag orders
for (p in 1:20) {
  # Estimate AR model using arburg()
  model <- burg(random_component1, p = p)
  
  # Retrieve the AIC value
  aicc <- model$aicc
  
  # Check if current model has lower AIC than the previous best model
  if (aicc < best_aicc) {
    model3 <- model
    best_p<-p
    best_aicc <- aicc
  }
}
print(model3)
print(best_p)



#################################################################################################################
# Tule Walker's Algorithm
#################################################################################################################

model4 <- NULL
best_aicc <- Inf

# Iterate over different lag orders
for (p in 1:20) {
  # Estimate AR model using arburg()
  model <- yw(random_component1, p = p)
  
  # Retrieve the AIC value
  aicc <- model$aicc
  
  # Check if current model has lower AIC than the previous best model
  if (aicc < best_aicc) {
    model4 <- model
    best_p<-p
    best_aicc <- aicc
  }
}
print(model4)
print(best_p)

#################################################################################################################
# Best Model residual test
#################################################################################################################

residual <- Resid(random_component1, M=NULL, a= model1)
test(residual)
adf_test <- adf.test(residual)
print(adf_test)

#################################################################################################################
# Forecast
#################################################################################################################
# A.
#################################################################################################################
# Set the size of the graphics device
#options(repr.plot.width = 20, repr.plot.height = 20)

# Generate the forecast plot
M <- c(NULL, "season", 4, "trend", 8)
x<-forecast(qgas, M = M, a = model1, h = 10, opt = 1,alpha=0.05)

qgas_forecasted <- append(qgas, x$pred)
plotc(qgas_forecasted,qgas)
#################################################################################################################
# B.
arar(qgas)
#################################################################################################################
# Fitting arima model
#################################################################################################################

# Fit the optimal ARIMA model using 'auto.arima()'
model5 <- auto.arima(random_component1)

# Print the fitted ARIMA model
print(model5)
#################################################################################################################
# Fitting inovation algorithm
#################################################################################################################


a = ia(random_component1,9)
print(a)