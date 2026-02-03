############################
# Supplemental text material
############################

# compare forecast accuracy

# Does S1 NTD improve over cases alone?

library(forecast)

# load the data
load(file = "data/combined_data.Rdata")

# remove nas
dat_state <- 
  dat_state %>%
  filter(!is.na(case_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w))

# model 1
# baseline
model1 <- arima(dat_state$case_incidence, order=c(1,0,0))
summary(model1)

# model 2
# arima with predictor
model2 <- auto.arima(dat_state$case_incidence, xreg=dat_state$ntd_pi_state_3w)

summary(model2)
AIC(model1)
AIC(model2)
BIC(model1)
BIC(model2)

# ntd pi is better time series model
#Fit an AR(2) model to each rolling origin subset
far2 <- function(x, h){forecast(Arima(x, order=c(1,0,0)), h=h)}
e1 <- tsCV(dat_state$case_incidence, far2, h=1)

#Example with exogenous predictors
far2_xreg <- function(x, h, xreg, newxreg) {
  forecast(Arima(x, order=c(2,0,0), xreg=xreg), xreg=newxreg)
}

y <- dat_state$case_incidence
xreg <- dat_state$ntd_pi_state_3w
e2 <- tsCV(y, far2_xreg, h=1, xreg=xreg)

dm.test(e1, e2, h = 1)


# try with missing nas

# Combine into a data frame for easier NA handling
error_data <- data.frame(e1, e2)

# Remove rows with NAs from the data frame
error_data_clean <- na.omit(error_data)

# Extract the clean error vectors
e1_clean <- error_data_clean$e1
e2_clean <- error_data_clean$e2

# Now run the dm.test on the clean data
library(forecast)
dm.test(e1_clean, e2_clean, h = 1) 

rmse(model1)
rmse(model2)


#-----------------------------------
# hospitalization incidence
# ----------------------------------
############################
# Supplemental text material
############################

# compare forecast accuracy

# Does S1 NTD improve over cases alone?

library(forecast)

# load the data
load(file = "data/combined_data.Rdata")

# remove nas
dat_state <- 
  dat_state %>%
  filter(!is.na(hosp_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w))

# model 1
# baseline
model1 <- arima(dat_state$hosp_incidence, order=c(1,0,0))
summary(model1)

# model 2
# arima with predictor
model2 <- auto.arima(dat_state$hosp_incidence, xreg=dat_state$ntd_pi_state_3w)

summary(model2)
AIC(model1)
AIC(model2)
BIC(model1)
BIC(model2)

# ntd pi is better time series model
#Fit an AR(2) model to each rolling origin subset
far2 <- function(x, h){forecast(Arima(x, order=c(1,0,0)), h=h)}
e1 <- tsCV(dat_state$hosp_incidence, far2, h=1)

#Example with exogenous predictors
far2_xreg <- function(x, h, xreg, newxreg) {
  forecast(Arima(x, order=c(2,0,0), xreg=xreg), xreg=newxreg)
}

y <- dat_state$hosp_incidence
xreg <- dat_state$ntd_pi_state_3w
e2 <- tsCV(y, far2_xreg, h=1, xreg=xreg)

dm.test(e1, e2, h = 1)


# try with missing nas

# Combine into a data frame for easier NA handling
error_data <- data.frame(e1, e2)

# Remove rows with NAs from the data frame
error_data_clean <- na.omit(error_data)

# Extract the clean error vectors
e1_clean <- error_data_clean$e1
e2_clean <- error_data_clean$e2

# Now run the dm.test on the clean data
library(forecast)
dm.test(e1_clean, e2_clean, h = 1) 

rmse(model1)
rmse(model2)
