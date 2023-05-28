library(ggplot2)
library(urca)
library(fpp3)
library(tidyverse)
library(lubridate)
library(tsibble)
library(seasonal)
library(ggfortify)
library(feasts)
library(strucchange)

####### Preliminary Steps #######

# Import raw data & shorten to retain only last 10 years
raw <- read.csv("tas_timeseries_monthly_cru_1901-2021_DNK.csv")
head(raw)
nrow(raw)
raw_short <- raw[112:121,]
nrow(raw_short)

# Transform data to obtain tsibble data
df <- data.frame()
for (i in 2:13){
  temp = raw_short[c(1,i)]
  temp$X <- paste(temp$X, i-1, sep="-")# %>% ym() %>% as.Date()
  colnames(temp) <- c("date", "temperature")
  df <- rbind(df, temp)  
}
df <- df%>%
  mutate("date" = yearmonth(date))%>%
  as_tsibble(index = date)
head(df)

# store also as Time series object
df_ts <- as.ts(df)
df_ts

# print summary statistics
summary(df_ts)

####### Exploratory Data Analysis & Preprocessing #######

# plot
df%>%autoplot(temperature) + labs(y= "Temperature(C°)", x="Month")# + ggtitle("Original plot before transformation") # include only during coding
# ACF
df%>%ACF%>%autoplot
# from both plots seasonality so strong that it's hard to see anything else --> therefore, look at decomposition

# STL decomposition (best suited if you expect seasonality)
dcmp_stl <- df %>%
  model(STL(temperature)) %>%
  components()
dcmp_stl %>% autoplot + xlab("Month") + ggtitle("STL Decomposition before transformation") # include only during coding
# comment on variations --> boxcox transform

# box cox --> makes forecasting easier
# The Box-Cox transformation transforms our data so that it closely resembles a normal distribution
# A good value of λ is one which makes the size of the seasonal variation about the same across the whole series, as that makes the forecasting model simpler.
lambda.g <- df%>%
  features(temperature, features = guerrero)%>%
  pull(lambda_guerrero)

# box cox transform data
df_bc <- df%>%
  mutate(temperature = box_cox(temperature, lambda.g))
head(df_bc)

df_bc%>%autoplot(temperature) + labs(y ="Box-Cox(Temperature in C°,1.72)",x = "Month") + ggtitle("Plot after transformation") # include only during coding

lambda.g
# --> value indicates power transformation

# STL decomposition
dcmp_stl_bc <- df_bc %>%
  model(STL(temperature)) %>%
  components()
dcmp_stl_bc %>% autoplot + xlab("Month")# + ggtitle("STL Decomposition after transformation") # include only during coding
# even less trend - cycle seem to fluctuate around constant mean

# Trend
summary(dcmp_stl_bc$trend)
dcmp_stl_bc%>%autoplot(trend)

df_trend_growth <- dcmp_stl_bc[c(2,4)] %>%
  mutate(growth_rate = trend - lag(trend))%>% #add column for growth rate
  drop_na(growth_rate) #dropping first entry as no growth rate computable
summary(df_trend_growth$growth_rate)
df_trend_growth%>%autoplot(growth_rate)
# unsteady but small from decomposition

df_growth%>%ACF%>%autoplot
# highly autocorrelated --> confirm with box pierce
df_growth%>%features(growth_rate, box_pierce, lag = 24, dof = 0)
# --> confirmation

# IF MORE PAGES NEEDED/ENOUGH SPACE - month wise
'monthly_growth <- raw_short
for (x in colnames(monthly_growth[,2:13])){
  monthly_growth%>%
    mutate()
}'

# OPTIONAL: Seasonality
df%>%gg_season(temperature)
# +Proof

# Stationarity
# ADF 
# Augmented Dickey-Fuller with drift as no clear trend observable
summary(ur.df(df_bc$temperature, type = "drift")) 
# both test statistics reject the 0 --> data potentially stationary --> test without drift

summary(ur.df(df_bc$temperature, type = "none")) 
# test statistics rejects the 0 --> data potentially stationary

# KPSS
summary(ur.kpss(df_bc$temperature, type = "mu"))
# can't reject the 0 --> data stationary

# Proof for taking zero difference (from feasts package):
unitroot_ndiffs(
  df_bc$temperature,
  alpha = 0.05,
  unitroot_fn = ~unitroot_kpss(.)["kpss_pvalue"],
  differences = 0:12
)

# structural breaks
df_bc_sb <- df_bc %>%
  mutate(
    Lag1 = lag(temperature)
  )

# QLR test with monthly lags
qlr <- Fstats(temperature ~ Lag1, data = df_bc_sb, from=0.15)
test <- sctest(qlr, type = "supF")
test

df_bc_sb%>%
  features(temperature, features=shift_level_max)

breakpoints(qlr, alpha = 0.05)

plot(qlr, alpha = 0.05, main = "QLR Test") + lines(breakpoints(qlr))
# This indicates that there are no structural breaks 

####### Model Creation #######

# train test split - FROM WHICH DATA
train_bc <- df_bc %>%
  filter_index(. ~ "2020 Jan")
test_bc <- df_bc %>%
  filter_index("2020 Jan" ~ .)

# ETS
# consistent error over time (doesn't depend on time level) --> additive 
# no consistent trend --> "None"
# consistent seasonality over time --> additive
models_1 <- train_bc %>%
  model(
    ets_guess = ETS(temperature ~ error("A") + trend("N") + season("A")),
    ets_auto = ETS(temperature),
    #ets_test = ETS(temperature ~ error("A") + trend("N") + season("M")),
  )

report(models_1%>%select(ets_guess))
#report(models_1%>%select(ets_test))
report(models_1%>%select(ets_auto))
# --> models are the same

# OPTIONAL
coef(models_1%>%select(ets_guess))

# Evaluate Residuals
models_1 %>%
  select(ets_guess) %>%
  gg_tsresiduals(type = "innovation")

# Ljung Box Test (we want HIGH p-value -> NULL not rejected, then resid. independent)
models_1 %>%
  select(ets_guess) %>%
  residuals()  %>%
  features(.resid, features = ljung_box, lag = 24, dof = 8)


# ARIMA
# as a standard rule, we apply seasonal differencing as we can see a clear seasonality in the data
train_sd <- train_bc %>%
  mutate(temperature = difference(temperature,lag=12))%>%
  drop_na(temperature)
train_sd%>%ACF(lag_max=24)%>%autoplot
train_sd%>%PACF(lag_max=24)%>%autoplot
# both ACF & PACF indicate overdifferencing on seasonal level --> add seasonal MA term
# first lag slightly significant in ACF --> either p=0 or p=1

models_2 <- train_bc %>%
  model(
    arima_guess1 = ARIMA(temperature ~ 1 + pdq(0,0,0)+PDQ(0,1,1)),
    arima_guess2 = ARIMA(temperature ~ 1 + pdq(1,0,0)+PDQ(0,1,1)),
    arima_auto = ARIMA(temperature)
  )

report(models_2%>%select(arima_guess1))
report(models_2%>%select(arima_guess2))
report(models_2%>%select(arima_auto))
# all IC slightly better for guessed model --> go with guessed model

# coefficient check to choose between guess 1 and 2
coef(models_2%>%select(arima_guess1))
coef(models_2%>%select(arima_guess2))
# ar1 component in guess 2 not significant --> difference in BIC only small, therefore go with guess1

# Residual Evaluation
models_2 %>%
  select(arima_guess1) %>%
  gg_tsresiduals(type = "innovation")

# Ljung Box Test (we want HIGH p-value -> NULL not rejected, then resid. independent)
models_2 %>%
  select(arima_guess1) %>%
  residuals()  %>%
  features(.resid, features = ljung_box, lag = 24, dof = 8)

####### Forecasting #######

# train final models on untransformed data to produce forecasts
train <- df %>%
  filter_index(. ~ "2020 Jan")
test <- df %>%
  filter_index("2020 Jan" ~ .)

# set up best models
models_3 <- train %>%
  model(
    arima_guess1 = ARIMA(temperature ~ 1 + pdq(0,0,0)+PDQ(0,1,1)),
    ets_guess = ETS(temperature ~ error("A") + trend("N") + season("A"))
  )

# Forecast of all models for test set (with train set)
models_3 %>%
  forecast(h=24)%>%
  autoplot(test,level = NULL)

# Forecast guessed ETS on test set
models_3 %>%
  select(ets_guess) %>%
  forecast(h=24) %>%
  autoplot(test) +
  ylim(-2,22)+
  labs(title = "Forecast from ETS model") + labs(y= "Temperature(C°)", x="Month")

# Forecast guessed ARIMA on test set
models_3 %>%
  select(arima_guess1) %>%
  forecast(h=24) %>%
  autoplot(test) +
  ylim(-2,22)+
  labs(title = "Forecast from ARIMA model") + labs(y= "Temperature(C°)", x="Month")

models_3 %>%
  forecast(test) %>%
  accuracy(df)

models_3 %>%
  forecast(train) %>%
  accuracy(df)

# final forecast with ARIMA
models_3 %>%
  select(ets_guess) %>%
  forecast(h=24) %>%
  autoplot(df) +
  labs(title = "Overall forecast from guessed ARIMA model")

# final model for 2022-2023 forecast
models_final <- df%>%
  model(ets_guess = ETS(temperature ~ error("A") + trend("N") + season("A")))

models_final %>%
  forecast(h=24)%>%
  autoplot(df) + labs(y= "Temperature(C°)", x="Month")



