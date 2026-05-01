library(quantmod)
library(zoo)
library(dplyr)
library(forecast)
library(tseries)
library(lubridate)

#Data acquistion
getSymbols(c("DGS10", "DCOILBRENTEU"), src = "FRED",
           from = "2000-01-01", to = "2024-12-31")

dgs10_df <- data.frame(
  date  = floor_date(as.Date(index(apply.monthly(DGS10, mean, na.rm = TRUE))), "month"),
  yield = as.numeric(apply.monthly(DGS10, mean, na.rm = TRUE))
)

brent_df <- data.frame(
  date  = floor_date(as.Date(index(apply.monthly(DCOILBRENTEU, mean, na.rm = TRUE))), "month"),
  brent = as.numeric(apply.monthly(DCOILBRENTEU, mean, na.rm = TRUE))
)

combined_df <- inner_join(dgs10_df, brent_df, by = "date") %>% arrange(date)

#Data transformation
yield_ts <- ts(combined_df$yield, start = c(2000, 1), frequency = 12)
brent_ts <- ts(combined_df$brent, start = c(2000, 1), frequency = 12)

diff_yield <- diff(log(yield_ts))
diff_brent <- diff(log(brent_ts))

diff_df <- data.frame(
  date       = combined_df$date[-1],
  diff_yield = as.numeric(diff_yield),
  diff_brent = as.numeric(diff_brent)
)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(yield_ts, col = "steelblue", lwd = 1.5,
     main = "Figure 1a: U.S. 10-Year Treasury Yield (2000-2024)",
     ylab = "Yield (%)", xlab = "")
plot(brent_ts, col = "darkorange", lwd = 1.5,
     main = "Figure 1b: Brent Crude Oil Price (2000-2024)",
     ylab = "USD per Barrel", xlab = "Date")

#Stationarity tests
adf.test(log(yield_ts))
adf.test(log(brent_ts))
adf.test(diff_yield)
adf.test(diff_brent)

#ACF/PACF of stationary series
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
acf(diff_yield,  lag.max = 36, main = "Figure 2a: ACF - Delta Log Yield")
pacf(diff_yield, lag.max = 36, main = "Figure 2b: PACF - Delta Log Yield")
acf(diff_brent,  lag.max = 36, main = "Figure 2c: ACF - Delta Log Brent")
pacf(diff_brent, lag.max = 36, main = "Figure 2d: PACF - Delta Log Brent")

#ARIMA eval
m010 <- Arima(diff_yield, order = c(0, 0, 0))
m110 <- Arima(diff_yield, order = c(1, 0, 0))
m011 <- Arima(diff_yield, order = c(0, 0, 1))
m111 <- Arima(diff_yield, order = c(1, 0, 1))
m210 <- Arima(diff_yield, order = c(2, 0, 0))
m012 <- Arima(diff_yield, order = c(0, 0, 2))

#AIC/BIC comparison
models      <- list(m010, m110, m011, m111, m210, m012)
model_names <- c("ARIMA(0,1,0)", "ARIMA(1,1,0)", "ARIMA(0,1,1)",
                 "ARIMA(1,1,1)", "ARIMA(2,1,0)", "ARIMA(0,1,2)")
data.frame(Model = model_names,
           AIC   = sapply(models, AIC),
           BIC   = sapply(models, BIC))

best_arima <- m011

#ARIMA Residual eval
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(residuals(best_arima), col = "steelblue", lwd = 1,
     main = "Figure 3a: Residuals of Best ARIMA Model",
     ylab = "Residual", xlab = "")
abline(h = 0, lty = 2, col = "gray")
acf(residuals(best_arima), lag.max = 36,
    main = "Figure 3b: ACF of ARIMA(0,1,1) Residuals")

Box.test(residuals(best_arima), lag = 20, type = "Ljung-Box")

#ARIMAX Lag eval
n          <- length(diff_yield)
brent_lag0 <- as.numeric(diff_brent)
brent_lag1 <- c(NA, as.numeric(diff_brent[-n]))
brent_lag2 <- c(NA, NA, as.numeric(diff_brent[-(n-0:1)]))
brent_lag3 <- c(NA, NA, NA, as.numeric(diff_brent[-(n-0:1)-1]))

fit_arimax <- function(xreg_raw) {
  valid <- !is.na(xreg_raw)
  Arima(diff_yield[valid], order = c(0, 0, 1), xreg = xreg_raw[valid])
}

arimax_lag0 <- fit_arimax(brent_lag0)
arimax_lag1 <- fit_arimax(brent_lag1)
arimax_lag2 <- fit_arimax(brent_lag2)
arimax_lag3 <- fit_arimax(brent_lag3)

#Lag comparison table
arimax_models  <- list(best_arima, arimax_lag0, arimax_lag1, arimax_lag2, arimax_lag3)
arimax_labels  <- c("ARIMA(0,1,1) Baseline", "ARIMAX (lag 0)",
                    "ARIMAX (lag 1)", "ARIMAX (lag 2)", "ARIMAX (lag 3)")
data.frame(Model = arimax_labels,
           AIC   = sapply(arimax_models, AIC),
           BIC   = sapply(arimax_models, BIC))

best_arimax <- arimax_lag0

#Param estimates
coef_vals   <- coef(best_arimax)
coef_se     <- sqrt(diag(vcov(best_arimax)))
coef_tstat  <- coef_vals / coef_se
coef_pval   <- 2 * pnorm(-abs(coef_tstat))
data.frame(Estimate = coef_vals, SE = coef_se,
           t = coef_tstat, p = coef_pval)

Box.test(residuals(best_arimax), lag = 20, type = "Ljung-Box")

train_end  <- which(diff_df$date == "2022-12-01")
test_start <- train_end + 1

train_yield <- diff_yield[1:train_end]
train_brent <- diff_brent[1:train_end]
test_yield  <- as.numeric(diff_yield[test_start:n])
test_brent  <- as.numeric(diff_brent[test_start:n])
h           <- length(test_yield)

m_naive  <- Arima(train_yield, order = c(0, 0, 0))
m_arima  <- Arima(train_yield, order = c(0, 0, 1))
m_arimax <- Arima(train_yield, order = c(0, 0, 1), xreg = train_brent)
m_hw     <- HoltWinters(train_yield, beta = FALSE, gamma = FALSE)

fc_naive  <- forecast(m_naive,  h = h)
fc_arima  <- forecast(m_arima,  h = h)
fc_arimax <- forecast(m_arimax, h = h, xreg = matrix(test_brent, ncol = 1))
fc_hw     <- forecast(m_hw,     h = h)

#RMSE and MAE
accuracy_table <- function(fc, actual, label) {
  pred <- as.numeric(fc$mean)
  data.frame(Model = label,
             RMSE  = sqrt(mean((pred - actual)^2)),
             MAE   = mean(abs(pred - actual)))
}

rbind(
  accuracy_table(fc_naive,  test_yield, "Random Walk"),
  accuracy_table(fc_arima,  test_yield, "ARIMA(0,1,1)"),
  accuracy_table(fc_arimax, test_yield, "ARIMAX (lag 0)"),
  accuracy_table(fc_hw,     test_yield, "Holt-Winters")
)

test_dates <- diff_df$date[test_start:nrow(diff_df)]

par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
plot(test_dates, test_yield, type = "l", col = "black", lwd = 2,
     main = "Figure 4: Model Forecasts vs Actual Delta Log Yield (2023-2024)",
     ylab = "Delta Log Yield", xlab = "Date",
     ylim = range(c(test_yield, fc_naive$mean, fc_arima$mean,
                    fc_arimax$mean, fc_hw$mean)))
lines(test_dates, fc_naive$mean,  col = "gray",      lwd = 1.5, lty = 2)
lines(test_dates, fc_arima$mean,  col = "steelblue", lwd = 1.5, lty = 2)
lines(test_dates, fc_arimax$mean, col = "darkorange", lwd = 1.5, lty = 2)
lines(test_dates, fc_hw$mean,     col = "darkgreen", lwd = 1.5, lty = 2)
legend("topright",
       legend = c("Actual", "Random Walk", "ARIMA(0,1,1)",
                  "ARIMAX (lag 0)", "Holt-Winters"),
       col    = c("black", "gray", "steelblue", "darkorange", "darkgreen"),
       lty    = c(1, 2, 2, 2, 2), lwd = c(2, 1.5, 1.5, 1.5, 1.5), bty = "n")

#12-month forward forecast

m_arimax_full  <- Arima(diff_yield, order = c(0, 0, 1), xreg = diff_brent)
brent_future   <- rep(mean(tail(as.numeric(diff_brent), 24)), 12)
fc_forward     <- forecast(m_arimax_full, h = 12,
                           xreg = matrix(brent_future, ncol = 1))
print(fc_forward)