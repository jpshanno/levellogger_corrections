source("code/levellogger_packages.R"); source("code/levellogger_functions.R")
loadd(mesonet_response, raw_barometric_data, levellogger_measurements)

f3311 <- mesonet_response[station == "F3311"]
f3311[, ex_pressure_error_cm := 2.7625 / qnorm(0.99)]

raw <- 
  read_xle_data(levellogger_measurements, "air") %>% 
  .[!between(sample_time, 
             as.POSIXct("2020-06-10 06:00", tz = "EST5EDT"), 
             as.POSIXct("2020-06-10 12:30", tz = "EST5EDT"))]

# Issue with large spikes in temp. Removed those points for now. For real data
# use 12-hour temperature

# Worked with dens sea level pressure -> temp -> delta_t_c

# Thinking about using dens_air_pressure ~ ex_abs_pressure. Want to correct the
# raw pressure and then use SLP for air and water compensation

# Comparing temperature compensation of var-dis and var-sim developed with 
# training data and tested on test-dat. Altimeter performed best
# Sea Level Pressure:
# serial_number         V1
# 1:       1066019 0.01234704
# 2:       1065861 0.01161183

# Altimeter:
# serial_number          V1
# 1:       1066019 0.005650081
# 2:       1065861 0.004877106

# Absolute Pressure:
# serial_number         V1
# 1:       1066019 0.01328713
# 2:       1065861 0.01252859

baro_dat <-
  raw %>%
  .[, dens_air_pressure_cm := 1000 * air_pressure_cm / water_density(air_temperature_c)] %>% 
  .[f3311, on = "sample_time", nomatch = NULL] %>% 
  .[, `:=`(temperature_difference_c = ex_air_temperature_c - air_temperature_c,
           raw_logger_error_cm = dens_air_pressure_cm - ex_abs_pressure_cm)] %>% 
  set_experiments("data/experimental_design.csv")

baro_training <- 
  baro_dat[experiment != "test-dat"]

baro_training[, air_temp2 := air_temperature_c]

# Center errors around linear relationship

vardis_means <- 
  baro_training[, 
                calculate_alignments(data = .SD, 
                                     full.data = baro_training, 
                                     reference.experiment = "var-dis",
                                     x = "air_temperature_c",
                                     y = "raw_logger_error_cm"),
                by = .(experiment)]

baro_training[, experimental_mean := mean(raw_logger_error_cm),
              by = .(serial_number, experiment)]

baro_training[, alignment_offset_cm := NULL]

baro_training[vardis_means, 
              alignment_offset_cm := alignment_offset_cm, 
              on = c("serial_number", "experiment")]

baro_training[, experimental_error := experimental_mean - alignment_offset_cm]

baro_training[, logger_error_cm := raw_logger_error_cm - experimental_error]

dat <- 
  baro_training[serial_number == "1065861"]

dat[, error_y_cm := combine_errors(pressure_error_cm, ex_pressure_error_cm)]

# library(brms)
# 
# # fit a simple error-in-variables model 
# fit1 <- brm(logger_error_cm ~ me(air_temperature_c, temperature_error_c), 
#             data = dat, 
#             save_mevars = TRUE)
# summary(fit1)
# 
# ggplot(data = cbind(dat[, .(air_temperature_c)], 
#                     fitted(fit1)), 
#        aes(x = air_temperature_c, y = Estimate)) +
#   geom_linerange(aes(ymin = Q2.5, ymax = Q97.5)) + 
#   geom_line()

library(rstan)

fit2 <- 
  stan("tmp/test.stan", 
       data = list(N = nrow(dat), 
                   tau = dat$temperature_error_c, 
                   x_meas = dat$air_temperature_c, 
                   y = dat$logger_error_cm,
                   tau_y = dat$error_y_cm), 
       # pars = c("alpha", "beta", "sigma"),
       # iter = 2500
       )

print(fit2, pars = c("alpha", "beta"))

# From:
# https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed

alpha_post <- extract(fit2)$alpha
beta_post <- extract(fit2)$beta

# Function for simulating y based on new x
predict_post <- function(x, alpha = 0.05) {
  set.seed(1234)
  alpha <- matrix(rep(sample(alpha_post, 1000, replace = TRUE),
                      times = length(x)), 
                  nrow = length(x))
  beta <- matrix(rep(sample(beta_post, 1000, replace = TRUE),
                     times = length(x)), 
                 nrow = length(x))
  mat <- 
    alpha + beta * x
  
  cbind(estimate = apply(mat, 1, median),
        lwr = apply(mat, 1, quantile, probs = 0.025),
        upr = apply(mat, 1, quantile, probs = 0.975))
}

trad_model <- 
  lm(logger_error_cm ~ air_temperature_c, data = dat)
  

# Run the function on x_test
plot(x = dat$air_temperature_c,
     y = dat$logger_error_cm, )
points(x = dat$air_temperature_c,
       y = predict_post(dat$air_temperature_c)[, "estimate"], 
       col = "blue",
       pch = 20)
segments(x0 = dat$air_temperature_c,
         y0 = predict_post(dat$air_temperature_c)[, "lwr"],
         y1 = predict_post(dat$air_temperature_c)[, "upr"],
         col = "red")
points(x = dat$air_temperature_c,
       y = predict(trad_model), 
       type = "l")
points(x = dat$air_temperature_c,
       y = predict(trad_model, interval = "prediction")[,"lwr"],
       type = "l", 
       col = "green")
points(x = dat$air_temperature_c,
       y = predict(trad_model, interval = "prediction")[,"upr"],
       type = "l", 
       col = "green")

