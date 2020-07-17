as.numeric(st_distance(st_sfc(st_point(c(-88.58083, 47.092)), st_point(c(-88.5742812, 47.119242)), crs = "+proj=longlat +datum=WGS84"))[1,2]/1000)`

source("code/levellogger_packages.R"); source("code/levellogger_functions.R")
loadd(raw_water_data, levellogger_measurements, mesonet_response, raw_barometric_data)

f3311 <- mesonet_response[station == "F3311"]
kcmx <- mesonet_response[station == "KCMX"]


# Water Correction --------------------------------------------------------


levellogger_measurements[, `:=`(water_depth_cm = cable_length_cm - water_height_cm)]

water_dat <- 
  raw_water_data[f3311, on = "sample_time", nomatch = NULL] %>% 
  .[, delta_t_c_min := c(NA_real_, diff(water_temperature_c)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
    by = .(serial_number)] %>% 
  set_experiments("data/experimental_design.csv") %>% 
  .[, `:=`(temperature_difference_c = ex_air_temperature_c - water_temperature_c, 
           abs_temperature_difference_c = abs(ex_air_temperature_c - water_temperature_c), 
           water_level_cm = water_pressure_cm - ex_sea_level_pressure_cm + 22500/826, 
           dens_water_level_cm = dens_water_pressure_cm - ex_sea_level_pressure_cm + 22500/826)]

water_dat[levellogger_measurements, 
         water_depth_cm := water_depth_cm, 
         on = "serial_number"]

# There are slight temperature difference errors when water is not density 
# adjusted

water_dat[, raw_logger_error_cm := dens_water_pressure_cm - (water_depth_cm + ex_sea_level_pressure_cm) + 22500/826]

water_dat[, error_range_cm := error_range(pressure_error_cm, rep(0.5, nrow(.SD)))]
water_dat[, elapsed_s := as.numeric(sample_time - min(sample_time)), 
         by = .(serial_number)]

water_training <- 
  water_dat[experiment != "test-dat"]

water_testing <- 
  water_dat[experiment == "test-dat"]

# Show plot of error using absolute pressure vs error using SLP
ggplot(water_training,
       aes(color = experiment, 
           y = raw_logger_error_cm,
           x = temperature_difference_c)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")

# Center errors around linear relationship

water_training[, temp_diff2 := temperature_difference_c]

statsim_means <- 
  water_training[, 
                .(statsim_mean = water_training[experiment == "stat-sim"][
                  .SD[,.(min_temp = 0.95 * min(temperature_difference_c), 
                         max_temp = 1.05 * max(temperature_difference_c))],
                  on = c("temp_diff2>=min_temp", "temp_diff2<=max_temp")][, 
                                                                        mean(raw_logger_error_cm)]),
                by = .(serial_number, experiment)]

water_training[, experimental_mean := mean(raw_logger_error_cm),
              by = .(serial_number, experiment)]

water_training[statsim_means, 
               statsim_mean := statsim_mean, 
               on = c("serial_number", "experiment")]

water_training[, experimental_error := experimental_mean - statsim_mean]

ggplot(water_training,
       aes(x = temperature_difference_c, 
           color = experiment)) + 
  geom_point(aes(y = raw_logger_error_cm),
             alpha = 0.25) +
  geom_point(aes(y = raw_logger_error_cm - experimental_error)) +
  facet_wrap(~serial_number, scales = "free")

water_training[, logger_error_cm := raw_logger_error_cm - experimental_error]


water_tdiff_models <-
  fit_correction_mod(data = water_training,
                     x = "temperature_difference_c",
                     y = "logger_error_cm",
                     by = "serial_number")

# Could justify one tdiff model here using #   fit_correction_mmod(data = water_training, x = "temperature_difference_c", y = "logger_error_cm", group = "serial_number", control = lmerControl(optimizer = "Nelder_Mead"))
# fit_correction_mod(data = water_training,
#                    x = "temperature_difference_c",
#                    y = "logger_error_cm",
#                    by = "serial_number") %>% 
#   {ggplot(., aes(x = slope)) + 
#       geom_density(trim = FALSE) + 
#       geom_density(data = data.frame(slope = rnorm(1000, 
#                                                    mean(.$slope), 
#                                                    sd(.$slope))), 
#                    color = "red")}

water_training <- 
  correction_residuals(water_tdiff_models,
                       original.data = water_training,
                       x = "temperature_difference_c",
                       y = "logger_error_cm",
                       by = "serial_number")

ggplot(water_training,
       aes(y = residuals,
           x = temperature_difference_c)) + 
  geom_point() + 
  facet_wrap(~serial_number)

water_training[water_tdiff_models, 
               tdiff_water_pressure_cm := dens_water_pressure_cm - slope * (temperature_difference_c - x_intercept),
               on = "serial_number"]

water_training[, 
               tdiff_logger_error_cm := tdiff_water_pressure_cm - (water_depth_cm + ex_sea_level_pressure_cm) + 22500/826]

ggplot(water_training,
       aes(y = tdiff_logger_error_cm,
           x = temperature_difference_c,
           color = water_temperature_c)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")

# Need to look at 2030899, 2064734 - may be something about temperature error?
# 2030899 had worst clock error
ggplot(water_training,
       aes(y = tdiff_logger_error_cm,
           x = water_temperature_c,
           color = temperature_difference_c)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")

# Temp Models
water_temperature_models <- 
  fit_correction_mod(water_training,
                     x = "water_temperature_c",
                     y = "tdiff_logger_error_cm",
                     by = "serial_number")

mcp_mod <- 
  mcp(list(tdiff_logger_error_cm ~ water_temperature_c, 
           ~ water_temperature_c), 
      data = water_training[serial_number == "2030899"])

plot(mcp_mod, q_predict = TRUE)

water_training <- 
  correction_residuals(water_temperature_models,
                       water_training,
                       x = "water_temperature_c",
                       y = "tdiff_logger_error_cm",
                       by = "serial_number")

ggplot(water_training,
       aes(y = residuals,
           x = water_temperature_c)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")

ggplot(water_training,
       aes(y = residuals,
           x = temperature_difference_c)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")

water_training[water_temperature_models, 
         temp_tdiff_water_pressure_cm := tdiff_water_pressure_cm - slope * (water_temperature_c - x_intercept),
         on = "serial_number"]

water_training[, 
         temp_tdiff_logger_error_cm := temp_tdiff_water_pressure_cm - (water_depth_cm + ex_sea_level_pressure_cm) + 22500/826]

# Need to add in sensor resolution
# Differs for Edge and Junior
# No loggers with average change > 0.05 C/min
# 
ggplot(water_training[abs(delta_t_c_min) > 0.003/2],
       aes(y = temp_tdiff_logger_error_cm,
           x = delta_t_c_min)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")


water_deltat_models <- 
  fit_correction_mod(water_training[abs(delta_t_c_min) > 0.003/2],
                     x = "delta_t_c_min",
                     y = "temp_tdiff_logger_error_cm",
                     by = "serial_number")

water_training <- 
  correction_residuals(water_deltat_models,
                       water_training,
                       x = "delta_t_c_min",
                       y = "temp_tdiff_logger_error_cm",
                       by = "serial_number")

water_training[water_deltat_models, 
               deltat_temp_tdiff_water_pressure_cm := temp_tdiff_water_pressure_cm - slope * (delta_t_c_min - x_intercept),
               on = "serial_number"]

water_training[, 
               deltat_temp_logger_error_cm := deltat_temp_tdiff_water_pressure_cm - (water_depth_cm + ex_sea_level_pressure_cm) + 22500/826]



# Final Offset
water_training[, `:=`(raw_offset_cm = mean(water_pressure_cm - ex_sea_level_pressure_cm - water_depth_cm),
                      rect_offset_cm = mean(deltat_temp_tdiff_water_pressure_cm - ex_sea_level_pressure_cm - water_depth_cm)), 
               by = .(serial_number)]
water_training[, `:=`(water_level_cm = water_pressure_cm - ex_sea_level_pressure_cm - raw_offset_cm,
                      rect_water_level_cm = deltat_temp_tdiff_water_pressure_cm - ex_sea_level_pressure_cm - rect_offset_cm)]

ggplot(water_training,
       aes(x = sample_time,
           y = rect_water_level_cm,
           color = experiment)) +
  geom_ribbon(aes(ymin = water_depth_cm - error_range_cm / 2,
                  ymax = water_depth_cm + error_range_cm / 2),
              fill = "gray80") +
  geom_line() +
  geom_line(aes(y = water_level_cm),
            linetype = "dashed") + 
  facet_wrap(~serial_number)

# Test-dat

water_testing[water_tdiff_models, 
              rect_water_pressure_cm := dens_water_pressure_cm - slope * (temperature_difference_c - x_intercept),
              on = "serial_number"]

water_testing[water_temperature_models, 
              rect_water_pressure_cm := rect_water_pressure_cm - slope * (water_temperature_c - x_intercept),
              on = "serial_number"]

water_testing[water_deltat_models, 
              full_rect_water_pressure_cm := rect_water_pressure_cm - slope * (delta_t_c_min - x_intercept),
              on = "serial_number"]

water_testing[, `:=`(raw_offset_cm = mean(water_pressure_cm - ex_sea_level_pressure_cm - water_depth_cm),
                     rect_offset_cm = mean(rect_water_pressure_cm - ex_sea_level_pressure_cm - water_depth_cm),
                     full_rect_offset_cm = mean(full_rect_water_pressure_cm - ex_sea_level_pressure_cm - water_depth_cm)), 
              by = .(serial_number)]
water_testing[, `:=`(water_level_cm = water_pressure_cm - ex_sea_level_pressure_cm - raw_offset_cm,
                     rect_water_level_cm = rect_water_pressure_cm - ex_sea_level_pressure_cm - rect_offset_cm,
                     full_rect_water_level_cm = full_rect_water_pressure_cm - ex_sea_level_pressure_cm - full_rect_offset_cm)]

ggplot(water_testing,
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = water_depth_cm - error_range_cm / 2,
                  ymax = water_depth_cm + error_range_cm / 2),
              fill = "gray80") +
  geom_line(aes(y = full_rect_water_level_cm)) +
  geom_line(aes(y = rect_water_level_cm),
            linetype = "dashed") + 
  geom_line(aes(y = water_level_cm),
            linetype = "dotted") + 
  facet_wrap(~serial_number)


# Baro Correction ---------------------------------------------------------
Logger error = dens_air_pressure_cm - ex_abs_pressure_cm
Tdiff model with only var-sim, var-dis and abs_temp_diff > 5
Temperature model with only var-sim & var-dis


# Calculate sea level pressure
# Apply density correction
# Correct for air temperature effect using models from var experiments
# Correct for 

# This method works for F3311, KCMX, and both datasets combined. The only
# difference comes in the slope of the temperature-compensation models

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
  .[, abs_pressure_hPa := 0.01 * air_pressure_cm * 9.80665 * 1000] %>%
  .[, dens_abs_pressure_hPa := 0.01 * dens_air_pressure_cm * 9.80665 * 1000] %>%
  .[, dens_altimeter_hPa := (dens_abs_pressure_hPa - 0.3) * (1 + ((1013.25^0.190284*0.0065)/288)*(225/((dens_abs_pressure_hPa-0.3)^0.190284)))^(1/0.190284)] %>% 
  .[, dens_altimeter_cm := 100 * dens_altimeter_hPa / (9.80665 * 1000)] %>% 
  .[, dens_sea_level_pressure_hPa := dens_abs_pressure_hPa * (1 - (0.0065 * 225) / (air_temperature_c + 273.15 + 0.0065 * 225))^-5.257] %>%
  .[, dens_sea_level_pressure_cm := 100 * dens_sea_level_pressure_hPa / (9.80665 * 1000)] %>%
  .[, `:=`(temperature_difference_c = ex_air_temperature_c - air_temperature_c,
           abs_temperature_difference_c = abs(ex_air_temperature_c - air_temperature_c),
           raw_logger_error_cm = dens_air_pressure_cm - ex_abs_pressure_cm - ((elevation_cm - 22500) / 826))] %>% 
  .[, `:=`(delta_t_c_min = c(NA_real_, diff(air_temperature_c)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
           delta_p_cm_min = c(NA_real_, diff(air_pressure_cm)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins")))),
    by = .(serial_number)] %>% 
  set_experiments("data/experimental_design.csv")

baro_dat[, error_range_cm := error_range(pressure_error_cm, rep(0, nrow(.SD)))]

baro_training <- 
  baro_dat[experiment != "test-dat"]

baro_testing <- 
  baro_dat[experiment == "test-dat"]

baro_training[, air_temp2 := air_temperature_c]

# Center errors around linear relationship

vardis_means <- 
  baro_training[, 
                .(vardis_mean = baro_training[experiment == "var-dis"][
                .SD[,.(min_temp = 0.95 * min(air_temperature_c), 
                       max_temp = 1.05 * max(air_temperature_c))],
                on = c("air_temp2>=min_temp", "air_temp2<=max_temp")][, 
                                                                      mean(raw_logger_error_cm)]),
                by = .(serial_number, experiment)]

baro_training[, experimental_mean := mean(raw_logger_error_cm),
              by = .(serial_number, experiment)]

baro_training[vardis_means, 
              vardis_mean := vardis_mean, 
              on = c("serial_number", "experiment")]

baro_training[, experimental_error := experimental_mean - vardis_mean]

ggplot(baro_training,
       aes(x = air_temperature_c, 
           color = experiment)) + 
  geom_point(aes(y = raw_logger_error_cm), 
             alpha = 0.25) +
  geom_point(aes(y = raw_logger_error_cm - experimental_error)) +
  facet_wrap(~serial_number, scales = "free")

baro_training[, logger_error_cm := raw_logger_error_cm - experimental_error]


# Temp Models
baro_temperature_models <- 
  fit_correction_mod(baro_training,
                     x = "air_temperature_c",
                     y = "logger_error_cm",
                     by = "serial_number")

baro_training <- 
  correction_residuals(baro_temperature_models,
                       baro_training,
                       x = "air_temperature_c",
                       y = "logger_error_cm",
                       by = "serial_number")

ggplot(baro_training,
       aes(y = residuals,
           x = air_temperature_c,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~serial_number)

ggplot(baro_training,
       aes(y = residuals,
           x = delta_t_c_min,
           color = experiment)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~serial_number)

baro_training[baro_temperature_models, 
              temp_dens_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
               on = "serial_number"]

baro_training[, 
               temp_dens_logger_error_cm := temp_dens_air_pressure_cm - ex_abs_pressure_cm - (30236 - 22500)/826]


ggplot(baro_training,
       aes(x = air_temperature_c, 
           y = logger_error_cm)) + 
  geom_point(aes(color = experiment)) +
  geom_line(aes(y = logger_error_cm - residuals)) +
  # geom_smooth(method = "lm") +
  facet_wrap(~serial_number, scales = "free")


ggplot(baro_training[abs(delta_t_c_min) > 0.1], 
       aes(x = delta_t_c_min, 
           y = temp_dens_logger_error_cm)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~serial_number)

# Delta_t mods
# Looking at the data, a delta_t_c_min of 0.2 highlights the problem areas. 
# Doing that centered on zero (0.1)
baro_delta_t_models <- 
  fit_correction_mod(baro_training[abs(delta_t_c_min) > 0.05],
                     x = "delta_t_c_min",
                     y = "temp_dens_logger_error_cm",
                     by = "serial_number")

baro_training <- 
  correction_residuals(baro_temperature_models,
                       baro_training,
                       x = "delta_t_c_min",
                       y = "temp_dens_logger_error_cm",
                       by = "serial_number")


baro_training[baro_delta_t_models,
              delta_t_adjustment := slope * (delta_t_c_min) - x_intercept,
              on = "serial_number"]

baro_training[baro_delta_t_models, 
              deltat_temp_dens_air_pressure_cm := temp_dens_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
              on = "serial_number"]

baro_training[baro_delta_t_models,
              deltat_temp_dens_air_pressure_cm := deltat_temp_dens_air_pressure_cm - ifelse(abs(delta_t_c_min) > 0.05, intercept, 0),
              on = "serial_number"]

baro_training[, 
              deltat_temp_dens_logger_error_cm := deltat_temp_dens_air_pressure_cm - ex_abs_pressure_cm - (30236 - 22500)/826]

ggplot(baro_training,
       aes(y = deltat_temp_dens_logger_error_cm,
           x = air_temperature_c,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~serial_number)

ggplot(baro_training,
       aes(y = deltat_temp_dens_logger_error_cm,
           x = delta_t_c_min,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~serial_number)

# Final Offset
# Unique by experiment because of shifts in logger hanging position
baro_training[, `:=`(raw_offset_cm = mean(air_pressure_cm - ex_abs_pressure_cm),
                     rect_offset_cm = mean(deltat_temp_dens_air_pressure_cm - ex_abs_pressure_cm)), 
              by = .(serial_number, experiment)]
baro_training[, `:=`(baro_level_cm = air_pressure_cm - raw_offset_cm,
                     rect_baro_level_cm = deltat_temp_dens_air_pressure_cm - rect_offset_cm)]

baro_offsets <- 
  baro_training[, lapply(.SD, mean),
                .SDcols = c("raw_offset_cm", "rect_offset_cm"),
                by = c("serial_number", "experiment")]

ggplot(baro_training,
       aes(x = sample_time,
           y = rect_baro_level_cm,
           color = experiment,
           fill = experiment)) +
  geom_ribbon(aes(ymin = ex_abs_pressure_cm - error_range_cm / 2,
                  ymax = ex_abs_pressure_cm + error_range_cm / 2),
              color = "gray80",
              alpha = 0.3) +
  geom_line() +
  geom_line(aes(y = baro_level_cm),
            linetype = "dotted",
            color = "gray40") + 
  facet_wrap(~ serial_number + experiment,
             scales = "free",
             nrow = 2)

# Test-dat

# baro_testing[baro_tdiff_models, 
#              rect_air_pressure_cm := dens_air_pressure_cm - slope * (temperature_difference_c - x_intercept),
#              on = "serial_number"]

baro_testing[baro_temperature_models, 
              rect_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
              on = "serial_number"]

baro_testing[baro_delta_t_models,
             full_rect_air_pressure_cm := rect_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
              on = "serial_number"]


# baro_testing[baro_delta_t_models,
#              full_rect_air_pressure_cm := rect_air_pressure_cm - ifelse(abs(delta_t_c_min) > 0.05, slope * delta_t_c_min, 0),
#              on = "serial_number"]

# Not using standard offsets from baro_training because it changes slightly each
# time the loggers are hung
baro_testing[, `:=`(raw_offset_cm = mean(air_pressure_cm - ex_abs_pressure_cm),
                    rect_offset_cm = mean(rect_air_pressure_cm - ex_abs_pressure_cm),
                    full_rect_offset_cm = mean(full_rect_air_pressure_cm - ex_abs_pressure_cm)), 
             by = .(serial_number)]

baro_testing[,`:=`(raw_baro_level_cm = air_pressure_cm - raw_offset_cm,
                   rect_baro_level_cm = rect_air_pressure_cm - rect_offset_cm,
                   full_rect_baro_level_cm = full_rect_air_pressure_cm - full_rect_offset_cm)]

baro_testing[, mean((rect_air_pressure_cm - ex_abs_pressure_cm) / ex_abs_pressure_cm), by = "serial_number"]

# Sea Level Pressure:
# serial_number         V1
# 1:       1066019 0.01230583
# 2:       1065861 0.0114791

# serial_number         V1
# 1:       1066019 0.01249084
# 2:       1065861 0.01174087

# Absolute Pressure:
# serial_number         V1
# 1:       1066019 0.01043319
# 2:       1065861 0.01009805

# Absolute Pressure w/ tdiff:
# serial_number          V1
# 1:       1066019 0.011926889
# 2:       1065861 0.008829354

baro_testing[, delta_t_adjustment := NULL]

baro_testing[baro_delta_t_models,
             delta_t_adjustment := slope * (delta_t_c_min - x_intercept) - intercept,
             on = "serial_number"]

baro_testing[abs(delta_t_c_min) <= 0.05, delta_t_adjustment := 0]

baro_testing[,cum_delta_t_adjustment := cumsum(delta_t_adjustment),
             by = "serial_number"]

ggplot(baro_testing, aes(x = sample_time,
                         y = cum_delta_t_adjustment)) +
  geom_line() + 
  facet_wrap(~serial_number, scales = "free")

{ggplot(baro_testing,
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = ex_abs_pressure_cm - error_range_cm / 2,
                  ymax = ex_abs_pressure_cm + error_range_cm / 2),
              fill = "gray80") + 
  geom_line(aes(y = ex_abs_pressure_cm),
            color = "gray70") +
  geom_line(aes(y = (full_rect_baro_level_cm + cum_delta_t_adjustment) - mean(full_rect_baro_level_cm + cum_delta_t_adjustment - ex_abs_pressure_cm)),
            linetype = "solid") +
  geom_line(aes(y = full_rect_baro_level_cm),
            linetype = "dashed") +
    geom_line(aes(y = raw_baro_level_cm),
              linetype = "dotted") +
  facet_wrap(~serial_number) +
  guides(linetype = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())} /
{ggplot(baro_testing,
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = ex_sea_level_pressure_cm - error_range_cm / 2,
                  ymax = ex_sea_level_pressure_cm + error_range_cm / 2),
              fill = "gray80") + 
  geom_line(aes(y = ex_sea_level_pressure_cm),
            color = "gray70") +
    geom_line(aes(y = slp(full_rect_baro_level_cm - cum_delta_t_adjustment, air_temperature_c) - mean(slp(full_rect_baro_level_cm - cum_delta_t_adjustment, air_temperature_c) - ex_sea_level_pressure_cm),
              linetype = "Temperature & Detla T & DetlaT Offset")) +
  geom_line(aes(y = slp((full_rect_baro_level_cm) , air_temperature_c) - mean(slp((full_rect_baro_level_cm) , air_temperature_c) - ex_sea_level_pressure_cm),
            linetype = "Temp & Delta T")) +
    geom_line(aes(y = slp(raw_baro_level_cm, air_temperature_c) - mean(slp(raw_baro_level_cm, air_temperature_c) - ex_sea_level_pressure_cm),
              linetype = "Raw")) +
  facet_wrap(~serial_number) +
  scale_linetype_manual(name = NULL,
                        values = c("dotted", "dashed", "solid"),
                        breaks = c("Raw", "Temp & Delta T", "Temperature & Detla T & DetlaT Offset")) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(legend.position = "bottom")}

baro_testing[, lapply(.SD, function(x){sum(!between(x, ex_abs_pressure_cm - error_range_cm / 2, ex_abs_pressure_cm + error_range_cm / 2))}),
          by = .(serial_number),
          .SDcols = c("raw_baro_level_cm",
                      "rect_baro_level_cm",
                      "full_rect_baro_level_cm")]

baro_testing[, lapply(.SD, function(x){mean(abs(ex_abs_pressure_cm - x))}),
          by = .(serial_number, experiment),
          .SDcols = c("raw_baro_level_cm",
                      "rect_baro_level_cm",
                      "full_rect_baro_level_cm")]

full_baro <- 
  f3311[raw, on = "sample_time"]

full_baro[, `:=`(dens_air_pressure_cm = 1000 * air_pressure_cm / water_density(air_temperature_c),
                 delta_t_c = c(NA_real_, diff(air_temperature_c)),
                 delta_time = c(NA_real_, as.numeric(diff(sample_time, unit = "mins")))),
          by = .(serial_number)]

full_baro[, delta_t_c_min :=  delta_t_c / delta_time]

full_baro[, ex_abs_pressure_cm := na.approx(ex_abs_pressure_cm, sample_time, na.rm = FALSE),
          by = .(serial_number)]

full_baro[, error_range_cm := error_range(pressure_error_cm, rep(0, nrow(.SD)))]

full_baro <- 
  set_experiments(full_baro, "data/experimental_design.csv")

full_baro[baro_temperature_models, 
          rect_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
          on = "serial_number"]

full_baro[baro_delta_t_models,
          full_rect_air_pressure_cm := rect_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
          on = "serial_number"]

# Not using standard offsets from baro_training because it changes slightly each
# time the loggers are hung
full_baro[, `:=`(raw_offset_cm = mean(air_pressure_cm - ex_abs_pressure_cm, na.rm = TRUE),
                 rect_offset_cm = mean(rect_air_pressure_cm - ex_abs_pressure_cm, na.rm = TRUE),
                 full_rect_offset_cm = mean(full_rect_air_pressure_cm - ex_abs_pressure_cm, na.rm = TRUE)), 
          by = .(serial_number, experiment)]

full_baro[,`:=`(raw_baro_level_cm = air_pressure_cm - raw_offset_cm,
                rect_baro_level_cm = rect_air_pressure_cm - rect_offset_cm,
                full_rect_baro_level_cm = full_rect_air_pressure_cm - full_rect_offset_cm)]

ggplot(full_baro[experiment == "test-dat"],
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = ex_abs_pressure_cm - error_range_cm / 2,
                  ymax = ex_abs_pressure_cm + error_range_cm / 2),
              fill = "gray80") + 
  geom_line(aes(y = ex_abs_pressure_cm),
            color = "gray70") +
  geom_line(aes(y = full_rect_baro_level_cm)) +
  geom_point(data = full_baro[experiment == "test-dat" & abs(delta_t_c_min) >= 0.2],
             aes(y = full_rect_baro_level_cm),
             color = "red") +
  geom_line(aes(y = raw_baro_level_cm),
            linetype = "dotted",
            color = "blue") +
  facet_wrap(~serial_number)


full_baro[, lapply(.SD, function(x){sum(!between(x, ex_abs_pressure_cm - error_range_cm / 2, ex_abs_pressure_cm + error_range_cm / 2))}),
             by = .(serial_number, experiment),
             .SDcols = c("raw_baro_level_cm",
                         "rect_baro_level_cm",
                         "full_rect_baro_level_cm")]

full_baro[, lapply(.SD, function(x){mean(abs(ex_abs_pressure_cm - x))}),
          by = .(serial_number, experiment),
          .SDcols = c("raw_baro_level_cm",
                      "rect_baro_level_cm",
                      "full_rect_baro_level_cm")]

# Full Compensation  ------------------------------------------------------

baro <- 
  raw %>%
  .[, `:=`(delta_t_c_min = round(c(NA_real_, diff(air_temperature_c)), 1),
           delta_p_cm = c(NA_real_, diff(air_pressure_cm)))] %>% 
  .[, dens_air_pressure_cm := 1000 * air_pressure_cm / water_density(air_temperature_c)] %>% 
  set_experiments("data/experimental_design.csv")

# baro[, dens_sea_level_pressure_cm := slp(dens_air_pressure_cm, air_temperature_c)]

baro[baro_temperature_models, 
          rect_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
          on = "serial_number"]

baro[baro_delta_t_models,
          full_rect_air_pressure_cm := rect_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
          on = "serial_number"]

baro <- 
  f3311[baro, 
        .(baro_sn = serial_number, 
          sample_time, 
          air_temperature_c, 
          raw_air_pressure_cm = air_pressure_cm, 
          dens_air_pressure_cm = dens_air_pressure_cm,
          rect_air_pressure_cm, 
          ex_sea_level_pressure_cm,
          ex_air_temperature_c,
          slp_air_pressure_cm = slp(rect_air_pressure_cm, air_temperature_c),
          baro_error_cm = pressure_error_cm),
        on = "sample_time",
        nomatch = NA]

water <-
  copy(water_dat)

# water[water_tdiff_models, 
#                   rect_water_pressure_cm := water_pressure_cm - slope * (temperature_difference_c - x_intercept),
#                   on = "serial_number"]
# 
# water[water_temperature_models, 
#       rect_water_pressure_cm := rect_water_pressure_cm - slope * (water_temperature_c - x_intercept),
#       on = "serial_number"]
# 
# water[water_deltat_models, 
#       rect_water_pressure_cm := rect_water_pressure_cm - slope * (delta_t_c_min - x_intercept),
#       on = "serial_number"]
# 
# water[, `:=`(raw_offset_cm = mean(logger_error_cm) - 22500 / 826,
#              rect_offset_cm = mean(rect_water_pressure_cm - (water_depth_cm + ex_sea_level_pressure_cm))), 
#       by = .(serial_number)]
# water[, `:=`(water_level_cm = water_pressure_cm - ex_sea_level_pressure_cm - raw_offset_cm,
#              rect_water_level_cm = rect_water_pressure_cm - ex_sea_level_pressure_cm - rect_offset_cm)]
# 
water <-
  water[, .(water_sn = serial_number,
            sample_time,
            water_temperature_c,
            water_depth_cm,
            delta_wt_t_c_min = delta_t_c_min,
            raw_water_pressure_cm = water_pressure_cm,
            dens_water_pressure_cm = dens_water_pressure_cm,
            # ex_rect_water_pressure_cm = rect_water_pressure_cm,
            water_error_cm = pressure_error_cm)]

combined <- 
  merge(water,
        baro,
        allow.cartesian = TRUE,
        by = "sample_time") %>% 
  set_experiments("data/experimental_design.csv")

combined[, error_range_cm := error_range(water_error_cm, baro_error_cm)]
combined[, temperature_difference_c := air_temperature_c - water_temperature_c]
combined[, ex_temperature_difference_c := ex_air_temperature_c - water_temperature_c]

combined[water_tdiff_models, 
         `:=`(rect_water_pressure_cm = dens_water_pressure_cm - slope * (temperature_difference_c - x_intercept),
              ex_rect_water_pressure_cm = dens_water_pressure_cm - slope * (ex_temperature_difference_c - x_intercept)),
         on = c(water_sn = "serial_number")]

combined[water_temperature_models, 
         `:=`(rect_water_pressure_cm = rect_water_pressure_cm - slope * (water_temperature_c - x_intercept),
              ex_rect_water_pressure_cm = ex_rect_water_pressure_cm - slope * (water_temperature_c - x_intercept)),
         on = c(water_sn = "serial_number")]

combined[water_deltat_models, 
         `:=`(rect_water_pressure_cm = rect_water_pressure_cm - slope * (delta_wt_t_c_min - x_intercept),
              ex_rect_water_pressure_cm = ex_rect_water_pressure_cm - slope * (delta_wt_t_c_min - x_intercept)),
         on = c(water_sn = "serial_number")]



# combined[, `:=`(ex_rect_wl_slp_rect_air_cm = ex_rect_water_pressure_cm - slp_air_pressure_cm,
#                 ex_rect_wl_ex_air_cm = ex_rect_water_pressure_cm - ex_sea_level_pressure_cm,
#                 raw_wl_raw_air_cm = raw_water_pressure_cm - raw_air_pressure_cm)]
# combined[, `:=`(ex_rect_wl_slp_rect_offset_cm = mean(ex_rect_wl_slp_rect_air_cm - water_depth_cm),
#                 ex_rect_wl_ex_air_offset_cm = mean(ex_rect_wl_ex_air_cm - water_depth_cm),
#                 raw_wl_raw_offset_cm = mean(raw_wl_raw_air_cm - water_depth_cm)), 
#          by = .(baro_sn, water_sn, experiment)]
# combined[, `:=`(ex_rect_wl_slp_rect_air_cm = ex_rect_wl_slp_rect_air_cm - ex_rect_wl_slp_rect_offset_cm,
#                 ex_rect_wl_ex_air_cm = ex_rect_wl_ex_air_cm - ex_rect_wl_ex_air_offset_cm,
#                 raw_wl_raw_air_cm = raw_wl_raw_air_cm - raw_wl_raw_offset_cm)]

combined[, `:=`(raw_water_level_cm = raw_water_pressure_cm - raw_air_pressure_cm,
                rect_water_level_cm = rect_water_pressure_cm - slp_air_pressure_cm,
                ex_raw_water_level_cm = raw_water_pressure_cm - ex_sea_level_pressure_cm,
                ex_rect_water_level_cm = ex_rect_water_pressure_cm - ex_sea_level_pressure_cm)]
combined[, `:=`(raw_offset_cm = mean(raw_water_level_cm - water_depth_cm),
                rect_offset_cm = mean(rect_water_level_cm - water_depth_cm),
                ex_raw_offset_cm = mean(ex_raw_water_level_cm - water_depth_cm),
                ex_rect_offset_cm = mean(ex_rect_water_level_cm - water_depth_cm)),
         by = .(baro_sn, water_sn, experiment)]
combined[, `:=`(raw_water_level_cm = raw_water_level_cm - raw_offset_cm,
                rect_water_level_cm = rect_water_level_cm - rect_offset_cm,
                ex_raw_water_level_cm = ex_raw_water_level_cm - ex_raw_offset_cm,
                ex_rect_water_level_cm = ex_rect_water_level_cm - ex_rect_offset_cm)]

ggplot(combined[experiment == "test-dat"],
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = water_depth_cm - error_range_cm / 2,
                  ymax = water_depth_cm + error_range_cm / 2),
              fill = "gray80") +
  geom_line(aes(y = rect_water_level_cm,
                color = baro_sn)) +
  geom_line(aes(y = ex_rect_water_level_cm,
                color = "f3311"),
            linetype = "dashed") +
  geom_line(aes(y = raw_water_level_cm,
                color = baro_sn),
            linetype = "dotted") +
  facet_wrap(~water_sn)


outside_bounds <- 
  combined[, lapply(.SD, function(x){sum(!between(x, water_depth_cm - error_range_cm / 2, water_depth_cm + error_range_cm / 2))}),
           by = .(water_sn, baro_sn, experiment),
           .SDcols = patterns("water_level_cm")] %>% 
  melt(c("water_sn", "baro_sn", "experiment"),
       variable.name = "correction_method",
       value.name = "n_errors")

boxplot(n_errors ~ correction_method,
        data = outside_bounds[experiment == "test-dat"])

median_errors <- 
  combined[, lapply(.SD, function(x){median(abs(water_depth_cm - x))}),
           by = .(water_sn, baro_sn, experiment),
           .SDcols = patterns("water_level_cm")] %>% 
  melt(c("water_sn", "baro_sn", "experiment"),
       variable.name = "correction_method",
       value.name = "mad_cm")

boxplot(mad_cm ~ correction_method,
        data = median_errors[experiment == "test-dat"])
abline(h = 0.1)

# Water With Solinst Baro -------------------------------------------------

# DO 3 WATER LEVEL CORRECTIONS. ONE WITH EACH BARO AND ONE WITH THE EXTERNAL 
# DATA. ERROR IN BARO MODELS PRECLUDE ACCURATE CORRECTION OF WATER MODELS. BUT
# USING THE RECTIFIED AIR PRESSURE AND THE MOST STABLE EXPERIMENTS MAY BE ENOUGH

comp_dat <- merge(raw_water_data[, .(water_sn = serial_number, sample_time, raw_water_pressure_cm = water_pressure_cm, dens_water_pressure_cm = 1000 * water_pressure_cm / water_density(water_temperature_c), water_temperature_c, water_error_cm = pressure_error_cm)],
                  baro,
                  allow.cartesian = TRUE,
                  by = "sample_time") %>% 
  set_experiments("data/experimental_design.csv")

comp_dat[, error_range_cm := error_range(water_error_cm, baro_error_cm)]

# Align barologgers
comp_dat[baro_offsets,
         rect_air_pressure_cm := rect_air_pressure_cm - rect_offset_cm,
         on = c(baro_sn = "serial_number", "experiment")]

comp_dat[, `:=`(temperature_difference_c = air_temperature_c - water_temperature_c, 
                abs_temperature_difference_c = abs(air_temperature_c - water_temperature_c),
                slp_dens_water_pressure_cm = slp(dens_water_pressure_cm, water_temperature_c),
                slp_rect_air_pressure_cm = slp(rect_air_pressure_cm, air_temperature_c))]

comp_dat[levellogger_measurements, 
          water_depth_cm := water_depth_cm, 
          on = c(water_sn = "serial_number")]


# Calculate Raw Error
comp_dat[, error_cm := slp_dens_water_pressure_cm - (water_depth_cm + slp_rect_air_pressure_cm)]
# comp_dat[, error_cm := error_cm - mean(error_cm), by = .(water_sn)]

combined_training <- 
  comp_dat[experiment != "test-dat"]

combined_testing <- 
  comp_dat[experiment == "test-dat"]

ggplot(combined_training,
       aes(color = experiment, 
           y = error_cm,
           x = water_temperature_c)) + 
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "black", formula = y~x) +
  facet_wrap(~water_sn, scales = "free")

# Need to add a dummy column becuase the non-equi join replaces the joined column
# with two columns (one each for min and max)
combined_training[, temp_diff2 := temperature_difference_c]

vardis_means_water_tdiff <- 
  combined_training[, 
                .(vardis_mean = combined_training[experiment == "var-dis"][
                  .SD[,.(min_temp = 0.90 * min(temperature_difference_c), 
                         max_temp = 1.10 * max(temperature_difference_c))],
                  on = c("temp_diff2>=min_temp", "temp_diff2<=max_temp")][, 
                                                                        mean(error_cm)]),
                by = .(water_sn, experiment)]

combined_training[, experimental_mean := mean(error_cm),
              by = .(water_sn, experiment)]

combined_training[vardis_means_water_tdiff, 
                  vardis_mean := vardis_mean, 
                  on = c("water_sn", "experiment")]

combined_training[, experimental_error := experimental_mean - vardis_mean]


ggplot(combined_training,
       aes(x = temperature_difference_c, 
           y = error_cm)) + 
  # geom_point(aes(color = experiment)) +
  geom_point(aes(y = error_cm - experimental_error,
                 color = baro_sn),
             # color = "black",
             shape = 21) +
  geom_smooth(aes(y = error_cm - experimental_error),
              method = "lm", se = FALSE) +
  facet_wrap(~water_sn, scales = "free")

combined_training[, error_cm := error_cm - experimental_error]



ggplot(combined_training,
       aes(color = experiment, 
           y = error_cm,
           x = water_temperature_c)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black")+ 
  facet_wrap(~water_sn, scales = "free")

ggplot(combined_training,
       aes(color = experiment, 
           y = error_cm,
           x = temperature_difference_c)) + 
  geom_point() + 
  facet_wrap(~water_sn, scales = "free")


# Tdiff Models

tdiff_models <- 
  fit_correction_mod(combined_training[experiment == "stat-dis"],
                     x = "temperature_difference_c",
                     y = "error_cm",
                     by = "water_sn")

combined_training <- 
  correction_residuals(tdiff_models,
                       combined_training,
                       x = "temperature_difference_c",
                       y = "error_cm",
                       by = "water_sn")

ggplot(combined_training,
       aes(y = residuals,
           x = temperature_difference_c,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~water_sn,
             scale = "free")

ggplot(combined_training,
       aes(y = residuals,
           x = water_temperature_c,
           color = baro_sn)) + 
  geom_point() + 
  facet_wrap(~water_sn,
             scale = "free")

combined_training[tdiff_models, 
                  tdiff_dens_water_pressure_cm := dens_water_pressure_cm - slope * (temperature_difference_c - x_intercept),
                  on = "water_sn"]

combined_training[, 
                  tdiff_dens_error_cm := tdiff_dens_water_pressure_cm - (water_depth_cm + rect_air_pressure_cm)]


ggplot(combined_training,
       aes(y = tdiff_dens_error_cm,
           x = water_temperature_c)) + 
  geom_point(aes(color = experiment)) + 
  geom_smooth(data = combined_training[experiment %in% c("var-dis", "var-sim")],
              method = "lm",
              se = FALSE) +
  facet_wrap(~water_sn,
             scale = "free")

# Temp Models

temp_models <- 
  fit_correction_mod(combined_training[experiment %in% c("var-sim")],
                     x = "water_temperature_c", 
                     y = "tdiff_dens_error_cm",
                     by = "water_sn")

combined_training <- 
  correction_residuals(temp_models,
                       combined_training,
                       x = "water_temperature_c",
                       y = "tdiff_dens_error_cm",
                       by = "water_sn")

ggplot(combined_training,
       aes(y = residuals,
           x = water_temperature_c,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~water_sn,
             scale = "free")

ggplot(combined_training,
       aes(y = residuals,
           x = temperature_difference_c,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~water_sn,
             scale = "free")

combined_training[temp_models, 
                  temp_tdiff_dens_water_pressure_cm := tdiff_dens_water_pressure_cm - slope * (water_temperature_c - x_intercept),
                  on = "water_sn"]

combined_training[, 
                  temp_tdiff_dens_error_cm := temp_tdiff_dens_water_pressure_cm - (water_depth_cm + rect_air_pressure_cm)]

ggplot(combined_training,
       aes(y = temp_tdiff_dens_error_cm,
           x = temperature_difference_c,
           color = experiment)) + 
  geom_point() + 
  facet_wrap(~water_sn,
             scales = "free")

# Final Offset
combined_training[, `:=`(rect_offset_cm = mean(temp_tdiff_dens_water_pressure_cm - rect_air_pressure_cm - water_depth_cm),
                         raw_offset_cm = mean(raw_water_pressure_cm - raw_air_pressure_cm - water_depth_cm)), 
                  by = .(water_sn, baro_sn, experiment)]
combined_training[, `:=`(raw_water_level_cm = raw_water_pressure_cm - raw_air_pressure_cm - raw_offset_cm,
                         rect_water_level_cm = temp_tdiff_dens_water_pressure_cm - rect_air_pressure_cm - rect_offset_cm)]

ggplot(combined_training,
       aes(x = sample_time,
           color = baro_sn)) +
  geom_ribbon(aes(ymin = water_depth_cm - error_range_cm / 2,
                  ymax = water_depth_cm + error_range_cm / 2),
              fill = "gray80",
              color = "gray80") +
  geom_line(aes(y = rect_water_level_cm)) +
  geom_line(aes(y = raw_water_level_cm),
  linetype = "dotted") +
  facet_wrap(~water_sn,
             scales = "free")

# Test-dat

# combined_testing[tdiff_models, 
#                  rect_water_pressure_cm := slp(dens_water_pressure_cm, water_temperature_c) - slope * (temperature_difference_c - x_intercept),
#                  on = "water_sn"]


combined_testing[tdiff_models, 
                 rect_water_pressure_cm := dens_water_pressure_cm - slope * (temperature_difference_c - x_intercept),
                 on = "water_sn"]

combined_testing[temp_models, 
                 rect_water_pressure_cm := rect_water_pressure_cm - slope * (water_temperature_c - x_intercept),
                 on = "water_sn"]


combined_testing[, `:=`(rect_offset_cm = mean(rect_water_pressure_cm - slp_rect_air_pressure_cm - water_depth_cm),
                        raw_offset_cm = mean(raw_water_pressure_cm - raw_air_pressure_cm - water_depth_cm)), 
                 by = .(water_sn, baro_sn)]
combined_testing[, `:=`(raw_water_level_cm = raw_water_pressure_cm - raw_air_pressure_cm - raw_offset_cm,
                        rect_water_level_cm = rect_water_pressure_cm - slp_rect_air_pressure_cm - rect_offset_cm)]

ggplot(combined_testing,
       aes(x = sample_time,
           color = baro_sn)) +
  geom_ribbon(aes(ymin = water_depth_cm - error_range_cm / 2,
                  ymax = water_depth_cm + error_range_cm / 2),
              fill = "gray80",
              color = "gray80") +
  geom_line(aes(y = water_depth_cm),
            color = "gray50") +
  geom_line(aes(y = rect_water_level_cm)) +
  geom_line(aes(y = raw_water_level_cm),
            linetype = "dotted") +
  # geom_line(aes(y = slp_water_level_cm),
  #           linetype = "dashed") + 
  facet_wrap(~water_sn)

# Apply to combined to compare with external baro

combined[, temperature_difference_c := air_temperature_c - water_temperature_c]

combined[tdiff_models, 
         sol_rect_water_pressure_cm := dens_water_pressure_cm - slope * (temperature_difference_c - x_intercept),
         on = "water_sn"]

combined[temp_models, 
         sol_rect_water_pressure_cm := sol_rect_water_pressure_cm - slope * (water_temperature_c - x_intercept),
         on = "water_sn"]


combined[, sol_rect_wl_slp_rect_air_cm := sol_rect_water_pressure_cm - slp(rect_air_pressure_cm, air_temperature_c)]
combined[, sol_rect_wl_slp_rect_air_cm := sol_rect_wl_slp_rect_air_cm - mean(sol_rect_wl_slp_rect_air_cm - water_depth_cm),
         by = .(baro_sn, water_sn, experiment)]


ggplot(combined[experiment == "test-dat"],
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = water_depth_cm - error_range_cm / 2,
                  ymax = water_depth_cm + error_range_cm / 2),
              fill = "gray80",
              color = "gray80") +
  geom_line(aes(y = water_depth_cm),
            color = "gray50") +
  geom_line(aes(y = sol_rect_wl_slp_rect_air_cm,
                color = baro_sn)) +
  geom_line(aes(y = ex_rect_wl_ex_air_cm,
                color = "f3311")) +
  facet_wrap(~water_sn)


# Compare Performance
error_percentage <- 
  combined_testing[, .(raw = 100 * sum(!between(water_level_cm, water_depth_cm - error_range_cm / 2, water_depth_cm + error_range_cm / 2))/.N,
                       rectified = 100 * sum(!between(rect_water_level_cm, water_depth_cm - error_range_cm / 2, water_depth_cm + error_range_cm / 2))/.N),
                   by = .(water_sn, baro_sn)]

library(ggridges)
boxplot(value ~ variable, melt(error_percentage, c("water_sn", "baro_sn")))
melt(error_percentage, c("water_sn", "baro_sn")) %>% 
  ggplot(aes(x = value, y = variable)) + 
  geom_density_ridges()

error_percentage[, difference := rectified - raw]
ggplot(error_percentage[order(baro_sn, water_sn)],
       aes(x = water_sn,
           y = difference,
           ymin = difference,
           ymax = difference,
           fill = difference > 0)) + 
  geom_col(show.legend = FALSE) +
  # ylim(-25, 25) +
  ylab("Change in Erroneous Points (percentage)") +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~baro_sn) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
