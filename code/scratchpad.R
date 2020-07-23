as.numeric(st_distance(st_sfc(st_point(c(-88.58083, 47.092)), st_point(c(-88.5742812, 47.119242)), crs = "+proj=longlat +datum=WGS84"))[1,2]/1000)`

source("code/levellogger_packages.R"); source("code/levellogger_functions.R")
loadd(raw_water_data, levellogger_measurements, mesonet_response, raw_barometric_data)

theme_set(ggthemes::theme_few())
brown_green_scale <- # Inspired by NOAA maps
  c('#97601c', '#a47138', '#b08353', '#bc956e', '#c6a78a', '#d0baa6', '#d8cdc3', '#e0e0e0', '#c6d1cf', '#adc2bf', '#93b3af', '#7aa49f', '#5f9690', '#438781', '#207972')

blue_orange_scale <-
  c('#3085c7', '#4a98d4', '#6babde', '#90bee3', '#b7cfe3', '#e0e0e0', '#eac79d', '#eaae66', '#e5943a', '#dd7a15', '#d55e00')

pale_pal <- 
  c(green = "#7DA050",
    orange = "#D19648",
    teal = "#329985",
    blue = "#5982A0",
    purple = "#966283",
    red = "#9B5249")

f3311 <- mesonet_response[station == "F3311"]
f3311[, ex_pressure_error_cm := 2.7625 / qnorm(0.99)]
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

# ex_sea_level_pressure_cm should be adjusted to experimental elevation by
# ex_sea_level_pressure_cm = slp(ex_abs_pressure_cm, ex_air_temperature_c, elevation = 225)
# lm(ex_sea_level_pressure_cm ~ slp(ex_abs_pressure_cm, ex_air_temperature_c, elevation_cm/100), data = f3311)
water_dat[, raw_logger_error_cm := dens_water_pressure_cm - (water_depth_cm + ex_sea_level_pressure_cm) + 22500/826]

water_dat[, error_range_cm := error_range(pressure_error_cm, ex_pressure_error_cm)]
water_dat[, elapsed_s := as.numeric(sample_time - min(sample_time)), 
         by = .(serial_number)]


# Figure 1
fig_1 <- 
  ggplot(water_dat[experiment == "test-dat" & serial_number == "1062520", 
                   .(sample_time, 
                     temperature_difference_c, 
                     water_temperature_c,
                     error_range_cm, 
                     raw_logger_error_cm = raw_logger_error_cm - mean(raw_logger_error_cm)), 
                   by = .(serial_number)][order(sample_time)],
         aes(y = raw_logger_error_cm,
             x = water_temperature_c)) + 
  geom_hline(aes(yintercept = -error_range_cm/2),
             linetype = "dashed",
             color = "gray30") +
  geom_hline(aes(yintercept = error_range_cm/2),
             linetype = "dashed",
             color = "gray30") +
  geom_path(aes(color = temperature_difference_c)) +
  geom_point(data = water_dat[experiment == "test-dat" & serial_number == "1062520", 
                              .(sample_time, 
                                temperature_difference_c, 
                                water_temperature_c,
                                error_range_cm, 
                                raw_logger_error_cm = raw_logger_error_cm - mean(raw_logger_error_cm)), 
                              by = .(serial_number)][order(sample_time)][1]) + 
  scale_color_gradientn(colors = blue_orange_scale) + 
  theme(legend.position = "bottom")

ggsave("output/figures/figure_1.png", fig_1)

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
  water_training[, calculate_alignments(data = .SD, 
                                        full.data = water_training, 
                                        reference.experiment = "stat-sim", 
                                        x = "temperature_difference_c", 
                                        y = "raw_logger_error_cm"),
                 by = .(experiment)]

water_training[, experimental_mean := mean(raw_logger_error_cm),
              by = .(serial_number, experiment)]

water_training[, alignment_offset_cm := NULL]

water_training[statsim_means, 
               alignment_offset_cm := alignment_offset_cm, 
               on = c("serial_number", "experiment")]

water_training[, experimental_error := experimental_mean - alignment_offset_cm]

ggplot(water_training,
       aes(x = temperature_difference_c, 
           color = experiment)) + 
  geom_point(aes(y = raw_logger_error_cm),
             alpha = 0.25) +
  geom_point(aes(y = raw_logger_error_cm - experimental_error)) +
  facet_wrap(~serial_number, scales = "free")

water_training[, logger_error_cm := raw_logger_error_cm - experimental_error]

# Figure 3

fig_3 <- 
  melt(water_training[serial_number == "1062528"], 
       id.vars = c("experiment", "logger_error_cm"),
       measure.vars = c("water_temperature_c", "temperature_difference_c"),
       variable.name = "driver",
       value.name = "value_c") %>% 
  .[, driver := factor(driver, levels = c("temperature_difference_c", "water_temperature_c"))] %>% 
  ggplot(aes(y = logger_error_cm,
             x = value_c)) + 
  geom_point(aes(color = experiment, 
                 shape = experiment)) + 
  scale_color_manual(values = as.character(pale_pal)) +
  facet_wrap(~driver,
             scales = "free_x",
             as.table = FALSE,
             strip.position = "bottom") +
  geom_text(data = data.frame(x = c(-10, 8.5),
                              y = c(-1.5, -1.5),
                              label = c("A", "B"),
                              driver = c("temperature_difference_c", "water_temperature_c")),
            aes(x = x, y = y, label = label)) +
  theme(legend.position = "top")

ggsave("output/figures/figure_3.png",
       fig_3)

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

# Figure 4
fig_4 <- 
  ggplot(water_training,
         aes(y = tdiff_logger_error_cm,
             x = water_temperature_c,
             color = experiment, 
             shape = experiment)) + 
  geom_point(alpha = 0.6) + 
  facet_wrap(~serial_number, scales = "free") +
  scale_color_manual(values = as.character(pale_pal)) +
  theme(legend.position = "top")

ggsave("output/figures/figure_4.png", fig_4)

# Temp Models
water_temperature_models <- 
  fit_correction_mod(water_training,
                     x = "water_temperature_c",
                     y = "tdiff_logger_error_cm",
                     by = "serial_number")

set.seed(1234)
mcp_mods <- 
  data.table(serial_number = c("2030899", "2064734"),
             mods = list(mcp(list(tdiff_logger_error_cm ~ water_temperature_c, 
                                       ~ water_temperature_c),
                             data = water_training[serial_number == "2030899"]),
                         mcp(list(tdiff_logger_error_cm ~ 1, 
                                  ~ water_temperature_c), 
                             data = water_training[serial_number == "2064734"])))

  # water_training[serial_number %in% c("2030899", "2064734"),
  #                .(mods = list(mcp(list(tdiff_logger_error_cm ~ 1 + water_temperature_c, 
  #                                       ~ 1 + water_temperature_c), 
  #                                  data = .SD))),
  #                by = .(serial_number)]

lapply(mcp_mods$mods, summary)
mcp_mods[, change_point := sapply(mods, extract_changepoint)]
mcp_mods

# Figure 5
fig_5 <- 
  map2(mcp_mods$mods,
       mcp_mods$serial_number,
       ~plot(.x,
             q_predict = TRUE,
             lines = FALSE) +
         ggtitle(.y)) %>% 
  reduce(`+`)

ggsave("output/figures/figure_5.png",
       fig_5)

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

ggplot(water_training,
       aes(y = temp_tdiff_logger_error_cm,
           x = water_temperature_c)) + 
  geom_point() + 
  facet_wrap(~serial_number, scales = "free")


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

baro_dat[, error_range_cm := error_range(pressure_error_cm, rep(2.7625, nrow(.SD)))]

baro_training <- 
  baro_dat[experiment != "test-dat"]

baro_testing <- 
  baro_dat[experiment == "test-dat"]

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

ggplot(baro_training,
       aes(x = air_temperature_c, 
           color = experiment)) + 
  geom_point(aes(y = raw_logger_error_cm), 
             alpha = 0.25) +
  geom_point(aes(y = raw_logger_error_cm - experimental_error)) +
  facet_wrap(~serial_number, scales = "free")

baro_training[, logger_error_cm := raw_logger_error_cm - experimental_error]

# Figure 2

fig_2 <- 
  {ggplot(baro_training[serial_number == "1065861"],
          aes(x = air_temperature_c, 
              color = experiment,
              shape = experiment)) + 
      geom_point(aes(y = raw_logger_error_cm)) +
      theme(legend.position = c(0.025, 0.95),
            legend.justification = c(0, 1)) +
      ggplot(baro_training[serial_number == "1065861"],
             aes(x = air_temperature_c, 
                 color = experiment,
                 shape = experiment)) + 
      geom_point(aes(y = logger_error_cm),
                 show.legend = FALSE)} *
  theme(aspect.ratio = 1) *
  scale_color_manual(values = as.character(pale_pal[1:4]))+
  plot_annotation(tag_levels = "A")

ggsave("output/figures/figure_2.png",
       fig_2)

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
              deltat_temp_dens_air_pressure_cm := temp_dens_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
              on = "serial_number"]

baro_training[, delta_t_adjustment := NULL]

baro_training[baro_delta_t_models,
             delta_t_adjustment := slope * (delta_t_c_min - x_intercept) - intercept,
             on = "serial_number"]

baro_training[abs(delta_t_c_min) <= 0.1, delta_t_adjustment := 0]

baro_training[,cum_delta_t_adjustment := cumsum(delta_t_adjustment),
             by = .(experiment, serial_number)]

baro_training[, deltat_temp_dens_air_pressure_cm := deltat_temp_dens_air_pressure_cm - cum_delta_t_adjustment]

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
baro_training[, `:=`(raw_baro_level_cm = air_pressure_cm - raw_offset_cm,
                     rect_baro_level_cm = deltat_temp_dens_air_pressure_cm - rect_offset_cm)]

## DO THESE FINAL OFFSETS MATCH THE vardis_mean OFFSETS ABOVE?
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
  geom_line(aes(y = raw_baro_level_cm),
            linetype = "dotted",
            color = "gray40") + 
  ggtitle("Absolute Pressure Test Data") +
  facet_wrap(~ serial_number + experiment,
             scales = "free",
             nrow = 2)

baro_training[, `:=`(raw_slp_pressure_cm = slp(air_pressure_cm, air_temperature_c),
                     rect_slp_pressure_cm = slp(deltat_temp_dens_air_pressure_cm, air_temperature_c))]
baro_training[, `:=`(raw_slp_offset_cm = mean(raw_slp_pressure_cm - ex_sea_level_pressure_cm),
                     rect_slp_offset_cm = mean(rect_slp_pressure_cm - ex_sea_level_pressure_cm)), 
              by = .(serial_number, experiment)]
baro_training[, `:=`(raw_slp_pressure_cm = raw_slp_pressure_cm - raw_slp_offset_cm,
                     rect_slp_pressure_cm = rect_slp_pressure_cm - rect_slp_offset_cm)]

ggplot(baro_training,
       aes(x = sample_time,
           y = rect_slp_pressure_cm,
           color = experiment,
           fill = experiment)) +
  geom_ribbon(aes(ymin = ex_sea_level_pressure_cm - error_range_cm / 2,
                  ymax = ex_sea_level_pressure_cm + error_range_cm / 2),
              color = "gray80",
              alpha = 0.3) +
  geom_line() +
  geom_line(aes(y = raw_slp_pressure_cm),
            linetype = "dotted",
            color = "gray40") + 
  ggtitle("Sea Level Pressure Test Data") +
  facet_wrap(~ serial_number + experiment,
             scales = "free",
             nrow = 2)

# Test-dat

# baro_testing[baro_tdiff_models, 
#              rect_air_pressure_cm := dens_air_pressure_cm - slope * (temperature_difference_c - x_intercept),
#              on = "serial_number"]

baro_testing[baro_temperature_models, 
              t_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
              on = "serial_number"]

baro_testing[baro_delta_t_models,
             dt_t_air_pressure_cm := t_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
              on = "serial_number"]

baro_testing[, delta_t_adjustment := NULL]

baro_testing[baro_delta_t_models,
             delta_t_adjustment := slope * (delta_t_c_min - x_intercept) - intercept,
             on = "serial_number"]

# Only apply cumulative impact if delta T was greater than 1 unit of resolution
baro_testing[abs(delta_t_c_min) <= 0.1, delta_t_adjustment := 0]

baro_testing[, cum_delta_t_adjustment := cumsum(delta_t_adjustment),
             by = .(serial_number, experiment)]

ggplot(data =baro_testing,
       aes(x = sample_time, 
           y = cum_delta_t_adjustment)) + 
  geom_line() +
  facet_wrap(~serial_number, scales = "free")

baro_testing[, cdt_dt_t_air_pressure_cm := dt_t_air_pressure_cm - cum_delta_t_adjustment]

# baro_testing[baro_delta_t_models,
#              full_rect_air_pressure_cm := rect_air_pressure_cm - ifelse(abs(delta_t_c_min) > 0.05, slope * delta_t_c_min, 0),
#              on = "serial_number"]

# Not using standard offsets from baro_training because it changes slightly each
# time the loggers are hung
baro_testing[, `:=`(raw_offset_cm = mean(air_pressure_cm - ex_abs_pressure_cm),
                    t_offset_cm = mean(t_air_pressure_cm - ex_abs_pressure_cm),
                    dt_t_offset_cm = mean(dt_t_air_pressure_cm - ex_abs_pressure_cm),
                    cdt_dt_t_offset_cm = mean(cdt_dt_t_air_pressure_cm - ex_abs_pressure_cm)), 
             by = .(serial_number, experiment)]

baro_testing[,`:=`(raw_baro_cm = air_pressure_cm - raw_offset_cm,
                   t_baro_cm = t_air_pressure_cm - t_offset_cm,
                   dt_t_baro_cm = dt_t_air_pressure_cm - dt_t_offset_cm,
                   cdt_dt_t_baro_cm = cdt_dt_t_air_pressure_cm - cdt_dt_t_offset_cm)]

baro_testing[, `:=`(sl_pressure_cm = slp(air_pressure_cm, air_temperature_c),
                    t_sl_pressure_cm = slp(t_air_pressure_cm, air_temperature_c),
                    dt_t_sl_pressure_cm = slp(dt_t_air_pressure_cm, air_temperature_c),
                    cdt_dt_t_sl_pressure_cm = slp(cdt_dt_t_air_pressure_cm, air_temperature_c))]

baro_testing[, `:=`(slp_offset_cm = mean(sl_pressure_cm - ex_sea_level_pressure_cm),
                    t_slp_offset_cm = mean(t_sl_pressure_cm - ex_sea_level_pressure_cm),
                    dt_t_slp_offset_cm = mean(dt_t_sl_pressure_cm - ex_sea_level_pressure_cm),
                    cdt_dt_t_slp_offset_cm = mean(cdt_dt_t_sl_pressure_cm - ex_sea_level_pressure_cm)), 
             by = .(serial_number, experiment)]

baro_testing[,`:=`(raw_slp_cm = sl_pressure_cm - slp_offset_cm,
                   t_slp_cm = t_sl_pressure_cm - t_slp_offset_cm,
                   dt_t_slp_cm = dt_t_sl_pressure_cm - dt_t_slp_offset_cm,
                   cdt_dt_t_slp_cm = cdt_dt_t_sl_pressure_cm - cdt_dt_t_slp_offset_cm)]

# baro_testing[, mean((rect_air_pressure_cm - ex_abs_pressure_cm) / ex_abs_pressure_cm), by = "serial_number"]

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


{ggplot(baro_testing,
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = ex_abs_pressure_cm - error_range_cm / 2,
                  ymax = ex_abs_pressure_cm + error_range_cm / 2),
              fill = "gray80") + 
  geom_line(aes(y = ex_abs_pressure_cm),
            color = "gray70") +
  geom_line(aes(y = cdt_dt_t_baro_cm),
            linetype = "solid") +
    geom_line(aes(y = raw_baro_cm),
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
    geom_line(aes(y = cdt_dt_t_slp_cm,
              linetype = "Corrected")) +
    geom_line(aes(y = raw_slp_cm,
              linetype = "Raw")) +
  facet_wrap(~serial_number) +
  scale_linetype_manual(name = NULL,
                        values = c("dotted", "solid"),
                        breaks = c("Raw", "Corrected")) +
  scale_x_datetime(date_labels = "%H:%M") +
  theme(legend.position = "bottom")}

baro_out_of_bounds <- 
  rbind(
    baro_testing[, lapply(.SD, function(x){sum(!between(x, ex_abs_pressure_cm - error_range_cm / 2, ex_abs_pressure_cm + error_range_cm / 2))}),
                 by = .(serial_number, experiment),
                 .SDcols = patterns("baro_cm")] %>% 
      melt(c("serial_number", "experiment"),
           variable.name = "correction_method",
           value.name = "n_errors"),
    baro_testing[, lapply(.SD, function(x){sum(!between(x, ex_sea_level_pressure_cm - error_range_cm / 2, ex_sea_level_pressure_cm + error_range_cm / 2))}),
                 by = .(serial_number, experiment),
                 .SDcols = patterns("slp_cm")] %>% 
      melt(c("serial_number", "experiment"),
           variable.name = "correction_method",
           value.name = "n_errors")
  )

plot(n_errors ~ correction_method,
        data = baro_out_of_bounds)

baro_mad <- 
  rbind(
    baro_testing[, lapply(.SD, function(x){median(abs(ex_abs_pressure_cm - x))}),
                 by = .(serial_number, experiment),
                 .SDcols = patterns("baro_cm")] %>% 
      melt(c("serial_number", "experiment"),
           variable.name = "correction_method",
           value.name = "mad_cm"),
    baro_testing[, lapply(.SD, function(x){median(abs(ex_sea_level_pressure_cm - x))}),
                 by = .(serial_number, experiment),
                 .SDcols = patterns("slp_cm")] %>% 
      melt(c("serial_number", "experiment"),
           variable.name = "correction_method",
           value.name = "mad_cm"))


boxplot(mad_cm ~ correction_method,
        data = baro_mad)
# Measurement resolution
abline(h = 0.1 + (0.01 * 34.531554269))


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
  .[, `:=`(delta_at_c_min = round(c(NA_real_, diff(air_temperature_c)), 1),
           delta_p_cm = c(NA_real_, diff(air_pressure_cm)))] %>% 
  .[, dens_air_pressure_cm := 1000 * air_pressure_cm / water_density(air_temperature_c)] %>% 
  set_experiments("data/experimental_design.csv")

# baro[, dens_sea_level_pressure_cm := slp(dens_air_pressure_cm, air_temperature_c)]

# baro[baro_temperature_models, 
#           rect_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
#           on = "serial_number"]
# 
# baro[baro_delta_t_models,
#           full_rect_air_pressure_cm := rect_air_pressure_cm - slope * (delta_t_c_min - x_intercept),
#           on = "serial_number"]

baro <- 
  f3311[baro, 
        .(baro_sn = serial_number, 
          sample_time, 
          air_temperature_c, 
          delta_at_c_min,
          raw_air_pressure_cm = air_pressure_cm, 
          dens_air_pressure_cm = dens_air_pressure_cm,
          # rect_air_pressure_cm, 
          ex_sea_level_pressure_cm,
          ex_air_temperature_c,
          # slp_air_pressure_cm = slp(rect_air_pressure_cm, air_temperature_c),
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

# Compensate Baro

combined[baro_temperature_models,
         rect_air_pressure_cm := dens_air_pressure_cm - slope * (air_temperature_c - x_intercept),
         on = c(baro_sn = "serial_number")]

combined[baro_delta_t_models,
         rect_air_pressure_cm := rect_air_pressure_cm - slope * (delta_at_c_min - x_intercept),
         on = c(baro_sn = "serial_number")]

combined[, delta_t_adjustment := NULL]

combined[baro_delta_t_models,
             delta_t_adjustment := slope * (delta_at_c_min - x_intercept) - intercept,
             on = c(baro_sn = "serial_number")]

combined[abs(delta_at_c_min) <= 0.01, delta_t_adjustment := 0]

combined[, cum_delta_t_adjustment := cumsum(delta_t_adjustment),
             by = .(baro_sn, water_sn, experiment)]

combined[, rect_air_pressure_cm := rect_air_pressure_cm - cum_delta_t_adjustment]

combined[, rect_slp_cm := slp(rect_air_pressure_cm, air_temperature_c)]

ggplot(data = combined[experiment == "test-dat"],
       aes(x = sample_time, 
           y = cum_delta_t_adjustment)) + 
  geom_line(stat = "summary",
            fun = "mean") +
  facet_wrap(~baro_sn, scales = "free")

ggplot(combined[experiment == "test-dat"],
       aes(x = sample_time)) +
  geom_ribbon(aes(ymin = 0 - error_range_cm / 2,
                  ymax = 0 + error_range_cm / 2),
              fill = "gray80") + 
  geom_line(aes(y = 0),
            color = "gray70") +
  geom_line(aes(y = rect_slp_cm - mean(rect_slp_cm - ex_sea_level_pressure_cm) - ex_sea_level_pressure_cm)) +
  facet_wrap(~baro_sn)

# Compensate Water
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


# combined[, tdiff_wt_adjustment := NULL]
# 
# combined[water_tdiff_models,
#          tdiff_wt_adjustment := slope * (delta_at_c_min - x_intercept) - intercept,
#          on = c(water_sn = "serial_number")]
# 
# combined[abs(delta_at_c_min) <= 0.05, tdiff_wt_adjustment := 0]
# 
# combined[, cum_tdiff_wt_adjustment := cumsum(tdiff_wt_adjustment),
#          by = .(baro_sn, water_sn, experiment)]
# 
# ggplot(combined[experiment == "test-dat"], 
#        aes(x = sample_time,
#            y = cum_tdiff_wt_adjustment)) + 
#   geom_line(stat = "summary",
#             fun = "mean")

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
                # The rectified water level does not know how to deal with the
                # noisy air temperatures. It is applying the regular
                # compensation to the noisy temperature difference, resulting in
                # large errors. So when rect_water_leve_cm =
                # rect_water_pressure_cm - rect_slp_cm its no good. But
                # rect_wateR_level_cm = raw_water_pressure - rect_slp_cm works
                # pretty well. Need to do something similar to the
                # cumsum(delta_at) adjustment (just subtracing cum_delta_t_adjustment)
                # from the rectified water level helps. Need to find some one to
                # combined the water_tdiff and the baro_deltat models
                rect_water_level_cm = rect_water_pressure_cm - rect_slp_cm,
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

ggplot(outside_bounds,
       aes(x = correction_method,
           y = n_errors)) + 
  geom_boxplot() + 
  facet_wrap(~experiment, 
             scales = "free",
             ncol = 1)

median_errors <- 
  combined[, lapply(.SD, function(x){median(abs(water_depth_cm - x))}),
           by = .(water_sn, baro_sn, experiment),
           .SDcols = patterns("water_level_cm")] %>% 
  melt(c("water_sn", "baro_sn", "experiment"),
       variable.name = "correction_method",
       value.name = "mad_cm")

ggplot(median_errors,
       aes(x = correction_method,
           y = mad_cm)) + 
  geom_boxplot() + 
  facet_wrap(~experiment, 
             scales = "free",
             ncol = 1)
abline(h = 0.1)

# Water With Solinst Baro -------------------------------------------------

# DO 3 WATER LEVEL CORRECTIONS. ONE WITH EACH BARO AND ONE WITH THE EXTERNAL 
# DATA. ERROR IN BARO MODELS PRECLUDE ACCURATE CORRECTION OF WATER MODELS. BUT
# USING THE RECTIFIED AIR PRESSURE AND THE MOST STABLE EXPERIMENTS MAY BE ENOUGH

raw_baro_data <- 
  raw

comp_dat <- 
  merge(raw_water_data[, .(water_sn = serial_number, sample_time, raw_water_pressure_cm = water_pressure_cm, dens_water_pressure_cm = 1000 * water_pressure_cm / water_density(water_temperature_c), water_temperature_c, water_error_cm = pressure_error_cm)],
        raw_baro_data[, .(baro_sn = serial_number, sample_time, raw_air_pressure_cm = air_pressure_cm, dens_air_pressure_cm = 1000 * air_pressure_cm / water_density(air_temperature_c), air_temperature_c, air_error_cm = pressure_error_cm)],
        allow.cartesian = TRUE,
        by = "sample_time") %>% 
  .[, `:=`(delta_wt_c_min = c(NA_real_, diff(water_temperature_c)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
           delta_at_c_min = c(NA_real_, diff(air_temperature_c)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins")))),
    by = .(water_sn, baro_sn)]  %>% 
  set_experiments("data/experimental_design.csv")

comp_dat[, instrument_error_cm := combine_errors(water_error_cm, air_error_cm)]

# # Align barologgers
# baro_offsets <- 
#   baro_dat[, .(raw_offset_cm = mean(dens_air_pressure_cm - ex_abs_pressure_cm)), 
#            by = .(serial_number, experiment)]
# 
# comp_dat[baro_offsets,
#          dens_air_pressure_cm := dens_air_pressure_cm - raw_offset_cm,
#          on = c(baro_sn = "serial_number", "experiment")]

# comp_dat[baro_offsets,
#          rect_air_pressure_cm := rect_air_pressure_cm - rect_offset_cm,
#          on = c(baro_sn = "serial_number", "experiment")]

comp_dat[, `:=`(temperature_difference_c = air_temperature_c - water_temperature_c, 
                abs_temperature_difference_c = abs(air_temperature_c - water_temperature_c),
                slp_dens_water_pressure_cm = slp(dens_water_pressure_cm, water_temperature_c),
                # slp_rect_air_pressure_cm = slp(rect_air_pressure_cm, air_temperature_c),
                slp_dens_air_pressure_cm = slp(dens_air_pressure_cm, air_temperature_c))]

comp_dat[levellogger_measurements, 
          water_depth_cm := water_depth_cm, 
          on = c(water_sn = "serial_number")]

# Calculate Raw Water Level
comp_dat[, raw_water_level_cm := dens_water_pressure_cm - dens_air_pressure_cm]

# Calculate Raw Error
comp_dat[, raw_error_cm := raw_water_level_cm - water_depth_cm]

combined_training <- 
  comp_dat[experiment != "test-dat"]

combined_testing <-
  comp_dat[experiment == "test-dat"]

# It looks like 1065861 is more sensitive to noise than 1066019. This means if
# I pull out the noisy periods with some criteria about delta AT then I can 
# probably do the compensations without any reference data apart from manual
# measurements
ggplot(combined_training[abs(delta_at_c_min) < 0.1 & baro_sn == "1065861"],
         aes(x = temperature_difference_c, 
             y = raw_error_cm)) + 
    geom_point(aes(color = sample_time <= as.POSIXct("2020-06-11 12:00", tz = "EST5EDT"))) +
    facet_wrap(~water_sn, scales = "free")

# combined_training <- 
#   combined_training[baro_sn == "1066019"]

ggplot(combined_training[abs(delta_at_c_min) > 0.1],
       aes(color = experiment, 
           y = raw_error_cm,
           x = delta_at_c_min)) + 
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", color = "black", formula = y~x) +
  geom_smooth(data = combined_training[abs(delta_at_c_min) > 0],
              se = FALSE, method = "lm", color = "black", linetype = "dashed", formula = y~x) +
  facet_wrap(~water_sn, scales = "free")

full_dat_mods <- 
  fit_correction_mod(combined_training,
                     "delta_at_c_min",
                     "raw_error_cm",
                     by = c("water_sn", "baro_sn"))

combined_training <- 
  correction_residuals(full_dat_mods, 
                       combined_training,
                       "delta_at_c_min", 
                       "raw_error_cm", 
                       by = c("water_sn", "baro_sn"))

setnames(combined_training, "residuals", "dat_resid")

ggplot(combined_training,
       aes(x = air_temperature_c, 
           y = dat_resid,
           color = baro_sn)) + 
  geom_point() +
  facet_wrap(~water_sn, scales = "free")

full_dwt_mods <- 
  fit_correction_mod(combined_training,
                     "delta_wt_c_min",
                     "dat_resid",
                     by = c("water_sn", "baro_sn"))

combined_training <- 
  correction_residuals(full_dwt_mods, 
                       combined_training,
                       "delta_wt_c_min", 
                       "dat_resid", 
                       by = c("water_sn", "baro_sn"))

setnames(combined_training, "residuals", "dwt_resid")

ggplot(combined_training,
       aes(x = air_temperature_c, 
           y = dwt_resid,
           color = baro_sn)) + 
  geom_point() +
  facet_wrap(~water_sn, scales = "free")

# combined_training[experiment == "var-dis", 
#                   `:=`(water_temp_mean = mean(water_temperature_c),
#                        water_temp_sd = sd(water_temperature_c)),
#                   by = .(water_sn, baro_sn)]
# combined_training[experiment == "var-dis", 
#                   `:=`(water_temp_lwr = water_temp_mean - water_temp_sd,
#                        water_temp_upr = water_temp_mean + water_temp_sd)]

ggplot(combined_training[experiment == "var-dis"],
       aes(x = air_temperature_c,
           y = dwt_resid,
           color = between(water_temperature_c, water_temp_lwr, water_temp_upr))) +
  geom_point() +
  facet_wrap(~water_sn, scales = "free")

full_at_mods <-
  fit_correction_mod(combined_training[experiment == "var-dis"],
                     "air_temperature_c",
                     "dwt_resid",
                     by = c("water_sn", "baro_sn"))

combined_training <-
  correction_residuals(full_at_mods,
                       combined_training,
                       "air_temperature_c",
                       "dwt_resid",
                       by = c("water_sn", "baro_sn"))
 
setnames(combined_training, "residuals", "at_resid")

full_wt_mods <-
  fit_correction_mod(combined_training[experiment %in% c("var-sim")],
                     "water_temperature_c",
                     "at_resid",
                     by = c("water_sn", "baro_sn"))

combined_training <- correction_residuals(full_wt_mods,
                                          combined_training,
                                          "water_temperature_c",
                                          "at_resid",
                                          by = c("water_sn", "baro_sn"))

setnames(combined_training, "residuals", "wt_resid")

ggplot(combined_training,
       aes(x = temperature_difference_c,
           y = wt_resid,
           color = experiment)) +
  geom_point() +
  facet_wrap(~water_sn, scales = "free")

# Compensate

combined_training[, rect_water_level_cm := raw_water_level_cm]

combined_training[full_dat_mods,
                  rect_water_level_cm := rect_water_level_cm - slope * (delta_at_c_min - x_intercept),
                  on = c("water_sn", "baro_sn")]

combined_training[full_dwt_mods,
                  rect_water_level_cm := rect_water_level_cm - slope * (delta_wt_c_min - x_intercept),
                  on = c("water_sn", "baro_sn")]

combined_training[full_at_mods, 
                  rect_water_level_cm := rect_water_level_cm - slope * (air_temperature_c - x_intercept),
                  on = c("water_sn", "baro_sn")]

combined_training[full_wt_mods, 
                 rect_water_level_cm := rect_water_level_cm - slope * (water_temperature_c - x_intercept),
                 on = c("water_sn", "baro_sn")]

# Final Offset
combined_training[, `:=`(rect_offset_cm = mean(rect_water_level_cm - water_depth_cm),
                         raw_offset_cm = mean(raw_water_level_cm - water_depth_cm)), 
                  by = .(water_sn, baro_sn, experiment)]
combined_training[, `:=`(raw_water_level_cm = raw_water_level_cm - raw_offset_cm,
                         rect_water_level_cm = rect_water_level_cm - rect_offset_cm)]

ggplot(combined_training,
       aes(x = sample_time,
           color = baro_sn)) +
  geom_ribbon(aes(ymin = water_depth_cm - instrument_error_cm * qnorm(0.99),
                  ymax = water_depth_cm + instrument_error_cm * qnorm(0.99)),
              fill = "gray80",
              color = "gray80") +
  geom_line(aes(y = rect_water_level_cm)) +
  # geom_line(aes(y = raw_water_level_cm),
  # linetype = "dotted") +
  facet_wrap(~water_sn,
             scales = "free") +
  ggtitle("td")

# Test-dat

combined_testing[, rect_water_level_cm := raw_water_level_cm]

combined_testing[full_dat_mods,
                  rect_water_level_cm := rect_water_level_cm - slope * (delta_at_c_min - x_intercept),
                  on = c("water_sn")]

combined_testing[full_dwt_mods,
                  rect_water_level_cm := rect_water_level_cm - slope * (delta_wt_c_min - x_intercept),
                  on = c("water_sn")]

combined_testing[full_at_mods, 
                  rect_water_level_cm := rect_water_level_cm - slope * (air_temperature_c - x_intercept),
                  on = c("water_sn")]

combined_testing[full_wt_mods, 
                  rect_water_level_cm := rect_water_level_cm - slope * (water_temperature_c - x_intercept),
                  on = c("water_sn")]

# Final Offset
combined_testing[, `:=`(rect_offset_cm = mean(rect_water_level_cm - water_depth_cm),
                         raw_offset_cm = mean(raw_water_level_cm - water_depth_cm)), 
                  by = .(water_sn, baro_sn, experiment)]

combined_testing[, `:=`(raw_water_level_cm = raw_water_level_cm - raw_offset_cm,
                        rect_water_level_cm = rect_water_level_cm - rect_offset_cm)]

ggplot(combined_testing,
       aes(x = sample_time,
           color = baro_sn)) +
  geom_ribbon(aes(ymin = water_depth_cm - instrument_error_cm * qnorm(0.99),
                  ymax = water_depth_cm + instrument_error_cm * qnorm(0.99)),
              fill = "gray80",
              color = "gray80") +
  geom_line(aes(y = water_depth_cm),
            color = "gray50") +
  geom_line(aes(y = rect_water_level_cm)) +
  geom_line(aes(y = raw_water_level_cm),
            linetype = "dotted") +
  facet_wrap(~water_sn,
             scales = "free")

# Apply to combined to compare with external baro compensations

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
  combined_testing[, .(raw = 100 * sum(!between(raw_water_level_cm, water_depth_cm - instrument_error_cm * qnorm(0.99), water_depth_cm + instrument_error_cm * qnorm(0.99)))/.N,
                       rectified = 100 * sum(!between(rect_water_level_cm, water_depth_cm - instrument_error_cm * qnorm(0.99), water_depth_cm + instrument_error_cm * qnorm(0.99)))/.N),
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

combined_mad <- 
  combined_testing[, .(raw = median(abs(raw_water_level_cm - water_depth_cm)),
                       rectified = median(abs(rect_water_level_cm - water_depth_cm)),
                       raw_pct = 100 * median(abs(raw_water_level_cm - water_depth_cm))/mean(water_depth_cm),
                       rectified_pct = 100 * median(abs(rect_water_level_cm - water_depth_cm))/mean(water_depth_cm)),
                   by = .(water_sn, baro_sn)]

boxplot(value ~ variable, 
        data = melt(combined_mad, id.vars = c("water_sn", "baro_sn"), measure.vars = c("raw", "rectified")))

boxplot(value ~ variable, 
        data = melt(combined_mad, id.vars = c("water_sn", "baro_sn"), measure.vars = c("raw_pct", "rectified_pct")))
