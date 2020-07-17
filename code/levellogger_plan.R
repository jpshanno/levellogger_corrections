# Read in all data
# Adjust for density
# Model & compensate barometric data using only var-dis. var-sim introduces
#   larger temperature differences between the Mesonet and Solinst air
#   temperatures. Can use var-dis & var-sim. There is a temperature difference
#   issue between the Solinst & Mesonet data. If abs(t-diff) < 2 then the effect
#   of t-diff is limited enough
# Use stat-sim to calculate offset for levelloggers
# Use stat-sim and var-sim to build water level models # VAR SIM MAY HAVE PROBLEMS with being in water
# Apply water level modes and use stat-dis and var-dis to build temperature 
#  difference models
# Run all models on test-dat and see how the correction turns out
# Test using external barometric data on var-sim & var-dis to check water temp
#   models, and then compare with test-dat

# There is a change in air pressure on the morning on June 10 that is captured
# by F3311, but not KCMX. It is the cause of a major problem. Might have to use
# F3311 because of the match of the data, even if there are fewer points.
# dy_graph(mesonet_response[order(station, sample_time), diff := c(NA_real_, diff(ex_abs_pressure_cm))], "diff", grouping = "station") %>% dyOptions(useDataTimezone = TRUE, connectSeparatedPoints = TRUE)
# dy_graph(baro_dat[experiment == "var-sim" & station == "KCMX"][order(serial_number, sample_time), diff := c(NA_real_, diff(air_pressure_cm))], "diff", grouping = "serial_number")
# It also shows up looking at baro_error ~ air_temperature for KCMX vs F3311. 
# Weirdly it appears in just 1065861

plan <- 
  drake_plan(
    
    # Load Data ---------------------------------------------------------------
    # The air pressure loggers were found wet at the end of the var-sim
    # experiment. After looking at the data extensively I determined that it 
    # must have fallen into the water at approximately 6:00 on 6/10/2020. The 
    # clearest way to see this was too look at a plot of water level error by 
    # water_temperature for just the var-sim experiment. Look at the plot in
    # output/figures/problem-baro-plot.png made using corrected air pressures.
    # with uncorrected water pressures. I compared that figure to one made with
    # the Mesonet altimeter data and found that the artifacts were removed
    
    levellogger_measurements = 
      fread(file_in("data/levellogger_measurements.csv"),
            colClasses = c("character", "character", "integer", "numeric", "numeric")),
    
    raw_water_data = 
      read_xle_data(levellogger_measurements, "water") %>% 
      .[, water_pressure_cm := align_loggers(., water_pressure_cm)] %>% 
      .[, dens_water_pressure_cm := adjust_density(water_pressure_cm, water_temperature_c)],
    
    raw_barometric_data = 
      read_xle_data(levellogger_measurements, "air") %>% 
      .[!between(sample_time, 
                 as.POSIXct("2020-06-10 06:00", tz = "EST5EDT"), 
                 as.POSIXct("2020-06-10 12:30", tz = "EST5EDT"))] %>%
      .[, `:=`(delta_t_c = c(NA_real_, diff(air_temperature_c)),
               delta_p_cm = c(NA_real_, diff(air_pressure_cm)),
               rolling_temp_c = frollmean(air_temperature_c, n = 720, align = "center"))] %>% 
      set_experiments("data/experimental_design.csv") %>%
      # Could align them to account for differences in hanger height (which I 
      # thought was equal). Though at times it did shift. Perhaps should align
      # them by experiment. There is something 
      # .[, air_pressure_cm := align_loggers(., air_pressure_cm, group = "serial_number")] %>%
      .[, dens_air_pressure_cm := adjust_density(air_pressure_cm, air_temperature_c)],
    
    logger_metadata = 
      get_metadata(raw_barometric_data, raw_water_data),
    
    mesonet_response = 
      fetch_external_baro(times = c(raw_barometric_data$sample_time)),
    
    # Compensate Barometric Data ----------------------------------------------
    # baro_temp_mod_data =
    #   merge(raw_barometric_data[, .(experiment, serial_number, sample_time, 
    #                                 dens_air_pressure_cm, air_temperature_c)],
    #         mesonet_response[, .(station, sample_time, ex_altimeter_cm,
    #                              ex_air_temperature_c)],
    #         by = "sample_time") %>%
    #   .[experiment %in% c("var-dis") & station %in% c("KCMX")] %>%
    #   .[.[, .(mean_error_cm = mean(ex_altimeter_cm - dens_air_pressure_cm)), by = .(serial_number)],
    #     on = "serial_number"] %>% 
    #   .[, `:=`(baro_error_cm = ex_altimeter_cm - dens_air_pressure_cm + mean_error_cm,
    #            temperature_difference_c = ex_air_temperature_c - air_temperature_c,
    #            abs_temperature_difference_c = abs(ex_air_temperature_c - air_temperature_c))],
    # 
    # # baro_temp_mod_data %>%
    # #   split(by = "serial_number") %>%
    # #   map(~ggplot(data = .x,
    # #               aes(x = air_temperature_c,
    # #                   y = baro_error_cm,
    # #                   color = abs(temperature_difference_c))) +
    # #         geom_point() +
    # #         ggtitle(unique(.x[["serial_number"]]))) %>%
    # #   reduce(`/`)
    # 
    # # Fit unrestricted model
    # 
    # # fit_correction_mod(baro_temp_mod_data, "air_temperature_c", "baro_error_cm", "serial_number,station")
    # # fit_correction_mod(baro_temp_mod_data, "air_temperature_c", "baro_error_cm", "serial_number")
    # 
    # unrestricted_baro_models = 
    #   fit_correction_mmod(baro_temp_mod_data, 
    #                      x = "air_temperature_c", 
    #                      y = "baro_error_cm",
    #                      group = "station", 
    #                      by = "serial_number"),
    
    # Check fit
    # correction_residuals(unrestricted_baro_models,
    #                      baro_temp_mod_data,
    #                      x = "air_temperature_c",
    #                      y = "baro_error_cm",
    #                      by = "serial_number") %>%
    #   ggplot(aes(x = air_temperature_c,
    #              y = baro_error_cm,
    #              color = abs(temperature_difference_c),
    #              shape = serial_number)) +
    #   geom_point() +
    #   geom_line(aes(y = baro_error_cm - residuals),
    #             color = "gray20")
    
    # Find Temperature Difference Threshold
    
    # The mcp package does not allow for random terms in the slope, only the
    # intercept or the variance. The two loggers have opposite responses to
    # large temperature differences, one shows large positive resids and one
    # shows large negative resids. The options are separate models by serial
    # number and find the change points, a single model looking for change point
    # of absolute value of the standardized resids, or a single model looking at
    # the change in variance rather than the change in slope. Quick checks yield
    # different change points (3.3 for abs(std.resid) 4.2 for var(std.resid)).
    # The resultant models are very similar diff(slope) < 0.0005; diff(intercept)
    # < 0.005. Going with variance because it makes sense intuitively when 
    # looking at standardized residuals, which should always have a variance of
    # 1. Additionally the direction of error seems to vary with logger, and a 
    # model that incorporates sigma() captures this even if the loggers are split
    # relatively even between positve and negative responses. A slope-dependent
    # model does not capture this
    
    # correction_residuals(unrestricted_baro_models,
    #                      baro_temp_mod_data,
    #                      x = "air_temperature_c",
    #                      y = "baro_error_cm",
    #                      by = "serial_number") %>%
    #   ggplot(aes(x = air_temperature_c,
    #              y = residuals,
    #              color = abs(temperature_difference_c))) +
    #   geom_point() +
    #   facet_wrap(~serial_number)
    # 
    # correction_residuals(unrestricted_baro_models,
    #                      baro_temp_mod_data,
    #                      x = "air_temperature_c",
    #                      y = "baro_error_cm",
    #                      by = "serial_number") %>%
    #   ggplot(aes(x = temperature_difference_c,
    #              y = abs_residuals,
    #              color = serial_number)) +
    #   geom_point() +
    #   facet_wrap(~serial_number)
    
    # The changepoint comes from temperature difference not absolute air temp
    # In this case 14.5 comes from 
    # {set.seed(1234)
    #   mcp(list(residuals ~ 1 + air_temperature_c,
    #            1 + (1|serial_number) ~ 1 + sigma(1)),
    #       data = .,
    #       cores = 3)} 
    # and 2.31 is taken from the model below. That means that the low-temp
    # difference model is the true model and the temp difference adjustment is
    # not necessary during corrections! The plot below shows the var-dis
    # experiment data with a big blob of low air temperature values sitting down
    # below the regression line. Looking at the time of those data points show
    # that they are from the earliest part of the experiment, where the loggers
    # were still adjusting to temperature the air temperature, and may have been
    # influenced by that fact, resulting in a larger temperature difference
    
    # ggplot(baro_temp_mod_data[air_temperature_c > 14.5],
    #        aes(x = air_temperature_c,
    #            y = baro_error_cm,
    #            color = abs_temperature_difference_c > 2.31)) +
    #   geom_point() +
    #   facet_wrap(~serial_number)
    
    # baro_changepoint_models =
    #   correction_residuals(unrestricted_baro_models,
    #                        baro_temp_mod_data,
    #                        x = "air_temperature_c",
    #                        y = "baro_error_cm",
    #                        by = "serial_number") %>%
    #   {set.seed(1234)
    #     mcp(list(residuals ~ 1 + abs_temperature_difference_c,
    #              1 ~ 1 + abs_temperature_difference_c + sigma(1)),
    #         data = .,
    #         cores = 3)},
    # 
    # # Fit restricted models
    # 
    # baro_models = 
    #   fit_correction_mod(baro_temp_mod_data[abs_temperature_difference_c < extract_changepoint(baro_changepoint_models)], 
    #                      x = "air_temperature_c", 
    #                      y = "baro_error_cm",
    #                      # group = "station", 
    #                      by = "serial_number") %>% 
    #   .[, change_point := extract_changepoint(baro_changepoint_models)],
    
    # baro_high_temp_models = 
    #   fit_correction_mmod(baro_temp_mod_data[air_temperature_c > extract_changepoint(baro_changepoint_models)], 
    #                      x = "air_temperature_c", 
    #                      y = "baro_error_cm",
    #                      group = "station", 
    #                      by = "serial_number") %>% 
    #   .[, change_point := extract_changepoint(baro_changepoint_models)],
 
    # Recheck residuals
    # rbind(correction_residuals(baro_low_temp_models,
    #                                baro_temp_mod_data[air_temperature_c <= extract_changepoint(baro_changepoint_models)],
    #                                x = "air_temperature_c",
    #                                y = "baro_error_cm",
    #                                by = "serial_number"),
    #           correction_residuals(baro_high_temp_models,
    #                                baro_temp_mod_data[air_temperature_c > extract_changepoint(baro_changepoint_models)],
    #                                x = "air_temperature_c",
    #                                y = "baro_error_cm",
    #                                by = "serial_number")) %>%
    # correction_residuals(baro_models,
    #                      baro_temp_mod_data[abs_temperature_difference_c < extract_changepoint(baro_changepoint_models)],
    #                      x = "air_temperature_c",
    #                      y = "baro_error_cm",
    #                      by = "serial_number") %>%
    #   ggplot(aes(x = air_temperature_c,
    #              y = residuals)) +
    #   geom_point() +
    #   facet_wrap(~serial_number)
  
    # Need to check that restricted model from here matches long-term model of deployed
    # loggers
    
    # Save final barometric data
    
    
    # barometric_data =
    #   # rbind(raw_barometric_data[baro_low_temp_models[, -c("mod")], 
    #   #                           .(serial_number, sample_time, slope = slope, x_intercept = x_intercept),
    #   #                           on = c("serial_number", "air_temperature_c <= change_point")], 
    #   #       raw_barometric_data[baro_high_temp_models[, -c("mod")], 
    #   #                           .(serial_number, sample_time, slope = slope, x_intercept = x_intercept),
    #   #                           on = c("serial_number", "air_temperature_c > change_point")]) %>% 
    #   # merge(raw_barometric_data, 
    #   #       by = c("serial_number", "sample_time")) %>% 
    #   raw_barometric_data[baro_models, on = "serial_number"] %>% 
    #   .[.[mesonet_response[station == "KCMX"], on = "sample_time", nomatch = 0][experiment == "var-dis", .(mean_error_cm = mean(ex_altimeter_cm - dens_air_pressure_cm)), by = .(serial_number)],
    #     on = "serial_number"] %>% 
    #   .[, temp_dens_air_pressure_cm := dens_air_pressure_cm - (slope * (x_intercept - air_temperature_c)) + mean_error_cm] %>% 
    #   # .[, temp_dens_air_pressure_cm := align_loggers(., temp_dens_air_pressure_cm)] %>% 
    #   # 0.13 from below
    #   # .[, off_temp_dens_air_pressure_cm := temp_dens_air_pressure_cm + 0.13] %>% 
    #   .[, .(baro_serial_number = serial_number,
    #         sample_time,
    #         experiment,
    #         air_temperature_c,
    #         air_pressure_cm,
    #         dens_air_pressure_cm,
    #         temp_dens_air_pressure_cm,
    #         # off_temp_dens_air_pressure_cm,
    #         air_pressure_error_cm = pressure_error_cm,
    #         air_temperature_error_c = temperature_error_c)],

    # There is no need for the offset becuase doing slope * (x_intercept - air_temperature_c) accounts for any offset
    # It would be worth aligning the loggers though
    # Offsets
    #  barometric_data[mesonet_response, on = "sample_time"][experiment == "var-dis" &  abs(ex_air_temperature_c - air_temperature_c) < extract_changepoint(baro_changepoint_models), .(offset = mean(ex_altimeter_cm - temp_dens_air_pressure_cm)), by = "serial_number"]
    
    # # Checking Data
    # # Looks good
    # loadd(barometric_data, mesonet_response, baro_changepoint_models)
    # test_dat <-
    #   merge(barometric_data,
    #         mesonet_response,
    #         by = "sample_time") %>%
    #   .[abs(ex_air_temperature_c - air_temperature_c) < extract_changepoint(baro_changepoint_models)] %>%
    #   .[, final_error := ex_altimeter_cm - temp_dens_air_pressure_cm]
    # 
    # # There's an error in var-dis 1065861 drops off of (0,1) line ~1001 cm
    # ggplot(test_dat,
    #        aes(x = ex_altimeter_cm,
    #            y = temp_dens_air_pressure_cm,
    #            color = baro_serial_number)) +
    #   geom_point(alpha = 0.5) +
    #   geom_abline(aes(intercept = 0, slope = 1)) +
    #   facet_wrap(~experiment,
    #              scales = "free")
    # # Some heteroskedacity shows up in the final dataset, which makes sense since
    # # it is two different models
    # ggplot(test_dat,
    #        aes(x = air_temperature_c,
    #            y = final_error,
    #            color = baro_serial_number)) +
    #   geom_point() +
    #   facet_wrap(~experiment,
    #              scales = "free")
    # # Offsets
    # # Look up to see if there is a known increase in barometric pressure
    # # in basements
    # ggplot(test_dat,
    #        aes(x = baro_serial_number,
    #            y = final_error,
    #            color = experiment)) +
    #   geom_boxplot()
    
    
    # stat-dis and stat-sim fall along the line of the models developed above. 
    # I tested and the model slopes and intercepts are similar, but the change
    # point shifts dramatically. This isn't too surprising as the variance 
    # associated with the indoor experiments is much higher than the outdoor
    # experiments

    # dat <-
    #   barometric_data[experiment != "test-dat",
    #                   c(.SD,
    #                     list(baro_error_cm = ex_altimeter_cm - dens_air_pressure_cm,
    #                          temperature_difference_c = ex_air_temperature_c - air_temperature_c))]
    # 
    # dat_mod_sigma <- 
    #   fit_correction_mod(dat,
    #                      x = "air_temperature_c",
    #                      y = "baro_error_cm",
    #                      group = "baro_serial_number") %>%
    #   correction_residuals(.,
    #                        dat,
    #                        x = "air_temperature_c",
    #                        y = "baro_error_cm",
    #                        group = "baro_serial_number") %>%
    #   # This is using slope instead of variance as changepoint
    #   mcp(list(std_residuals ~ temperature_difference_c,
    #            ~ temperature_difference_c + sigma(1)),
    #       data = .)
    # 
    # plot(dat_mod, q_predict = TRUE, lines = FALSE)
    # 
    # fit_correction_mod(dat[temperature_difference_c > extract_changepoint(dat_mod)],
    #                    x = "air_temperature_c",
    #                    y = "baro_error_cm",
    #                    group = "baro_serial_number")
    
    # ggplot(barometric_data,
    #        aes(x = air_temperature_c,
    #            y = ex_altimeter_cm - dens_air_pressure_cm,
    #            color = experiment)) +
    #   geom_point() +
    #   facet_wrap(~baro_serial_number, scales = "free")
    
    
    # Calculate Logger Offsets ------------------------------------------------
    # merge(raw_barometric_data[, .(serial_number, sample_time, dens_air_pressure_cm, air_temperature_c)],
    #       mesonet_response[, .(station, sample_time, ex_altimeter_cm)],
    #       by = "sample_time") %>% 
    #   set_experiments("data/experimental_design.csv") %>% 
    #   .[experiment == "var-dis"] %>% 
    #   .[, baro_error_cm := ex_altimeter_cm - dens_air_pressure_cm] %>% 
    #   .[,.(mod = list(lmer(baro_error_cm ~ air_temperature_c + (air_temperature_c | station),
    #                        control = lmerControl(optimizer = "Nelder_Mead")))),
    #     by = .(serial_number)] %>% 
    #   .[, `:=`(slope = vapply(mod, function(x){fixef(x)[2]}, numeric(1)),
    #            intercept = vapply(mod, function(x){fixef(x)[1]}, numeric(1)))] %>% 
    #   .[, x_intercept := -intercept / slope],
    

    # Compensate for Water Temperature ----------------------------------------
    
    #,
    # combined = 
    #   merge_pressures(raw_water_data, raw_barometric_data) %>% 
    #   .[mesonet_response, on = "sample_time", nomatch = 0] %>% 
    #   adjust_density() %>% 
    #   set_experiments(file_in("data/experimental_design.csv")) %>% 
    #   set(j = "naive_water_level_cm",
    #       value = .[, water_pressure_cm - air_pressure_cm]) %>%
    #   set(j = "ex_naive_water_level_cm",
    #       value = .[, water_pressure_cm - ex_air_pressure_cm]) %>%
    #   set(j = "ex_densW_water_level_cm",
    #       value = .[, dens_water_pressure_cm - ex_air_pressure_cm]) %>%
    #   set(j = "densA_water_level_cm",
    #       value = .[, water_pressure_cm - dens_air_pressure_cm]) %>% 
    #   set(j = "densAW_water_level_cm",
    #       value = .[, dens_water_pressure_cm - dens_air_pressure_cm])
  )
