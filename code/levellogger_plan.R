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

# Load Data ---------------------------------------------------------------
# The air pressure loggers were found wet at the end of the var-sim
# experiment. After looking at the data extensively I determined that it 
# must have fallen into the water at approximately 6:00 on 6/10/2020. The 
# clearest way to see this was too look at a plot of water level error by 
# water_temperature for just the var-sim experiment. Look at the plot in
# output/figures/problem-baro-plot.png made using corrected air pressures.
# with uncorrected water pressures. I compared that figure to one made with
# the Mesonet altimeter data and found that the artifacts were removed

plan <- 
  drake_plan(
    
    # Read in manual measurements of cable length and water depth
    levellogger_measurements = 
      read_manual_measurements(file_in("data/levellogger_measurements.csv")),
    
    # Load xle data from water transducers
    raw_water_data =
      read_xle_data(levellogger_measurements, 
                    "water"),
    
    # Load xle data from barometric transducers
    raw_barometric_data =
      read_xle_data(levellogger_measurements, 
                    "air"),
    
    # L
    instrument_metadata =
      get_metadata(raw_barometric_data, 
                   raw_water_data),
     
    combined_data = 
      combine_datasets(raw_water_data, 
                       raw_barometric_data, 
                       levellogger_measurements, 
                       file_in("data/experimental_design.csv")),
                       
    bootstrap_models = 
      fit_models(combined_data),
     
    fitted = 
      calculate_fitted_values(bootstrap_models, combined_data)
    
    # predicted =
    #   calculate_predicted_values(boostrap_models, combined_data)
    
    # All models predicted for their own data
    # uncertainty_rdses =   
    #   {pred_mat <-
    #     bootstrap_models[, .(B = list(matrix(c(intercept,
    #                                            air_temperature_c,
    #                                            water_temperature_c,
    #                                            delta_at_01c_min),
    #                                          dimnames = list(c("i", 
    #                                                            "b_at",
    #                                                            "b_wt",
    #                                                            "b_dat"),
    #                                                          NULL),
    #                                          nrow = 4)),
    #                          sigma = list(matrix(rep(sigma, nobs),
    #                                              dimnames = list(NULL, "sigma"),
    #                                              ncol = 1))),
    #                      keyby = .(water_sn, baro_sn, experiment, rep)]
    #   #   
    #   #   
    #   x_mat <- Reduce(rbind, combined_data)[,  .(S = list(matrix(sample_time, 
    #                                                              dimnames = list(NULL, "sample_time"),
    #                                                              ncol = 1)),
    #                                              X = list(matrix(c(rep(1, .N),
    #                                                                air_temperature_c,
    #                                                                water_temperature_c,
    #                                                                delta_at_01c_min),
    #                                                              dimnames = list(NULL, c("intercept", 
    #                                                                                      "air_temperature_c",
    #                                                                                      "water_temperature_c",
    #                                                                                      "delta_at_01c_min")),
    #                                                              ncol = 4))),
    #                                         keyby = .(water_sn, baro_sn, experiment)]
    #   
    #   pred_mat[x_mat,
    #            `:=`(S = i.S,
    #                 X = i.X)]
    #   
    #   rm(bootstrap_models, x_mat)
    #   # matrix(rnorm(nrow(i), y_hat, s), dimnames = list(NULL, "predicted_error_cm"))}     
    #   pred_mat[, E_hat := Map(function(x, b, s){matrix(rnorm(n = nrow(x), 
    #                                                          mean = x %*% b, 
    #                                                          sd = s), 
    #                                                    ncol = 1,
    #                                                    dimnames = list(NULL, 
    #                                                                    "predicted_error_cm"))}, 
    #                           x = X,
    #                           b = B,
    #                           s = sigma)]
    #   
    #   for(Ex in c("var-sim", "var-dis", "stat-sim", "stat-dis", "test-dat")){
    #     saveRDS(pred_mat[experiment == Ex], 
    #             paste0("output/tabular/bootstrap_prediction_matrices_", Ex, ".rds"))}
    #   
    #   # paste0("output/tabular/bootstrap_prediction_matrices_", c("var-sim", "var-dis", "stat-sim", "stat-dis", "test-dat"), ".rds")
    #   }
    
    # predicted = {
    #   posterior_pred <-
    #     disk.frame("output/tabular/bootstrap_predictions_test_data.df")
    # 
    #   predictions <-
    #     posterior_pred[, .(predicted_error_cm = quantile(predicted_error_cm, probs = 0.5),
    #                        lower_bound_cm = quantile(predicted_error_cm, probs = 0.025),
    #                        upper_bound_cm = quantile(predicted_error_cm, probs = 0.975)),
    #                    by = .(water_sn, baro_sn, sample_time),
    #                    keep = c("water_sn", "baro_sn", "sample_time", "predicted_error_cm")]
    #   
    #   predictions[, `:=`(water_sn = as.character(water_sn),
    #                      baro_sn = as.character(baro_sn))]
    #   
    #   setkey(combined_data$testing, "water_sn", "baro_sn", "sample_time")
    #   setkey(predictions, "water_sn", "baro_sn", "sample_time")
    #   
    #   predictions[combined_data$testing,
    #               `:=`(experiment = i.experiment,
    #                    water_depth_cm = i.water_depth_cm,
    #                    raw_water_level_cm = i.raw_water_level_cm,
    #                    raw_error_cm = i.raw_error_cm,
    #                    instrument_error_cm = i.instrument_error_cm)]
    #   
    #   predictions[, rect_water_level_cm := raw_water_level_cm - predicted_error_cm]
    #   predictions[, `:=`(raw_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
    #                      rect_water_level_cm = rect_water_level_cm - mean(rect_water_level_cm - water_depth_cm)),
    #               by = .(water_sn, baro_sn, experiment)]
    # }
      
    # Runs out of memory. Ran once interactively to build disk.frame. Need to 
    # sort this out
    # {
    #   predictors <- 
    #     dcast(rbind(bootstrap_models$fit_metrics[, .(model_id, water_sn, baro_sn, rep,
    #                                                  term = "sigma", estimate = sigma)], 
    #                 bootstrap_models$coefficients[, .(model_id, water_sn, baro_sn, 
    #                                                   rep, term, estimate)], 
    #                 fill = TRUE),
    #           model_id + water_sn + baro_sn + rep ~ term,
    #           value.var = "estimate")
    #   
    #   predictors[combined_data$testing[, .(nobs = .N), .(water_sn, baro_sn)],
    #              nobs := i.nobs,
    #              on = c("water_sn", "baro_sn")]
    #   
    #   setnames(predictors, 
    #            c("(Intercept)", "rep"),
    #            c("intercept", "replication"))
    #   
    #   predictors[, `:=`(air_temperature_c = nafill(air_temperature_c, "const", 0),
    #                     water_temperature_c = nafill(water_temperature_c, "const", 0),
    #                     delta_at_01c_min = nafill(delta_at_01c_min, "const", 0))]
    #   
    # Test data predicted with training models
    # test_predictions = 
    
    # ggplot(bootstrap_models$coefficients[term == "air_temperature_c"],
    #        aes(x = estimate,
    #            y = water_sn,
    #            fill = baro_sn)) +
    #   stat_halfeye(alpha = 0.5)
    # ggplot(bootstrap_models$fit_metrics, 
    #        aes(x = rmse, 
    #            y = water_sn, 
    #            fill = baro_sn)) + 
    #   stat_halfeye(alpha = 0.5)
    
    # Compensate Water Levels -------------------------------------------------
    
    # delta_air_temp_mods =
    #   fit_correction_mod(data = combined_data[experiment != "test-dat"],
    #                      x = "delta_at_c_min",
    #                      y = "raw_error_cm",
    #                      by = c("water_sn", "baro_sn")),
    # 
    # delta_water_temp_mods =
    #   correction_residuals(mods = delta_air_temp_mods, 
    #                        original.data = combined_data[experiment != "test-dat"],
    #                        x = "delta_at_c_min", 
    #                        y = "raw_error_cm", 
    #                        by = c("water_sn", "baro_sn")) %>% 
    #   fit_correction_mod(x = "delta_wt_c_min",
    #                      y = "residuals",
    #                      by = c("water_sn", "baro_sn")),
    # 
    # air_temp_mods =
    #   correction_residuals(mods = delta_air_temp_mods, 
    #                        original.data = combined_data[experiment != "test-dat"],
    #                        x = "delta_at_c_min", 
    #                        y = "raw_error_cm", 
    #                        by = c("water_sn", "baro_sn")) %>% 
    #   correction_residuals(mods = delta_water_temp_mods, 
    #                        original.data = .,
    #                        x = "delta_wt_c_min", 
    #                        y = "residuals", 
    #                        by = c("water_sn", "baro_sn")) %>% 
    #   fit_correction_mod(x = "air_temperature_c",
    #                      y = "residuals",
    #                      by = c("water_sn", "baro_sn")),
    # 
    # water_temp_mods =
    #   correction_residuals(mods = delta_air_temp_mods, 
    #                        original.data = combined_data[experiment != "test-dat"],
    #                        x = "delta_at_c_min", 
    #                        y = "raw_error_cm", 
    #                        by = c("water_sn", "baro_sn")) %>% 
    #   correction_residuals(mods = delta_water_temp_mods, 
    #                        original.data = .,
    #                        x = "delta_wt_c_min", 
    #                        y = "residuals", 
    #                        by = c("water_sn", "baro_sn")) %>% 
    #   correction_residuals(mods = air_temp_mods, 
    #                        original.data = .,
    #                        x = "air_temperature_c", 
    #                        y = "residuals", 
    #                        by = c("water_sn", "baro_sn")) %>% 
    #   fit_correction_mod(x = "water_temperature_c",
    #                      y = "residuals",
    #                      by = c("water_sn", "baro_sn")),
    # 
    # compensated_data =
    #   combined_data %>% 
    #   .[, rect_water_level_cm := raw_water_level_cm] %>% 
    #   .[delta_air_temp_mods,
    #     rect_water_level_cm := rect_water_level_cm - slope * (delta_at_c_min - x_intercept),
    #     on = c("water_sn", "baro_sn")] %>% 
    #   .[delta_water_temp_mods,
    #     rect_water_level_cm := rect_water_level_cm - slope * (delta_wt_c_min - x_intercept),
    #     on = c("water_sn", "baro_sn")] %>% 
    #   .[air_temp_mods, 
    #     rect_water_level_cm := rect_water_level_cm - slope * (air_temperature_c - x_intercept),
    #     on = c("water_sn", "baro_sn")] %>% 
    #   .[water_temp_mods, 
    #     rect_water_level_cm := rect_water_level_cm - slope * (water_temperature_c - x_intercept),
    #     on = c("water_sn", "baro_sn")] %>% 
    #   .[, `:=`(raw_water_level_cm = raw_water_level_cm - median(raw_water_level_cm - water_depth_cm),
    #            rect_water_level_cm = rect_water_level_cm - median(rect_water_level_cm - water_depth_cm)),
    #     by = .(water_sn, baro_sn, experiment)] %>% 
    #   .[, `:=`(raw_abs_error_cm = abs(raw_water_level_cm - water_depth_cm),
    #            rect_abs_error_cm = abs(rect_water_level_cm - water_depth_cm))]
  )
