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
      read_xle_data(levellogger_measurements, "water"),
    
    raw_barometric_data = 
      read_xle_data(levellogger_measurements, "air") %>% 
      .[!between(sample_time, 
                 as.POSIXct("2020-06-10 06:00", tz = "EST5EDT"), 
                 as.POSIXct("2020-06-10 12:30", tz = "EST5EDT"))] %>%
      .[, `:=`(delta_t_c = c(NA_real_, diff(air_temperature_c)),
               delta_p_cm = c(NA_real_, diff(air_pressure_cm)),
               rolling_temp_c = frollmean(air_temperature_c, n = 720, align = "center"))],
    
    instrument_metadata = 
      get_metadata(raw_barometric_data, raw_water_data),
    
    mesonet_response = 
      fetch_external_baro(times = c(raw_barometric_data$sample_time)),
    
    combined_data = {
      combined <- 
        merge(raw_water_data[, .(water_sn = serial_number, sample_time, raw_water_pressure_cm = water_pressure_cm, dens_water_pressure_cm = 1000 * water_pressure_cm / water_density(water_temperature_c), water_temperature_c, water_error_cm = pressure_error_cm, water_temp_error_c = temperature_error_c)],
              raw_barometric_data[, .(baro_sn = serial_number, sample_time, raw_air_pressure_cm = air_pressure_cm, dens_air_pressure_cm = 1000 * air_pressure_cm / water_density(air_temperature_c), air_temperature_c, air_error_cm = pressure_error_cm, air_temp_error_c = temperature_error_c)],
              allow.cartesian = TRUE,
              by = "sample_time")
      
      combined[levellogger_measurements, 
               water_depth_cm := cable_length_cm - water_height_cm, 
               on = c(water_sn = "serial_number")]
      
      combined[, `:=`(delta_wt_01c_min = 100 * c(NA_real_, diff(water_temperature_c)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
                      delta_wp_cm_min = c(NA_real_, diff(raw_water_pressure_cm)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
                      delta_at_01c_min = 100 * c(NA_real_, diff(air_temperature_c)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
                      delta_at_error_01c_min = 100 * combine_errors(air_temp_error_c, air_temp_error_c),
                      delta_ap_cm_min = c(NA_real_, diff(raw_air_pressure_cm)) / c(NA_real_, as.numeric(diff(sample_time, unit = "mins"))),
                      instrument_error_cm = combine_errors(water_error_cm, air_error_cm),
                      raw_water_level_cm = dens_water_pressure_cm - dens_air_pressure_cm,
                      raw_error_cm = (dens_water_pressure_cm - dens_air_pressure_cm) - water_depth_cm),
               by = .(water_sn, baro_sn)]
        
      combined <- 
        set_experiments(combined, "data/experimental_design.csv")
      
      split(combined, by = "dataset")},
    
    lambdas =
      find_optimal_lambda(combined_data[["training"]], 
                          by = c("water_sn", "baro_sn"),
                          n = 100),
    
    # Should break this step up to save the models and then save the fit info
    # list as a separate target
    
    # Add samples to fit_metrics. Move lambda to fit_metrics
    bootstrap_models = {
      set.seed(1234)
      
      dat <- 
        split(Reduce(rbind, combined_data),
              by = c("water_sn", "baro_sn", "experiment"))
      
      mods <- mclapply(dat, 
                     function(x){
                       
                       i <- 1L
                       
                       reps <- 
                         1000L
                       
                       # models <- 
                       #   vector("list", 2)
                       
                       wsn <- 
                         unique(x$water_sn)
                       
                       bsn <- 
                         unique(x$baro_sn)
                       
                       model_id <- 
                         paste(wsn, bsn, unique(x$experiment), sep = "_")
                       
                       nobs <- 
                         nrow(x)
                       
                       models <- 
                         data.table(model_id, # = rep(model_id, reps), 
                                    rep = 1:reps,
                                    nobs = nobs,
                                    intercept = NA_real_,
                                    air_temperature_c = NA_real_,
                                    water_temperature_c = NA_real_,
                                    delta_at_01c_min = NA_real_, 
                                    sigma = NA_real_,
                                    adj_r2 = NA_real_,
                                    rmse = NA_real_)
                       
                       # fit_cols <-
                       #   names(models[, intercept:rmse])
                       
                       while(i < reps + 1){
                         
                         samples <- 
                           sample.int(n = nobs, size = nobs, replace = TRUE)
                         
                         inst_err <- 
                           rnorm(nobs, 0, x$instrument_error_cm)
                         
                         at_err <- 
                           rnorm(nobs, 0, x$air_error_c)
                         
                         wt_err <- 
                           rnorm(nobs, 0, x$water_error_c)
                         
                         d_at_err <- 
                           rnorm(nobs, 0, x$delta_at_error_01c_min)
                         
                         Y <-
                           matrix(x$raw_error_cm[samples] + inst_err,
                                  ncol = 1,
                                  dimnames = list(NULL, "raw_error_cm"))
                         
                         X <-
                           matrix(c(rep(1, nobs),
                                    x$air_temperature_c[samples] + at_err,
                                    x$water_temperature_c[samples] + wt_err,
                                    x$delta_at_01c_min[samples] + d_at_err),
                                  nrow = nobs,
                                  dimnames = list(NULL,
                                                  c("intercept",
                                                    "air_temperature_c",
                                                    "water_temperature_c",
                                                    "delta_at_01c_min")))
                         
                         p_B <- rep(999, 4)
                         
                         while(!all(p_B <= 0.05)){
                         
                           # Remove non-significant variables
                           X <- 
                             X[, which(p_B <= 0.05 | p_B == 999), drop = FALSE]
                           
                           # https://web.stanford.edu/~mrosenfe/soc_meth_proj3/matrix_OLS_NYU_notes.pdf
                           XX <- 
                             solve(crossprod(X))
                           
                           B <- 
                             XX %*% t(X) %*% Y
                           
                           fits <- 
                             X %*% B
                           
                           k <- 
                             nrow(B) - 1
                           
                           adj_r2 <- 
                             as.numeric(adjusted_r2(fits, Y, nobs, k))
                           
                           rs <- 
                             (fits - Y)^2
                           
                           rss <- 
                             sum(rs)
                           
                           rmse <- 
                             sqrt(mean(rs))
                           
                           df <- 
                             (nobs - k - 1)
                           
                           sigma <- 
                             sqrt(rss / df)
                           
                           # https://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression
                           se_B <- 
                             sqrt(diag(sigma^2 * XX))
                           
                           p_B <- 
                             c(intercept = 0, (2*(1 - pt(abs(B / se_B), df)))[-c(1), ])
                           
                           # B <- 
                           #   B * as.numeric(p_B <= 0.05)
                           
                         }
                         
                         write_fst(data.table(model_id, rep = i, samples, inst_err, at_err, wt_err, d_at_err),
                                   paste0("output/models/bootstrap_samples/", model_id, "_", i, ".fst"))
                         
                         model_fit <- 
                           setNames(c(t(B), sigma, adj_r2, rmse),
                                    c(rownames(B), "sigma", "adj_r2", "rmse"))
                         
                         for(I in names(model_fit)){
                           set(models, i, I, model_fit[[I]])
                         }
                         # 
                         # models[i, ] <- 
                         #   c(model_id, i, model_fit)
                         # 
                         # models[[i]] <- 
                         #   list(coefficients = data.table(model_id, 
                         #                                  rep = i, 
                         #                                  t(B), 
                         #                                  sigma,
                         #                                  adj_r2,
                         #                                  rmse))
                         
                         i <- i + 1L
                       }
                       models[, `:=`(air_temperature_c = nafill(air_temperature_c, "const", 0),
                                     water_temperature_c = nafill(water_temperature_c, "const", 0),
                                     delta_at_01c_min = nafill(delta_at_01c_min, "const", 0))]
                     }
                     
      )
      
      rbindlist(mods)[, c("water_sn", "baro_sn", "experiment") := tstrsplit(model_id, split = "_")][]
    }
    
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
