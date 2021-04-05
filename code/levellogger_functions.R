# High Level Functions ----------------------------------------------------

read_manual_measurements <- 
  function(file){
    fread(file,
          colClasses = c("character", "character", "integer", "numeric", "numeric"))
  }

read_xle_data <- 
  function(levellogger_info, type){
    
    stopifnot(length(type) == 1)
    
    files <- 
      make_xle_files(levellogger_info, type = type)
    
    dat <- 
      lapply(files, ingest_xle)
    
    dat <- 
      lapply(dat, as.data.table)
    
    dat <- 
      lapply(dat, 
             setnames, 
             old = c("temperature_deg_c", "temperature_c"), 
             new = paste0(c(type, type), "_temperature_c"),
             skip_absent = TRUE)
    
    dat <- 
      rbindlist(dat)
    
    dat[, sample_millisecond := NULL]
    
    dat[level_cm < 950, level_cm := level_cm + 950]
    
    setnames(dat, 
             old = "level_cm", 
             new = paste0(type, "_pressure_cm"))
    
    meta_dat <- 
      lapply(files, ingest_header)
    
    meta_dat <- 
      lapply(meta_dat, as.data.table)
    
    meta_dat <- 
      rbindlist(meta_dat, fill = TRUE)
    
    errors <- 
      data.table(instrument_type = c("LT_EDGE", "LT_EDGE_JR", "LT_Jr"),
                 pressure_error_cm = c(0.1, 0.5, 0.5) / qnorm(0.99),
                 temperature_error_c = c(0.05, 0.1, 0.1) / qnorm(0.99))
    
    meta_dat <- 
      meta_dat[errors, 
               nomatch = NULL,
               on = "instrument_type"]
    
    dat[, sample_time := as.POSIXct(paste(sample_date, sample_time),
                                    tz = "EST5EDT")]
    
    meta_dat[, `:=`(logger_start = as.POSIXct(logger_start, tz = "EST5EDT"),
                    logger_stop = as.POSIXct(logger_stop, tz = "EST5EDT"),
                    download_time = as.POSIXct(download_time, tz = "EST5EDT"))]
    
    dat[meta_dat,
        on = "input_source"]
  }

get_metadata <- 
  function(...){
    dat <- 
      rbindlist(lapply(list(...), 
                       function(x){
                         unique(x[, .(serial_number, 
                                      instrument_type, 
                                      model_number, 
                                      firmware, 
                                      software_version, 
                                      logger_start, 
                                      logger_stop, 
                                      n_records, 
                                      pressure_error_cm = pressure_error_cm, 
                                      temperature_error_c = temperature_error_c)])}))
    
    resolutions <- 
      data.table(instrument_type = c("LT_EDGE", "LT_EDGE_JR", "LT_Jr"),
                 pressure_resolution_cm = c(0.00012, 0.14, 0.14),
                 temperature_resolution_cm = c(0.003, 0.1, 0.1))
    
    dat[resolutions,
        on = "instrument_type"]
  }

combine_datasets <- 
  function(water.data, baro.data, 
           manual.measurements, design.file){
    combined <-
      merge(water.data[, .(water_sn = serial_number, sample_time, raw_water_pressure_cm = water_pressure_cm, dens_water_pressure_cm = 1000 * water_pressure_cm / water_density(water_temperature_c), water_temperature_c, water_error_cm = pressure_error_cm, water_temp_error_c = temperature_error_c)],
            baro.data[, .(baro_sn = serial_number, sample_time, raw_air_pressure_cm = air_pressure_cm, dens_air_pressure_cm = 1000 * air_pressure_cm / water_density(air_temperature_c), air_temperature_c, air_error_cm = pressure_error_cm, air_temp_error_c = temperature_error_c)],
            allow.cartesian = TRUE,
            by = "sample_time")
    
    combined[manual.measurements,
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
      set_experiments(combined, design.file)
    
    split(combined, by = "dataset")
  }

fit_models <- 
  function(data){
    set.seed(1234)
    
    dat <- 
      split(Reduce(rbind, data),
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
                                    model_rep = 1:reps,
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
                         
                         # fwrite(x = data.table(model_id, model_rep = i, samples, inst_err, at_err, wt_err, d_at_err),
                         #        file = "output/tabular/bootstrap_samples_and_errors.csv",
                         #        append = TRUE)
                          
                         write_fst(x = data.table(model_id, model_rep = i, samples, inst_err, at_err, wt_err, d_at_err),
                                   path = paste0("output/models/bootstrap_samples.df/", model_id, "_", i, ".fst"),
                                   compress = 100)
                         
                         model_fit <- 
                           setNames(c(t(B), sigma, adj_r2, rmse),
                                    c(rownames(B), "sigma", "adj_r2", "rmse"))
                         
                         for(I in names(model_fit)){
                           set(models, i, I, model_fit[[I]])
                         }

                         i <- i + 1L
                       }
                       models[, `:=`(air_temperature_c = nafill(air_temperature_c, "const", 0),
                                     water_temperature_c = nafill(water_temperature_c, "const", 0),
                                     delta_at_01c_min = nafill(delta_at_01c_min, "const", 0))]
                     }
                     
    )
    
    rbindlist(mods)[, c("water_sn", "baro_sn", "experiment") := tstrsplit(model_id, split = "_")][]
  }

calculate_fitted_values <- 
  function(models, data) {
  
    cat("\nCalculating Fitted Values\n")
  
    data <- 
      Reduce(rbind, data)
    
    pred_mat <-
      models[, .(B = list(matrix(c(intercept,
                                   air_temperature_c,
                                   water_temperature_c,
                                   delta_at_01c_min),
                                 dimnames = list(c("i",
                                                   "b_at",
                                                   "b_wt",
                                                   "b_dat"),
                                                 NULL),
                                 nrow = 4)),
                 sigma = list(matrix(rep(sigma, nobs),
                                     dimnames = list(NULL, "sigma"),
                                     ncol = 1))),
             keyby = .(water_sn, baro_sn, experiment, model_rep)]
    
    x_mat <-
      data[,  .(S = list(matrix(sample_time,
                                dimnames = list(NULL, "sample_time"),
                                ncol = 1)),
                X = list(matrix(c(rep(1, .N),
                                  air_temperature_c,
                                  water_temperature_c,
                                  delta_at_01c_min),
                                dimnames = list(NULL, c("intercept",
                                                        "air_temperature_c",
                                                        "water_temperature_c",
                                                        "delta_at_01c_min")),
                                ncol = 4))),
           keyby = .(water_sn, baro_sn, experiment)]
    
    predictions <-
      data[, .(water_sn,
               baro_sn,
               experiment,
               sample_time,
               predicted_error_cm = NA_real_,
               lower_bound_cm = NA_real_,
               upper_bound_cm = NA_real_)]
    
    setkey(predictions, "water_sn", "baro_sn", "experiment", "sample_time")
    
    # Keep the for loops. Unnesting all of these predictions very quickly runs
    # out of memory. It would probably be possible to fix it with careful joining
    # and unnesting, but this works
    
    pb <-
      txtProgressBar(min = 0, max = 5 * 36)
    
    c <- 0
    
    for(Ex in unique(data$experiment)){
      
      dat <- 
        pred_mat[experiment == Ex]
      
      pred_mat <- 
        pred_mat[experiment != Ex]
      
      gc()
      
      dat[x_mat,
          `:=`(S = i.S,
               X = i.X)]
      
      dat[, E_hat := Map(function(x, b, s){matrix(rnorm(n = nrow(x),
                                                        mean = x %*% b,
                                                        sd = s),
                                                  ncol = 1,
                                                  dimnames = list(NULL,
                                                                  "predicted_error_cm"))},
                         x = X,
                         b = B,
                         s = sigma)]
      
      # saveRDS(dat, fit_matrix_files()[Ex])
      
      setkey(dat, "water_sn", "baro_sn", "experiment")
      
      for(i in unique(dat$water_sn)){
        for(j in unique(dat$baro_sn)){
          post_pred <-
            dat[CJ(i,j),
                .(sample_time = as.POSIXct(S[[1]][, 1], tz = "EST5EDT", origin = "1970-01-01"), 
                  predicted_error_cm = E_hat[[1]][, 1]),
                by = .(water_sn, baro_sn, experiment, model_rep)]
          
          write_fst(post_pred,
                    path = paste0("output/tabular/bootstrap_fit_values.df/",
                                  i, "_", j, "_", Ex, ".fst"))
          
          ints <-
            post_pred[, .(pred_lwr = quantile(predicted_error_cm, probs = 0.025),
                          fit = quantile(predicted_error_cm, probs = 0.5),
                          pred_upr = quantile(predicted_error_cm, probs = 0.975)),
                      keyby = .(water_sn, baro_sn, experiment, sample_time)]
          
          predictions[ints,
                      `:=`(predicted_error_cm = fit,
                           lower_bound_cm = pred_lwr,
                           upper_bound_cm = pred_upr)]
          
          c <- c + 1/36
          
          setTxtProgressBar(pb, c)
        }
      }
      
      setkey(data, "water_sn", "baro_sn", "experiment", "sample_time")
      setkey(predictions, "water_sn", "baro_sn", "experiment", "sample_time")
      
      predictions[data,
                  `:=`(water_depth_cm = i.water_depth_cm,
                       raw_water_level_cm = i.raw_water_level_cm,
                       raw_error_cm = i.raw_error_cm,
                       instrument_error_cm = i.instrument_error_cm)]
      
      predictions[, rect_water_level_cm := raw_water_level_cm - predicted_error_cm]
      predictions[, `:=`(centered_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
                         rect_water_level_cm = rect_water_level_cm - mean(rect_water_level_cm - water_depth_cm)),
                  by = .(water_sn, baro_sn, experiment)]
      
      # Weird problem where between() (data.table or dplyr) was returning TRUE when FALSE
      predictions[, `:=`(instrument_lower = water_depth_cm - qnorm(0.975) * instrument_error_cm,
                         instrument_upper = water_depth_cm + qnorm(0.975) * instrument_error_cm,
                         propagated_lower = water_depth_cm - (predicted_error_cm - lower_bound_cm),
                         propagated_upper = water_depth_cm + (upper_bound_cm - predicted_error_cm))][]
    }
    
    predictions
  }

calculate_predicted_values <- 
  function(models, data) {
    
    cat("\nCalculating Predicted Values\n")
    
    data <- 
      data$testing
    
    pred_mat <-
      models[experiment == "test-dat", 
             .(B = list(matrix(c(intercept,
                                 air_temperature_c,
                                 water_temperature_c,
                                 delta_at_01c_min),
                               dimnames = list(c("i",
                                                 "b_at",
                                                 "b_wt",
                                                 "b_dat"),
                                               NULL),
                               nrow = 4)),
               sigma),
             by = .(water_sn, baro_sn, model_rep)]
    
    setkey(pred_mat, water_sn, baro_sn)
    
    x_mat <-
      data[,  .(S = list(matrix(sample_time,
                                dimnames = list(NULL, "sample_time"),
                                ncol = 1)),
                X = list(matrix(c(rep(1, .N),
                                  air_temperature_c,
                                  water_temperature_c,
                                  delta_at_01c_min),
                                dimnames = list(NULL, c("intercept",
                                                        "air_temperature_c",
                                                        "water_temperature_c",
                                                        "delta_at_01c_min")),
                                ncol = 4))),
           keyby = .(water_sn, baro_sn, experiment)]
    
    predictions <-
      data[, .(water_sn,
               baro_sn,
               experiment,
               sample_time,
               predicted_error_cm = NA_real_,
               lower_bound_cm = NA_real_,
               upper_bound_cm = NA_real_)]
    
    setkey(predictions, "water_sn", "baro_sn", "experiment", "sample_time")
    
    # Keep the for loops. Unnesting all of these predictions very quickly runs
    # out of memory. It would probably be possible to fix it with careful joining
    # and unnesting, but this works
    
    pb <-
      txtProgressBar(min = 0, max = 36)
    
    c <- 0
    
    for(Ex in unique(data$experiment)){
      
      dat <- 
        x_mat[experiment == Ex]
      
      x_mat <- 
        x_mat[experiment != Ex]
      
      gc()
      
      dat <- 
        dat[pred_mat]
      
      dat[, sigma := Map(function(s, x){matrix(rep(s, nrow(x)), 
                                               dimnames = list(NULL, "sigma"), 
                                               ncol = 1)}, 
                         s = sigma, 
                         x = X)]

      dat[, E_hat := Map(function(x, b, s){matrix(rnorm(n = nrow(x),
                                                        mean = x %*% b,
                                                        sd = s),
                                                  ncol = 1,
                                                  dimnames = list(NULL,
                                                                  "predicted_error_cm"))},
                         x = X,
                         b = B,
                         s = sigma)]
      
      setkey(dat, "water_sn", "baro_sn", "experiment")
      
      for(i in unique(dat$water_sn)){
        for(j in unique(dat$baro_sn)){
          post_pred <-
            dat[CJ(i,j),
                .(sample_time = as.POSIXct(S[[1]][, 1], tz = "EST5EDT", origin = "1970-01-01"), 
                  predicted_error_cm = E_hat[[1]][, 1]),
                by = .(water_sn, baro_sn, experiment, model_rep)]
          
          write_fst(post_pred,
                    path = paste0("output/tabular/bootstrap_predict_values.df/",
                                  i, "_", j, "_", Ex, ".fst"))
          
          ints <-
            post_pred[, .(pred_lwr = quantile(predicted_error_cm, probs = 0.025),
                          fit = quantile(predicted_error_cm, probs = 0.5),
                          pred_upr = quantile(predicted_error_cm, probs = 0.975)),
                      keyby = .(water_sn, baro_sn, experiment, sample_time)]
          
          predictions[ints,
                      `:=`(predicted_error_cm = fit,
                           lower_bound_cm = pred_lwr,
                           upper_bound_cm = pred_upr)]
          
          c <- c + 1/36
          
          setTxtProgressBar(pb, c)
        }
      }
      
      setkey(data, "water_sn", "baro_sn", "experiment", "sample_time")
      setkey(predictions, "water_sn", "baro_sn", "experiment", "sample_time")
      
      predictions[data,
                  `:=`(water_depth_cm = i.water_depth_cm,
                       raw_water_level_cm = i.raw_water_level_cm,
                       raw_error_cm = i.raw_error_cm,
                       instrument_error_cm = i.instrument_error_cm)]
      
      predictions[, rect_water_level_cm := raw_water_level_cm - predicted_error_cm]
      predictions[, `:=`(centered_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
                         rect_water_level_cm = rect_water_level_cm - mean(rect_water_level_cm - water_depth_cm)),
                  by = .(water_sn, baro_sn, experiment)]
      
      # Weird problem where between() (data.table or dplyr) was returning TRUE when FALSE
      predictions[, `:=`(instrument_lower = water_depth_cm - qnorm(0.975) * instrument_error_cm,
                         instrument_upper = water_depth_cm + qnorm(0.975) * instrument_error_cm,
                         propagated_lower = water_depth_cm - (predicted_error_cm - lower_bound_cm),
                         propagated_upper = water_depth_cm + (upper_bound_cm - predicted_error_cm))][]
    }
    
    predictions
  }

export_correction_models <- 
  function(models, instrument.metadata, file.out, generic.ltjr.file = NULL){
    
    median_mods <- 
      models[experiment == "test-dat", 
             lapply(.SD, median), 
             by = .(water_sn, baro_sn),
             .SDcols = c("intercept", "air_temperature_c", 
                         "water_temperature_c", "delta_at_01c_min")]
    
    median_mods[, instrument_error_cm := combine_errors(instrument.metadata[serial_number == .BY[[1]], pressure_error_cm], 
                                                     instrument.metadata[serial_number == .BY[[2]], pressure_error_cm]), 
                by = .(water_sn, baro_sn)]
    
    fwrite(median_mods, 
           file.out)
    
    if(!is.null(generic.ltjr.file)){
      generic_mods <- 
        models[instrument.metadata[instrument_type == "LT_Jr",
                               .(serial_number)],
               on = c("water_sn" = "serial_number"), 
               nomatch = NULL]
      
      generic_mods <- 
        generic_mods[experiment == "test-dat", 
                     lapply(.SD, median),
                     by = .(baro_sn), 
                     .SDcols = c("intercept", "air_temperature_c", 
                                 "water_temperature_c", "delta_at_01c_min")]
      
      fwrite(generic_mods, 
             generic.ltjr.file)
    }
    
    median_mods
  }

# Figures -----------------------------------------------------------------


create_drivers_panel <- 
  function(data){
    
    dat_1 <- 
      data$testing[experiment == "var-dis" & baro_sn == "1066019" & water_sn == "1062452"]
    
    fig1a <- 
      ggplot(dat_1,
             aes(x = air_temperature_c,
                 y = raw_error_cm)) + 
      geom_point(color = "#329985") +
      geom_smooth(method = "lm",
                  formula = "y~x",
                  se = FALSE,
                  color = "black") +
      labs(x = expression(paste("Air Temperature, ", degree, "C")),
           y = "Raw Error, cm")
    
    fig1b <- 
      ggplot(dat_1,
             aes(x = water_temperature_c,
                 y = raw_error_cm,
                 color = baro_sn)) + 
      geom_point(color = "#329985") +
      geom_smooth(method = "lm",
                  formula = "y~x",
                  se = FALSE,
                  color = "black") +
      labs(x = expression(paste("Water Temperature, ", degree, "C")),
           y = "Raw Error, cm")
    
    fig1c <- 
      ggplot(dat_1,
             aes(x = water_temperature_c,
                 y = raw_error_cm - fitted(lm(raw_error_cm ~ air_temperature_c + delta_at_01c_min, data = dat_1)),
                 color = baro_sn)) + 
      geom_point(color = "#329985") +
      geom_smooth(method = "lm",
                  formula = "y~x",
                  se = FALSE,
                  color = "black") +
      labs(x = expression(paste("Water Temperature, ", degree, "C")),
           y = "Raw Error minus Air Temperature Effect, cm") +
      theme_few()
    
    {fig1a + fig1b + fig1c + plot_annotation(tag_levels = "A")}
  }

create_bootstrap_timeseries <- 
  function(fits, models, exp, water.sn, baro.sn){
    dat2 <- 
      fits[experiment == exp & baro_sn == baro.sn & water_sn == water.sn]
    
    model_file <- 
      glue("output/tabular/bootstrap_fit_values.df/{water.sn}_{baro.sn}_{exp}.fst")
    
    mods <- 
      read_fst(model_file,
               as.data.table = TRUE)
    
    r2_quartiles <- 
      quantile(models[water_sn == water.sn &
                                  baro_sn == baro.sn & 
                                  experiment == exp, 
                                adj_r2], 
               seq(0, 1, 0.5),
               type = 1)
    
    quartile_reps <- 
      models[water_sn == water.sn &
                         baro_sn == baro.sn & 
                         experiment == exp & 
                         adj_r2 %in% r2_quartiles, 
                       model_rep]
    
    mods2 <- 
      mods[model_rep %in% quartile_reps]
    
    mods2[dat2, 
          `:=`(raw_water_level_cm = i.raw_water_level_cm,
               centered_water_level_cm = i.centered_water_level_cm,
               orig_rect_water_level_cm = i.rect_water_level_cm,
               water_depth_cm = i.water_depth_cm,
               instrument_lower = i.instrument_lower,
               instrument_upper = i.instrument_upper,
               propagated_lower = i.propagated_lower,
               propagated_upper = i.propagated_upper),
          on = c("sample_time")]
    
    mods2[, rect_water_level_cm := (raw_water_level_cm - predicted_error_cm) - mean(raw_water_level_cm - predicted_error_cm - water_depth_cm),
          by = .(model_rep)]
    
    mods2[, sample := "no"]
    mods2[between(sample_time, 
                  as.POSIXct("2020-06-10 18:00", tz = "EST5EDT"), 
                  as.POSIXct("2020-06-10 20:00", tz = "EST5EDT")),
          sample := "yes"]
    
    ggplot(mods2, 
           aes(x = sample_time)) +
      geom_ribbon(aes(ymin = propagated_lower,
                      ymax = propagated_upper),
                  fill = "gray85") +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  fill = "gray75") +
      geom_line(aes(y = rect_water_level_cm,
                    color = as.factor(model_rep)),
                show.legend = FALSE) +
      geom_line(aes(y = orig_rect_water_level_cm),
                color = "gray5",
                size = 1) +
      geom_line(aes(y = centered_water_level_cm),
                color = "gray5",
                linetype = "dotted") +
      coord_cartesian(expand = FALSE) +
      labs(y = "Water Level", 
           x = "Sample Time") +
      facet_zoom(x = sample == "yes",
                 zoom.size = 0.5, 
                 horizontal = FALSE,
                 ylim = range(mods2[sample == "yes", rect_water_level_cm])) +
      theme_few()+ 
      theme(strip.background = element_rect(color = "black",
                                            size = rel(1.5)))
  }

create_coefficients_panel <- 
  function(models){
    
    dat4a <- 
      models[experiment == "test-dat", 
                       .(y = mean(air_temperature_c),
                         ymin = quantile(air_temperature_c, 0.025), 
                         ymax = quantile(air_temperature_c, 0.975)), 
                       by = .(water_sn, baro_sn)] 
    
    dat4a[, outlier_sn := ifelse(water_sn == "2025928",
                                 "yes",
                                 "no")]  
    
    fig4a <- 
      ggplot(dat4a,
             aes(y = y, ymin = ymin, ymax = ymax,
                 x = baro_sn,
                 group = water_sn,
                 color = outlier_sn,
                 shape = outlier_sn)) +
      geom_pointrange(position = position_dodge(width = 0.8),
                      show.legend = FALSE,
                      fill = "white") +
      labs(y = expression(hat(beta)[Air~Temperature]),
           x = "Barometric Transducer Serial Number") +
      scale_shape_manual(values = c(19, 21)) +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = 1))
    
    dat4b <- 
      models[experiment == "test-dat", 
                       .(y = mean(water_temperature_c),
                         ymin = quantile(water_temperature_c, 0.025), 
                         ymax = quantile(water_temperature_c, 0.975)), 
                       by = .(water_sn, baro_sn)] 
    
    dat4b[, outlier_sn := ifelse(water_sn == "2025928",
                                 "yes",
                                 "no")]  
    
    fig4b <- 
      ggplot(dat4b,
             aes(y = y, ymin = ymin, ymax = ymax,
                 x = water_sn,
                 group = baro_sn,
                 color = outlier_sn,
                 shape = outlier_sn)) +
      geom_pointrange(position = position_dodge(width = 0.8),
                      show.legend = FALSE,
                      fill = "white") +
      labs(y = expression(hat(beta)[Water~Temperature]),
           x = "Water Transducer Serial Number") +
      scale_shape_manual(values = c(19, 21)) +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = 1))
    
    dat4c <- 
      models[experiment == "test-dat", 
             .(y = mean(delta_at_01c_min),
               ymin = quantile(delta_at_01c_min, 0.025), 
               ymax = quantile(delta_at_01c_min, 0.975)), 
             by = .(water_sn, baro_sn)]
    
    dat4c[, outlier_sn := ifelse(water_sn == "2025928",
                                 "yes",
                                 "no")]
    
    fig4c <- 
      ggplot(dat4c,
             aes(y = y, ymin = ymin, ymax = ymax,
                 x = baro_sn,
                 group = water_sn,
                 color = outlier_sn,
                 shape = outlier_sn)) +
      geom_pointrange(position = position_dodge(width = 0.8),
                      show.legend = FALSE,
                      fill = "white") +
      labs(y = expression(hat(beta)[Delta~Air~Temperature]),
           x = "Barometric Transducer Serial Number") +
      scale_shape_manual(values = c(19, 21)) +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = 1))
    
    fig4a + fig4b + fig4c + plot_annotation(tag_levels = "A")
  }

create_case_study_panel <- 
  function(data, out.path){
    
    panel_data <- 
      data[between(sample_time, 
                         as.POSIXct("2018-08-16 00:00:00", tz = "EST5EDT"), 
                         as.POSIXct("2018-08-17 23:45:00", tz = "EST5EDT")),
                 .(sample_time,
                   sample_date = as.Date(sample_time, tz = "EST5EDT"),
                   water_temperature_c,
                   air_temperature_c,
                   raw_compensated_level_cm,
                   corrected_compensated_level_cm)]
    
    panel_data[, 
               `:=`(raw_white_cm = predict(lm(raw_compensated_level_cm ~ sample_time,
                                              data = .SD[hour(sample_time) <= 7]),
                                           newdata = .SD),
                    raw_white_slope = coef(lm(raw_compensated_level_cm ~ sample_time,
                                              data = .SD[hour(sample_time) <= 7]))[2],
                    corrected_white_cm = predict(lm(corrected_compensated_level_cm ~ sample_time,
                                                    data = .SD[hour(sample_time) <= 7]),
                                                 newdata = .SD),
                    corrected_white_slope = coef(lm(corrected_compensated_level_cm ~ sample_time,
                                                    data = .SD[hour(sample_time) <= 7]))[2]),
               by = .(sample_date)]
    
    panel_g_labels <- 
      panel_data[sample_time %in% c(as.POSIXct("2018-08-17 12:00", tz = "EST5EDT"),
                                    as.POSIXct("2018-08-17 18:00", tz = "EST5EDT")),
                 .(sample_date,
                   sample_time,
                   label_y = c(raw_white_cm[1], corrected_white_cm[2]),
                   label = c("Uncorrected~G[`in`]", "Corrected~G[`in`]"))]
    
    panel_et_labels <- 
      rbind(panel_data[sample_time == as.POSIXct("2018-08-16 23:45", tz = "EST5EDT"),
                       .(sample_date,
                         sample_time = sample_time + 3600,
                         label_y = raw_compensated_level_cm + 0.75*(raw_white_cm - raw_compensated_level_cm),
                         label = "Uncorrected~ET")],
            panel_data[sample_time == as.POSIXct("2018-08-16 23:45", tz = "EST5EDT"),
                       .(sample_date,
                         sample_time = sample_time + 9000,
                         label_y = (corrected_white_cm + corrected_compensated_level_cm) / 2,
                         label = "Corrected~ET")])
    
    panel <- 
      {ggplot(data, aes(x = sample_time)) + 
          geom_rect(aes(xmin = min(panel_data$sample_time),
                        xmax = max(panel_data$sample_time),
                        ymin = min(panel_data$corrected_compensated_level_cm),
                        ymax = max(panel_data$corrected_compensated_level_cm)),
                    fill = 'gray80',
                    color = NA) +
          geom_line(aes(y = corrected_compensated_level_cm, 
                        linetype = 'Corrected',
                        color = 'Corrected')) +
          geom_line(aes(y = raw_compensated_level_cm, 
                        linetype = 'Uncorrected',
                        color = 'Uncorrected')) + 
          scale_color_manual(name = NULL, 
                             values = c(Corrected = 'black',
                                        Uncorrected = 'blue')) +
          scale_linetype_manual(name = NULL, 
                                values = c(Corrected = 'solid',
                                           Uncorrected = 'dashed')) +
          ylab("Water Level (cm)")} / 
      {ggplot(panel_data) +
          aes(x = sample_time) +
          geom_line(aes(y = corrected_compensated_level_cm)) + 
          geom_line(aes(y = corrected_white_cm),
                    color = 'black',
                    linetype = 'dotted') +
          geom_line(aes(y = raw_compensated_level_cm), 
                    color = 'blue',
                    linetype = 'dashed') +
          geom_line(aes(y = raw_white_cm),
                    color = 'blue',
                    linetype = 'dotted') +
          geom_segment(data = panel_data[, last(.SD), by = sample_date],
                       aes(y = raw_compensated_level_cm,
                           yend = raw_white_cm,
                           x = sample_time + 3600,
                           xend = sample_time + 3600),
                       linetype = 'longdash',
                       arrow = arrow(angle = 90, ends = 'both', length = unit(0.25, "lines")),
                       color = xaringanthemer::lighten_color('blue', strength = 0.6)) +
          geom_segment(data = panel_data[, last(.SD), by = sample_date],
                       aes(y = corrected_compensated_level_cm,
                           yend = corrected_white_cm,
                           x = sample_time + 9000,
                           xend = sample_time + 9000),
                       arrow = arrow(angle = 90, ends = 'both', length = unit(0.25, "lines")),
                       color = 'gray20') +
          geom_text(data = panel_g_labels,
                    aes(x = sample_time,
                        y = label_y,
                        label = label),
                    angle = c(9, -6),
                    vjust = 0,
                    size = 14*5/14,
                    parse = TRUE) +
          geom_text(data = panel_et_labels,
                    aes(x = sample_time,
                        y = label_y,
                        label = label),
                    angle = 90,
                    vjust = -0.15,
                    size = 14*5/14,
                    parse = TRUE) +
          facet_wrap(~sample_date,
                     scales = "free",
                     nrow = 1,
                     strip.position = 'bottom') + 
          ylab("Water Level (cm)") +
          scale_x_datetime(breaks = seq(as.POSIXct("2018-08-13 06:00:00", tz = "EST5EDT"), 
                                        as.POSIXct("2018-08-20 18:00:00", tz = "EST5EDT"), 
                                        by = 12*3600),
                           date_labels = "%H:%M", 
                           expand = expansion(mult = c(0, 0),
                                              add = c(0, 3600)))} +
      plot_annotation(tag_levels = "A") &
      theme_minimal(base_size = 20) +
      theme(strip.placement = 'outside',
            axis.title.x = element_blank(),
            legend.position = c(0.08, 0.22),
            legend.background = element_rect(fill = 'white',
                                             color = NA),
            legend.margin = margin(0, 0, 0, 0))
    
    ggsave(plot = panel,
           filename = out.path,
           width = 16,
           height = 8)
  }

# Case Study Functions ----------------------------------------------------

load_case_study <- 
  function(data) {
    dat <- 
      fread(data,
            colClasses = c("character", 
                           "character",
                           "POSIXct",
                           "numeric",
                           "numeric",
                           "numeric",
                           "numeric",
                           "numeric",
                           "numeric",
                           "Date",
                           "numeric"))
    
    dat[, sample_time := setattr(sample_time, "tzone", "EST5EDT")]
    
    setkey(dat, "sample_time")
    
    dat
  }

correct_data <- 
  function(data, models){
    
    data[models, 
         `:=`(error_cm = air_temperature_c * i.air_temperature_c + water_temperature_c * i.water_temperature_c + delta_at_01c_min * i.delta_at_01c_min,
              instrument_error_cm = i.instrument_error_cm),
         on = c("water_sn", "baro_sn")]
    
    data[, corrected_water_level_cm := raw_water_level_cm - error_cm]
    
    data
  }

compensate_data <- 
  function(data, calibration.data){
    
    c_dat <- 
      fread(calibration.data)
    
    attributes(c_dat)$tzone <- "EST5EDT"
    
    data[c_dat, 
         on = "sample_time"]
    
    c_dat <- 
      na.omit(data[c_dat, on = "sample_time"], 
              cols = c("raw_water_level_cm", "corrected_water_level_cm"))
    
    c_dat <- 
      c_dat[, .(raw_correction_cm = mean(raw_water_level_cm - (well_height_cm - manual_level_cm)),
                corrected_correction_cm = mean(corrected_water_level_cm - (well_height_cm - manual_level_cm))),
            by = .(field_season)]
    
    data[, field_season := year(sample_time)]
    
    data[c_dat,
         `:=`(raw_compensated_level_cm = raw_water_level_cm - i.raw_correction_cm,
              corrected_compensated_level_cm = corrected_water_level_cm - i.corrected_correction_cm),
         on = "field_season"]
    
    data[, field_season := NULL]
    
    data
  }

smooth_data <- 
  function(data, n, ...){
    
    for(col in c(...)){
      set(data, j = col, value = frollmean(data[[col]], n = n, align = "center"))
    }
    
    # data[, `:=`(raw_compensated_level_cm = frollmean(raw_compensated_level_cm, n = n, align = "center"),
    #             corrected_compensated_level_cm = frollmean(corrected_compensated_level_cm, n = n, align = "center"))]
    
    data
  }

calculate_sy <- 
  function(data){
    
    # ESy function & min.esy taken from climate_impacts
    esy_function <- 
      function (wl = NULL, min.esy = 1.00046) 
        pmax(min.esy, 9.86792148868664 - (9.86792148868664 - 2.39189118793206) * 
               exp(0.00983144846311064 * wl))
    
    data[, `:=`(corrected_sy = 1/esy_function(corrected_compensated_level_cm),
                raw_sy = 1/esy_function(raw_compensated_level_cm))]
   
    data
  }

subset_year <- 
  function(data, year){
    data[year(sample_time) == year]
  }

subset_case_study <- 
  function(data, start, end){
    data[between(sample_time, 
                 as.POSIXct(start, tz = "EST5EDT"), 
                 as.POSIXct(end, tz = "EST5EDT"))]
  }

calculate_detrended_g <- 
  function(data){
    
    recharge_period <- 
      0:7
    
    # nested <-
    #   map_dfr(unique(data$sample_date)[-c(length(unique(data$sample_date)))],
    #            # c(list(c(0, 1, 2)), replicate(length(unique(data$sample_date))-2, -1:1, simplify = FALSE)),
    #            ~data.table(sample_date = .x,
    #                        dat = list(data[hour(sample_time) %in% recharge_period & sample_date %in% c(.x + -1:1)])))
    # 
    # # nested[sample_date %in% as.Date(paste0("2018-08-", 10:23)), dat] %>%
    # #   map(~ggplot(data = .x,
    # #               aes(x = sample_time, y = corrected_compensated_level_cm)) +
    # #         ggtitle(unique(.x$sample_date)) +
    # #         geom_point() +
    # #         geom_smooth(method = "lm", se = FALSE)) %>%
    # #   reduce(`+`)
    # 
    # nested[, `:=`(raw_detrend_mod = map(dat, function(x){lm(raw_compensated_level_cm ~ sample_time, data = x)}),
    #               corrected_detrend_mod = map(dat, function(x){lm(corrected_compensated_level_cm ~ sample_time, data = x)}))]
    # 
    # nested[, `:=`(raw_detrend_m_cm_s = map_dbl(raw_detrend_mod,
    #                                    ~coef(.x)[["sample_time"]]),
    #               raw_detrend_b_cm_s = map_dbl(raw_detrend_mod,
    #                                    ~coef(.x)[["(Intercept)"]]),
    #               corrected_detrend_m_cm_s = map_dbl(corrected_detrend_mod,
    #                                          ~coef(.x)[["sample_time"]]),
    #               corrected_detrend_b_cm_s = map_dbl(corrected_detrend_mod,
    #                                          ~coef(.x)[["(Intercept)"]]))]
    # 
    # data[nested,
    #      `:=`(raw_detrend_m_cm_s = i.raw_detrend_m_cm_s,
    #           raw_detrend_b_cm_s = i.raw_detrend_b_cm_s,
    #           corrected_detrend_m_cm_s = i.corrected_detrend_m_cm_s,
    #           corrected_detrend_b_cm_s = i.corrected_detrend_b_cm_s),
    #      on = "sample_date"]
    # 
    # data[,`:=`(raw_detrended_cm = raw_compensated_level_cm - (raw_detrend_b_cm_s + raw_detrend_m_cm_s * as.numeric(sample_time)),
    #            corrected_detrended_cm = corrected_compensated_level_cm - (corrected_detrend_b_cm_s + corrected_detrend_m_cm_s * as.numeric(sample_time)))]

    # daily detrend:
    nested <-
      data[hour(sample_time) %in% recharge_period, 
           .(raw_i_cm = first(raw_compensated_level_cm),
             corrected_i_cm = first(corrected_compensated_level_cm)),
           by = .(sample_date)]

    nested[, `:=`(raw_f_cm = shift(raw_i_cm, -1),
                 corrected_f_cm = shift(corrected_i_cm, -1))]

    nested[, `:=`(raw_detrend_m_cm_s = (raw_f_cm - raw_i_cm) / 86400,
                  corrected_detrend_m_cm_s = (corrected_f_cm - corrected_i_cm) / 86400)]

    nested[, `:=`(raw_detrend_b_cm_s = raw_i_cm,
                  corrected_detrend_b_cm_s =corrected_i_cm)]

    data[nested,
         `:=`(raw_detrend_m_cm_s = i.raw_detrend_m_cm_s,
              raw_detrend_b_cm_s = i.raw_detrend_b_cm_s,
              corrected_detrend_m_cm_s = i.corrected_detrend_m_cm_s,
              corrected_detrend_b_cm_s = i.corrected_detrend_b_cm_s),
         on = "sample_date"]

    data[,`:=`(raw_detrended_cm = raw_compensated_level_cm - (raw_detrend_b_cm_s + raw_detrend_m_cm_s * as.numeric(sample_time - min(sample_time))),
               corrected_detrended_cm = corrected_compensated_level_cm - (corrected_detrend_b_cm_s + corrected_detrend_m_cm_s * as.numeric(sample_time - min(sample_time)))),
         by = .(sample_date)]

    data[!(hour(sample_time) %in% recharge_period),
         `:=`(raw_detrend_m_cm_s = NA_real_,
              raw_detrend_b_cm_s = NA_real_,
              corrected_detrend_m_cm_s = NA_real_,
              corrected_detrend_b_cm_s = NA_real_)]

    data[, `:=`(raw_detrend_m_cm_s = zoo::na.spline(raw_detrend_m_cm_s),
                raw_detrend_b_cm_s = zoo::na.spline(raw_detrend_b_cm_s), 
                corrected_detrend_m_cm_s = zoo::na.spline(corrected_detrend_m_cm_s), 
                corrected_detrend_b_cm_s = zoo::na.spline(corrected_detrend_b_cm_s))]

    data[, `:=`(raw_dwtDT_dt = three_point_slope(sample_time, raw_detrended_cm),
                corrected_dwtDT_dt = three_point_slope(sample_time, corrected_detrended_cm))]
    
    # Gamma_dat <- 
    #   data[, .(sample_time,
    #            sample_date,
    #            raw_compensated_level_cm,
    #            corrected_compensated_level_cm,
    #            raw_detrended_cm,
    #            corrected_detrended_cm,
    #            raw_dwtDT_dt,
    #            corrected_dwtDT_dt)]
    # 
    # Gamma_dat <- 
    #   Gamma_dat[hour(sample_time) %in% recharge_period]

    # Need to create Gamma(wt) for each day separately. Then need to find a 
    # way to smoothly move between the two days to create a function. If it were
    # just slope it would be easy to use linear interpretation between the two 
    # recharge period models. But it isn't, need to think about how to transition
    # slope and intercept
    # The assumption of an approximately linear relationship for Gamma(wt) 
    # between recharge periods does not hold
    
    # Gamma_dat <- 
    #   map_dfr(unique(Gamma_dat$sample_date)[-c(length(unique(Gamma_dat$sample_date)))], 
    #           ~data.table(sample_date = .x, 
    #                       data = list(Gamma_dat[sample_date %in% c(.x)])))
    # 
    # Gamma_dat[sample_date %in% as.Date(paste0("2018-08-", 10:23)), data] %>%
    #   map(~ggplot(data = .x,
    #               aes(x = corrected_detrended_cm, y = corrected_dwtDT_dt)) +
    #         ggtitle(unique(.x$sample_date)) +
    #         geom_point(aes(color = as.factor(sample_date)),
    #                    show.legend = FALSE) +
    #         geom_smooth(method = "lm", se = FALSE, formula = "y ~ x")) %>%
    #   reduce(`+`)
    # 
    # Gamma_dat[sample_date %in% as.Date(paste0("2018-08-", 10:23)), data] %>%
    #   map(~ggplot(data = .x,
    #               aes(x = sample_time, y = three_point_slope(sample_time, corrected_dwtDT_dt))) +
    #         ggtitle(unique(.x$sample_date)) +
    #         geom_point(aes(color = as.factor(sample_date)),
    #                    show.legend = FALSE) +
    #         geom_smooth(method = "lm", se = FALSE, formula = "y ~ x")) %>%
    #   reduce(`+`)
    
    # Gamma_dat[sample_date %in% as.Date(paste0("2018-08-", 10:23)), data][1:13] %>% 
    #   rbindlist() %>% 
    #   ggplot(aes(x = corrected_detrended_cm, 
    #              y = corrected_dwtDT_dt, 
    #              color = as.factor(sample_date))) +
    #   geom_point() + scale_color_viridis_d() +
    #   geom_smooth(method = "lm", se = FALSE)
    
    # For daily data need to remove days with no valid observations
    # null_dates <- 
    #   Gamma_dat[, .(sample_date, 
    #                 nobs = map_dbl(data, nrow), 
    #                 na_dwtDT_dt = map_dbl(data, ~sum(is.na(.x$corrected_dwtDT_dt))), 
    #                 na_WTDT = map_dbl(data, ~sum(is.na(.x$corrected_detrended_cm))))]
    # 
    # null_dates <- 
    #   null_dates[nobs == na_dwtDT_dt | nobs == na_WTDT, sample_date]
    #     
    # Gamma_wt <- 
    #   Gamma_dat[!(sample_date %in% null_dates), 
    #             c(as.list(coef(lm(corrected_dwtDT_dt ~ corrected_detrended_cm, data = data[[1]]))),
    #               as.list(coef(lm(raw_dwtDT_dt ~ raw_detrended_cm, data = data[[1]])))), 
    #             by = .(sample_date)]
    # 
    # setnames(Gamma_wt, c("sample_date", "corrected_b", "corrected_m", "raw_b", "raw_m"))
    # 
    # data[Gamma_wt,
    #      `:=`(raw_gamma_m = i.raw_m,
    #           raw_gamma_b = i.raw_b,
    #           corrected_gamma_m = i.corrected_m,
    #           corrected_gamma_b = i.corrected_b),
    #      on = "sample_date"]
    
    # data[!(hour(sample_time) %in% recharge_period),
    #      `:=`(raw_gamma_m = NA_real_,
    #           raw_gamma_b = NA_real_,
    #           corrected_gamma_m = NA_real_,
    #           corrected_gamma_b = NA_real_)]
    # 
    # data[, `:=`(raw_gamma_m = zoo::na.spline(raw_gamma_m),
    #             raw_gamma_b = zoo::na.spline(raw_gamma_b), 
    #             corrected_gamma_m = zoo::na.spline(corrected_gamma_m), 
    #             corrected_gamma_b = zoo::na.spline(corrected_gamma_b))]
    # 
    # data[,`:=`(raw_Gamma_cm_s = raw_gamma_b + raw_gamma_m * raw_detrended_cm,
    #            corrected_Gamma_cm_s = corrected_gamma_b + corrected_gamma_m * corrected_detrended_cm)]

    data[, `:=`(raw_gamma_m = three_point_slope(sample_time, raw_dwtDT_dt),
                corrected_gamma_m = corrected_dwtDT_dt + three_point_slope(sample_time, corrected_dwtDT_dt))]
    
    data[!(hour(sample_time) %in% recharge_period),
         `:=`(raw_Gamma_cm_s = NA_real_,
              corrected_Gamma_cm_s = NA_real_)]
    
    data[!is.na(raw_gamma_m) & hour(sample_time) %in% recharge_period, 
         raw_Gamma_cm_s := fitted(lm(raw_gamma_m ~ sample_time)),
         by = .(sample_date)]
    
    data[!is.na(corrected_gamma_m) & hour(sample_time) %in% recharge_period,
         corrected_Gamma_cm_s := fitted(lm(corrected_gamma_m ~ sample_time)),
         by = .(sample_date)]
    
    # data[, `:=`(raw_gamma_m = zoo::na.approx(raw_gamma_m, rule = 2),
    #             corrected_gamma_m = zoo::na.approx(corrected_gamma_m, rule = 2))]
    
    data[, `:=`(raw_Gamma_cm_s = zoo::na.spline(raw_Gamma_cm_s),
                corrected_Gamma_cm_s = zoo::na.spline(corrected_Gamma_cm_s))]
    
    # data[,`:=`(raw_Gamma_cm_s = raw_gamma_m * as.numeric(sample_time),
    #            corrected_Gamma_cm_s = corrected_gamma_m * as.numeric(sample_time))]
        
    data[, `:=`(raw_net_in_cm_s = raw_sy * (raw_Gamma_cm_s - raw_detrend_m_cm_s),
                corrected_net_in_cm_s = corrected_sy * (corrected_Gamma_cm_s + corrected_detrend_m_cm_s))]
    
    data[]
  }

calculate_delta_s <- 
  function(data){
    
    data[, `:=`(raw_dwt_cm_s = rolling_slope(sample_time, raw_compensated_level_cm, 3),
                corrected_dwt_cm_s = rolling_slope(sample_time, corrected_compensated_level_cm, 3))]
   
    data <- 
      smooth_data(data, "raw_dwt_cm_s", "corrected_dwt_cm_s", n = 13)
    
    data[] 
    # data[, `:=`(raw_dwt_cm_s = (shift(raw_compensated_level_cm, -1) - raw_compensated_level_cm) / 900,
    #             corrected_dwt_cm_s = (shift(corrected_compensated_level_cm, -1) - corrected_compensated_level_cm) / 900)]
  }

# The decline of water levels following rain is steeper than anticpated by the 
# basic detrending procedure. This means that for 1-2 days after rain there is
# a faster decline in water levels than would be anticipated. A simple/post-hoc
# solution may be to use a left-aligned window to detrend (which will cause 
# problems for days with preceding rainfall). A more mechanistic solution would
# be to develop a recession curve (which would likely vary with water level) and
# apply that in addition to the detrending slope. Or perhaps there is an event-
# based recession curve, but also a recession curve that can be applied to non-
# precipitation events, which is just how quickly the wetland drains, which 
# would certainly vary with water level.

# @watras-2017 was able to apply the @loheideii-2008 method because of the 
# small magnitude of seasonal water level fluctations (an the corresponding 
# small responses to preciptation)

# We have an outflow stream! Its really calculating net loss, not just ET. May
# need a PQ ratio for these wetlands to subtract out surface losses

calculate_et <- 
  function(data){
    data[, `:=`(raw_et_cm_s = raw_net_in_cm_s - raw_sy * raw_dwt_cm_s,
                corrected_et_cm_s = corrected_net_in_cm_s - corrected_sy * corrected_dwt_cm_s)]
    
    data <- 
      smooth_data(data, "raw_et_cm_s", "corrected_et_cm_s", n = 13)
    
    data[]
    
  }

adjust_water_balance <- 
  function(data){
    data[raw_et_cm_s < 0,
         `:=`(raw_et_cm_s = raw_et_cm_s - raw_et_cm_s,
              raw_net_in_cm_s = raw_net_in_cm_s - raw_et_cm_s)]
    
    data[corrected_et_cm_s < 0,
         `:=`(corrected_et_cm_s = corrected_et_cm_s - corrected_et_cm_s,
              corrected_net_in_cm_s = corrected_net_in_cm_s - corrected_et_cm_s)]
    
    data
  }

drop_trailing_data <- 
  function(data, drop.date){
    data[sample_date != as.Date(drop.date, tz = "EST5EDT")]
  }

# Helper Functions --------------------------------------------------------

three_point_slope <- 
  function(x, y){
    if(inherits(x, "POSIXt")){
      x <- as.numeric(x)
    }
    
    m1 <- 
      (shift(y, -1) - y) / (shift(x, -1) - x)
      
    m2 <- 
      (y - shift(y, 1)) / (x - shift(x, 1))
    
    (m1 + m2) / 2
  }

rolling_slope <- 
  function(x, y, n){
    if(inherits(x, "POSIXt")){
      x <- as.numeric(x)
    }
    
    if(n %% 2 == 0){
     stop("n must be odd")
    }
    
    n0 <- (n-1)/2
    
    m1 <- 
      (shift(y, -n0) - y) / (shift(x, -n0) - x)
    
    m2 <- 
      (y - shift(y, n0)) / (x - shift(x, n0))
    
    (m1 + m2) / 2
  }

loadd_raw_data <- 
  function(){
    loadd(levellogger_measurements, raw_water_data, raw_barometric_data, logger_metadata, mesonet_response,
          envir = .GlobalEnv)
  }

make_xle_files <- 
  function(levellogger_info, type = c("water", "air")){
    
    stopifnot(c("location", "serial_number") %in% names(levellogger_info))
    
    serial_numbers <- 
      levellogger_info[location %in% type, serial_number]

    path("data", serial_numbers, ext = "xle")
  }

align_loggers <- 
  function(data, x, group = "serial_number"){
    x <- 
      deparse(substitute(x))
   
    raw_mean <- 
      median(data[[x]])
    
    offsets <- 
      data[, .(offset = median(raw_mean - get(x))),
          by = group]
    
    offsets <- 
      offsets[data, on = group][["offset"]]
    
    data[[x]] + offsets
  }




merge_pressures <- 
  function(x = raw_water_data, 
           y = raw_barometric_dat){
    merge(x[, .(sample_time, 
                             water_pressure_cm,
                             water_temperature_c,
                             water_serial_number = serial_number,
                             water_pressure_error_cm = pressure_error_cm,
                             water_temperature_error_c = temperature_error_c)],
          y[, .(sample_time, 
                                  air_pressure_cm,
                                  air_temperature_c,
                                  air_serial_number = serial_number,
                                  air_pressure_error_cm = pressure_error_cm,
                                  air_temperature_error_c = temperature_error_c)],
          by = "sample_time",
          allow.cartesian = TRUE)
  }


fetch_external_baro <- 
  function(stations = c("F3311","KCMX"), 
           vars = c("air_temp", "pressure", "altimeter", "sea_level_pressure"),
           times = NULL){
    
    start_time <- 
      paste0(format(min(times), "%Y%m%d"), "0000")
    
    end_time <- 
      paste0(format(max(times), "%Y%m%d"), "2359")
    
    url <- 
      paste("https://api.synopticdata.com/v2/stations/timeseries?obtimezone=UTC",
            paste0("stid=", paste0(stations, collapse = ",")),
            paste0("vars=", paste0(vars, collapse = ",")),
            paste0("token=", Sys.getenv("SYNOPTIC_KEY")),
            paste0("start=", 202006080000),
            paste0("end=", 202006152359),
            "units=temp|C,pres|pa",
            sep = "&")
    
    response <- 
      POST(url)
    
    if(status_code(response) != 200){
      stop("Bad response from Mesonet API\n", 
           http_status(response)$message, "\n",
           url,
           call. = FALSE)
    }
    
    query_content <- 
      content(response, 
              type = "application/json")
    
    dat <- 
      lapply(seq_along(query_content[["STATION"]]), 
             function(x){
               dat <- query_content[["STATION"]][[x]]
               
               cols <- 
                 c("date_time", "air_temp_set_1", "pressure_set_1d", "altimeter_set_1", "sea_level_pressure_set_1d")
               
               # There are null values that can't just be unlisted, have to 
               # convert the NULLs to NAs first
               obs <- 
                 as.data.table(lapply(dat$OBSERVATIONS[cols], 
                                      function(z){Reduce(function(x, y){c(x, ifelse(is.null(y), NA, y))}, z)}))
               
               info <- 
                 data.table(station = dat$STID,
                            elevation_dem_ft = as.numeric(dat$ELEV_DEM),
                            elevation_ft = as.numeric(dat$ELEVATION))
               
               cbind(info, obs)})
    
    dat <- 
      rbindlist(dat)
    
    dat <- 
      dat[, .(station,
              sample_time = as.POSIXct(date_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
              elevation_cm = elevation_ft * 30.48037,
              elevation_dem_cm = elevation_dem_ft * 30.48037,
              ex_air_temperature_c = air_temp_set_1,
              ex_abs_pressure_cm = (pressure_set_1d / 98.0665) - 0.00121 * ((elevation_dem_ft * 30.48037) - 22500),
              ex_altimeter_cm = (altimeter_set_1 / 98.0665) - 0.00121 * 22500,
              ex_sea_level_pressure_cm = (sea_level_pressure_set_1d / 98.0665),
              ex_air_pressure_cm = (sea_level_pressure_set_1d / 98.0665) - 0.00121 * 22500)]
    
    attributes(dat$sample_time)$tzone <- "EST5EDT"
    
    dat[!is.na(ex_air_temperature_c) & !is.na(ex_altimeter_cm)]
  }

fit_correction_mod <- 
  function(data, 
           x = NULL,
           y = NULL,
           by = NULL,
           intercept = TRUE){
    
    if(intercept){
      form <- 
        formula(paste(y, "~", paste(x, collapse = "+")))
    } else {
      form <- 
        formula(paste(y, "~ 0 + ", paste(x, collapse = "+")))
    }
    
    
    mods <- 
      data[, .(mod = list(lm(form, data = .SD))),
           by = by]
    
    mods[, intercept := vapply(mod, get_intercept, numeric(1))]
    
    slopes <- 
      rbindlist(lapply(mods$mod, 
                       function(m) {
                         data.table(slope = list(matrix(get_slope(m), 
                                                         nrow = 1, 
                                                         dimnames = list(1, x))))
                         }))
    
    mods <- 
      cbind(mods, slopes)
    
    mods[, r2 := vapply(mod, get_r2, numeric(1))]
    
    p_values <- 
      rbindlist(lapply(mods$mod,
                       function(m) {
                         data.table(slope_p = list(matrix(get_slope_p(m),
                                                          nrow = 1,
                                                          dimnames = list(1, x))))
                       }))
    
    mods <- 
      cbind(mods, p_values)
    
    if(!intercept){
      mods[, intercept := 0]
    }
    
    if(length(mods$slope[[1]]) == 1){
      mods[, x_intercept := -intercept / slope][]
      if(!intercept){
        mods[, x_intercept := 0]
      }
    }
    
    mods
    
  }

fit_correction_rq <- 
  function(data, 
           x = NULL,
           y = NULL,
           by = NULL,
           intercept = TRUE,
           tau = NULL){
    
    if(intercept){
      mods <- 
        data[, .(mod = list(rq(get(y) ~ get(x), data = .SD, tau = tau))),
             by = by]
    } else {
      mods <- 
        data[, .(mod = list(rq(get(y) ~ 0 + get(x), data = .SD, tau = tau))),
             by = by]
    }
    
    mods[, `:=`(slope = vapply(mod, get_slope, numeric(1)),
                intercept = vapply(mod, get_intercept, numeric(1)),
                r2 = vapply(mod, get_r2, numeric(1)), 
                p_value = vapply(mod, get_slope_p, numeric(1)))]
    
    mods[, x_intercept := -intercept / slope][]
    
    if(!intercept){
      mods[, `:=`(intercept = 0,
                  x_intercept = 0)]
    }
    
    mods
    
  }

fit_correction_mmod <- 
  function(data, 
           x = NULL,
           y = NULL,
           group = NULL,
           by = NULL,
           ...){
    
    form <- 
      formula(paste0(y, " ~ ", x, " + (0 + ", x, " | ", group, ")"))
    
    mods <- 
      data[, .(mod = list(lmer(form, data = .SD, ...))),
           by = by]
    
    mods[, `:=`(slope = vapply(mod, get_slope, numeric(1), x = x),
                intercept = vapply(mod, get_intercept, numeric(1)),
                r2 = vapply(mod, get_r2, numeric(1)), 
                p_value = vapply(mod, get_slope_p, numeric(1)))]
    
    mods[, x_intercept := -intercept / slope][]
  }

get_slope <- 
  function(model){
    switch(class(model),
           lmerMod = lme4::fixef(model)[x],
           coef(model)[-c(1)])
  }

get_intercept <- 
  function(model){
    switch(class(model),
           lmerMod = lme4::fixef(model)["(Intercept)"],
           coef(model)["(Intercept)"])
  }

get_slope_p <- 
  function(model, x = "get(x)"){
    if(class(model) %in% c("lmerMod", "rq")){
      message("Objects of class lmerMod do not have p-values from tidy()")
      return(NA_real_)
    }
    
    res <- broom::tidy(model)
    
    p_values <- res[["p.value"]][which(res$term != "(Intercept)")]
    
    names(p_values) <- 
      paste0("p_", res[["term"]][which(res$term != "(Intercept)")])
    
    p_values
    }

get_r2 <- 
  function(model){
    if(class(model) %in% c("lmerMod", "rq")){
      message("Objects of class lmerMod do not have r-squared values from tidy()")
      return(NA_real_)
    }
    
    broom::glance(model)[["adj.r.squared"]]
  }

correction_residuals <- 
  function(mods, 
           original.data,
           x = NULL,
           y = NULL,
           by = NULL){
    
    
    
    
    
    if(is.null(by)){
      if(nrow(mods) != 1){
        stop("`by` must be supplied if more than one model is present")
      }
      original.data <- 
        cbind(original.data, 
              mods[, .(intercept, slope, p_value)])
    } else {
      original.data[mods, 
                    `:=`(intercept = intercept,
                         slope = slope,
                         p_value = p_value), 
                    on = strsplit(by, ",")[[1]]]
    }
    
    original.data[p_value > 0.05,
                  `:=`(intercept = intercept,
                       slope = slope)]
    
    original.data[, 
                  residuals := get(y) - (intercept + slope * get(x))]
    
    original.data[, std_residuals := as.numeric(scale(residuals)),
                  by = by]
    
    original.data[, `:=`(abs_residuals = abs(residuals),
                          abs_std_residuals = abs(std_residuals))][]
    
    original.data[, `:=`(intercept = NULL, 
                         slope = NULL,
                         p_value = NULL)]
  }

apply_correction_models <- 
  function(mods, 
           original.data,
           x = NULL,
           y = NULL,
           by = "serial_number"){
    
    mods[p_value > 0.05,
         slope := 0]
    
    original.data[mods, 
                  residuals := get(y) - (intercept + slope * get(x)), 
                  on = strsplit(by, ",")[[1]]]
    
    original.data[, std_residuals := as.numeric(scale(residuals)),
                  by = by]
    
    original.data[, abs_residuals := abs(std_residuals)][]
  }

set_value <- 
  function(data, i = NULL, j, value){
    val <- with(data, value)
    set(data, i, j, val)
  }

error_range <- 
  function(error.x, error.y){
    
    stopifnot(length(error.x) == length(error.y))
    
    df <- 
      data.frame(x = 0,
                 error.x = error.x,
                 y = 0, 
                 error.y = error.y)
    
    error.x <- 
      set_errors(df$x, 
                 df$error.x)
    
    error.y <- 
      set_errors(df$y, 
                 df$error.y)
    
    2*errors(error.x + error.y)
  }

combine_errors <- 
  function(error.x, error.y){
    
    stopifnot(length(error.x) == length(error.y))
    
    df <- 
      data.frame(x = 0,
                 error.x = error.x,
                 y = 0, 
                 error.y = error.y)
    
    error.x <- 
      set_errors(df$x, 
                 df$error.x)
    
    error.y <- 
      set_errors(df$y, 
                 df$error.y)
    
    errors(error.x + error.y)
  }

model_changepoints <- 
  function(data, ...){
    set.seed(1234)
 
    mcp(list(abs_residuals ~ 1 + temperature_difference_c,
             ~ 0 + temperature_difference_c),
        data = data,
        ...)
    
  }

extract_changepoint <- 
  function(model){
    res <- mcp::fixef(model)
    res <- subset(res, name == "cp_1", drop = TRUE)
    res[["mean"]]
  }

dy_graph <- 
  function(data, ..., grouping = NULL){
    
    if(!is.null(grouping)){
      form <- 
        as.formula(paste("sample_time ~", grouping))
      
      data_list <- 
        lapply(list(...), 
               function(x){
                 dcast(data, 
                       formula = form, 
                       value.var = x)})
      
      data_list <- 
        mapply(function(x, y){setnames(x, 
                                       names(x)[-1], 
                                       paste(names(x)[-1], y, sep = "_"))}, 
               data_list, 
               list(...),
               SIMPLIFY = FALSE)
      
      data <- 
        Reduce(merge, data_list)
      
      dy <- 
        dygraph(data, 
                names(data)[-c(1)],
                group = "all-plots")
      
    } else {
      dy <- 
        dygraph(data, ..., group = "all-plots")
    }
    dy %>% 
      dyOptions(useDataTimezone = TRUE,
                connectSeparatedPoints = TRUE) %>% 
      dyLegend(labelsSeparateLines = TRUE) %>% 
      dyRangeSelector()
  }

water_density <- 
  function(t, threshold = TRUE){
    # From [@tanaka-2001]
    
    if(threshold){
      t <- ifelse(t < -10, -10, t)
    }
    
    
    a1 <- -3.983035
    a2 <- 301.797
    a3 <- 522528.9
    a4 <- 69.34881
    a5 <- 999.97495
    
    a5 * (1 - ((t + a1)^2*(t+a2)) / (a3*(t+a4)))
  }

# adjust_density <- 
#   function(data){
#     data[,`:=`(dens_water_pressure_cm = water_pressure_cm * 1000 / water_density(water_temperature_c),
#                dens_air_pressure_cm = air_pressure_cm * 1000 / water_density(air_temperature_c))][]
#   }

adjust_density <- 
  function(pressure.cm, temperature.c){
    pressure.cm * 1000 / water_density(temperature.c)
  }


set_experiments <- 
  function(data, design.file){
    # The air pressure loggers were found wet at the end of the var-sim
    # experiment. After looking at the data extensively I determined that it 
    # must have fallen into the water at approximately 6:00 on 6/10/2020. The 
    # clearest way to see this was too look at a plot of water level error by 
    # water_temperature for just the var-sim experiment. Look at the plot in
    # output/figures/problem-baro-plot.png made using corrected air pressures.
    # with uncorrected water pressures. I compared that figure to one made with
    # the Mesonet altimeter data and found that the artifacts were removed
    
    design <- 
      fread(design.file)
    
    design <- 
      design[, data.table(sample_time = seq(as.POSIXct(start_edt, tz = "EST5EDT"), 
                                            as.POSIXct(end_edt, tz = "EST5EDT"), 
                                            by = 60)), 
             by = .(experiment)]
    
    dat <- 
      merge(data, design, by = "sample_time")
    
    dat[, dataset := fifelse(experiment == "test-dat", "training", "testing")][]
  }

slp <- 
  function(p_cm, temp_c, elevation = 225){
    # p_hPa <-
    #   0.01 * p_cm * 9.80665 * 1000
    
    p_hPa <-
      p_cm * 98.0665
    
    slp_hPa <- 
      p_hPa * (1 - (0.0065 * elevation) / (temp_c + 273.15 + 0.0065 * elevation))^-5.257
    
    # 100 * slp_hPa / (9.80665 * 1000)
    slp_hPa / 98.0665
  }

# Calculate alignments between experiments
calculate_alignments <- 
  function(data, full.data, reference.experiment, x, y = "raw_logger_error_cm", by = "serial_number"){
    reference_dat <- 
      copy(full.data[experiment == reference.experiment])
    
    reference_dat[["x2"]] <- 
      reference_dat[[x]]
    
    temp_ranges <- 
      data[,.(min_x = 0.95 * min(.SD[[x]]), 
              max_x = 1.05 * max(.SD[[x]])),
           by = by]
    
    merged_data <- 
      reference_dat[temp_ranges, 
                 on = c(by, "x2>=min_x", "x2<=max_x")]
    
    merged_data[,.(alignment_offset_cm = mean(.SD[[y]])), by = by]
  }

find_optimal_lambda <- 
  function(data, by = NULL, n = 100){
    set.seed(1234)
    
    dat <- 
      split(data, by = by)
    
    names(dat) <- 
      sub("\\.", "_", names(dat))
    
    sampled <- 
      lapply(dat,
             generate_lambda_min,
             n = n)
    
    optimal <- 
      vapply(sampled, mean, numeric(1))
    
    list(sampled = sampled,
         optimal = optimal)
  }


generate_lambda_min <- 
  function(x, n){
    Y <- 
      x$raw_error_cm
    
    X <- 
      as.matrix(x[, .(air_temperature_c, water_temperature_c, delta_at_01c_min)])
    
    replicate(n,
              cv.glmnet(y = Y, 
                        x = X, 
                        alpha = 1, 
                        lambda = c(10^seq(3, -6, by = -.1), 0))$lambda.min)
  }

# rmse <- 
#   function(y_hat, y){
#     sqrt(mean((y_hat - y)^2))
#   }

adjusted_r2 <- 
  function(y_hat, y, nobs, k){
    if(isTRUE(all.equal(as.numeric(y_hat), as.numeric(y), check.atributes = FALSE))){
      return(0)
    }
    
    1 - (1-cor(y_hat, y)^2) * (nobs - 1) / (nobs - k - 1)
  }



# experiments <- 
#   function(){
#     c("var-sim", "var-dis", "stat-sim", "stat-dis", "test-dat")
#   }

EXPERIMENTS <- 
  c("var-sim", "var-dis", "stat-sim", "stat-dis", "test-dat")

fit_matrix_files <-
  function(x = NULL){
    if(is.null(x)){
      x <- experiments()
    }
    
    files <- paste0("output/tabular/bootstrap_prediction_matrices_", x, ".rds")
    setNames(files, x)
  }

