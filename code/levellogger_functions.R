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
                         
                         # fwrite(x = data.table(model_id, rep = i, samples, inst_err, at_err, wt_err, d_at_err),
                         #        file = "output/tabular/bootstrap_samples_and_errors.csv",
                         #        append = TRUE)
                          
                         write_fst(x = data.table(model_id, rep = i, samples, inst_err, at_err, wt_err, d_at_err),
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
           keyby = .(water_sn, baro_sn, experiment, rep)]
  
  x_mat <-data[,  .(S = list(matrix(sample_time,
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
  
  for(Ex in EXPERIMENTS){
    
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
              c(list(water_sn = water_sn, baro_sn = baro_sn, experiment = experiment),
                lapply(.SD, function(x){do.call(rbind, x)})),
              .SDcols = c("S", "E_hat")]
        
        setnames(post_pred,
                 names(post_pred),
                 sub("^.*\\.", "", names(post_pred)))
        
        post_pred[, sample_time := as.POSIXct(sample_time, tz = "EST5EDT",
                                              origin = "1970-01-01")]
        
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

# Helper Functions --------------------------------------------------------


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
    
    res <- tidy(model)
    
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
    
    glance(model)[["adj.r.squared"]]
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
    res <- fixef(model)
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

