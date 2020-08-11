loadd_raw_data <- 
  function(){
    loadd(levellogger_measurements, raw_water_data, raw_barometric_data, logger_metadata, mesonet_response,
          envir = .GlobalEnv)
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

rmse <- 
  function(y_hat, y){
    sqrt(mean((y_hat - y)^2))
  }

adjusted_r2 <- 
  function(y_hat, y, nobs, k){
    1 - (1-cor(y_hat, y)) * (nobs - 1) / (nobs - k - 1)
  }
