source("code/levellogger_packages.R")
source("code/levellogger_functions.R")

loadd(bootstrap_models, combined_data)

pred_mat <-
    bootstrap_models[experiment == "test-dat", 
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
                       sigma = list(matrix(rep(sigma, nobs),
                                           dimnames = list(NULL, "sigma"),
                                           ncol = 1))),
                     keyby = .(water_sn, baro_sn, rep)]


x_mat <- 
  combined_data$testing[,  .(S = list(matrix(sample_time,
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
                         keyby = .(water_sn, baro_sn)]

pred_mat[x_mat,
         `:=`(S = i.S,
              X = i.X)]

# matrix(rnorm(nrow(i), y_hat, s), dimnames = list(NULL, "predicted_error_cm"))}
pred_mat[, E_hat := Map(function(x, b, s){matrix(rnorm(n = nrow(x),
                                                       mean = x %*% b,
                                                       sd = s),
                                                 ncol = 1,
                                                 dimnames = list(NULL,
                                                                 "predicted_error_cm"))},
                        x = X,
                        b = B,
                        s = sigma)]

setkey(pred_mat, "water_sn", "baro_sn")

predictions <-
  combined_data$testing[, .(water_sn,
                            baro_sn,
                            sample_time,
                            predicted_error_cm = NA_real_,
                            lower_bound_cm = NA_real_,
                            upper_bound_cm = NA_real_)]

setkey(predictions, "water_sn", "baro_sn", "sample_time")

pb <- 
  txtProgressBar(min = 0, max = 36)

c <- 0

for(i in unique(pred_mat$water_sn)){
  for(j in unique(pred_mat$baro_sn)){
    post_pred <- 
      pred_mat[CJ(i,j), 
                c(list(water_sn = water_sn, baro_sn = baro_sn),
                  lapply(.SD, function(x){do.call(rbind, x)})),
                .SDcols = c("S", "E_hat")]
    
    setnames(post_pred,
             names(post_pred),
             sub("^.*\\.", "", names(post_pred)))
    
    post_pred[, sample_time := as.POSIXct(sample_time, tz = "EST5EDT", 
                                          origin = "1970-01-01")]
    
    write_fst(post_pred,
              path = paste0("output/tabular/bootstrap_predictions_test_data.df/",
                            i, "_", j, ".fst"))
    
    ints <-
      post_pred[, .(pred_lwr = quantile(predicted_error_cm, probs = 0.025),
                    fit = quantile(predicted_error_cm, probs = 0.5),
                    pred_upr = quantile(predicted_error_cm, probs = 0.975)),
                keyby = .(water_sn, baro_sn, sample_time)]
    
    predictions[ints,
                `:=`(predicted_error_cm = fit,
                     lower_bound_cm = pred_lwr,
                     upper_bound_cm = pred_upr)]
    
    rm(post_pred, ints)
    
    c <- c + 1
    
    setTxtProgressBar(pb, c)
  }
}


setkey(combined_data$testing, "water_sn", "baro_sn", "sample_time")

predictions[combined_data$testing,
            `:=`(experiment = i.experiment,
                 water_depth_cm = i.water_depth_cm,
                 raw_water_level_cm = i.raw_water_level_cm,
                 raw_error_cm = i.raw_error_cm,
                 instrument_error_cm = i.instrument_error_cm)]

predictions[, rect_water_level_cm := raw_water_level_cm - predicted_error_cm]
predictions[, `:=`(raw_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
                   rect_water_level_cm = rect_water_level_cm - mean(rect_water_level_cm - water_depth_cm)),
            by = .(water_sn, baro_sn, experiment)]

# Weird problem where between() (data.table or dplyr) was returning TRUE when FALSE
predictions[, `:=`(instrument_lower = water_depth_cm - qnorm(0.975) * instrument_error_cm,
                   instrument_upper = water_depth_cm + qnorm(0.975) * instrument_error_cm,
                   propagated_lower = water_depth_cm - (predicted_error_cm - lower_bound_cm),
                   propagated_upper = water_depth_cm + (upper_bound_cm - predicted_error_cm))]

fwrite(predictions, 
       "output/tabular/bootstrap_prediction_values.csv")
