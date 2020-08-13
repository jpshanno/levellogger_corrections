source("code/levellogger_packages.R"); source("code/levellogger_functions.R")

loadd(combined_data)

combined_data <- 
  Reduce(rbind, combined_data)

pb <-
  txtProgressBar(min = 0, max = 5)

c <- 0

predictions <-
  combined_data[, .(water_sn,
                    baro_sn,
                    experiment,
                    sample_time,
                    predicted_error_cm = NA_real_,
                    lower_bound_cm = NA_real_,
                    upper_bound_cm = NA_real_)]

setkey(predictions, "water_sn", "baro_sn", "experiment", "sample_time")

for(Ex in unique(predictions$experiment)){
  test_pred <-
    readRDS(paste0("output/tabular/bootstrap_prediction_matrices_", Ex, ".rds"))
  
  setkey(test_pred, "water_sn", "baro_sn", "experiment")
  
  for(i in unique(test_pred$water_sn)){
    for(j in unique(test_pred$baro_sn)){
      post_pred <-
        test_pred[CJ(i,j),
                  c(list(water_sn = water_sn, baro_sn = baro_sn, experiment = experiment),
                    lapply(.SD, function(x){do.call(rbind, x)})),
                  .SDcols = c("S", "E_hat")]
      
      setnames(post_pred,
               names(post_pred),
               sub("^.*\\.", "", names(post_pred)))
      
      post_pred[, sample_time := as.POSIXct(sample_time, tz = "EST5EDT",
                                            origin = "1970-01-01")]
      
      write_fst(post_pred,
                path = paste0("output/tabular/bootstrap_predictions_uncertainty_test.df/",
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
      
      rm(post_pred, ints)
      
      c <- c + 1/36
      
      setTxtProgressBar(pb, c)
    }
  }
}

setkey(combined_data, "water_sn", "baro_sn", "experiment", "sample_time")
setkey(predictions, "water_sn", "baro_sn", "experiment", "sample_time")

predictions[combined_data,
            `:=`(experiment = i.experiment,
                 water_depth_cm = i.water_depth_cm,
                 raw_water_level_cm = i.raw_water_level_cm,
                 raw_error_cm = i.raw_error_cm,
                 instrument_error_cm = i.instrument_error_cm)]

fwrite(predictions,
       "output/tabular/bootstrap_predictions_uncertainty_test.csv")
