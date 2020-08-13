source("code/levellogger_packages.R")
source("code/levellogger_functions.R")
loadd(bootstrap_models, combined_data)

pred_mat <-
    bootstrap_models[, .(B = list(matrix(c(intercept,
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
  #
  #
  x_mat <- Reduce(rbind, combined_data)[,  .(S = list(matrix(sample_time,
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

  pred_mat[x_mat,
           `:=`(S = i.S,
                X = i.X)]

  rm(bootstrap_models, x_mat)
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

  for(Ex in c("var-sim", "var-dis", "stat-sim", "stat-dis", "test-dat")){
    saveRDS(pred_mat[experiment == Ex],
            paste0("output/tabular/bootstrap_prediction_matrices_", Ex, ".rds"))}