library(data.table)
library(glmnet)
library(broom)
library(errors)
library(parallel)
options(mc.cores = 4)
set.seed(1234)

drake::loadd(compensated_data)

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

training_data <- 
  compensated_data[experiment == "test-dat"]

training_data[, delta_at_error_c_min := combine_errors(100 * air_temp_error_c, 100 * air_temp_error_c)]
training_data[, delta_at_c_min := 100 * delta_at_c_min]

dat <-
  split(training_data,
        by = c("water_sn", "baro_sn"))

lambdas <- 
  {set.seed(1234)
    lapply(split(training_data, by = "baro_sn"), 
           function(x){vapply(1:100, 
                              function(y){cv.glmnet(y = x$raw_error_cm, 
                                                    x = as.matrix(x[, .(air_temperature_c, water_temperature_c, delta_at_c_min)]), 
                                                    alpha = 1, 
                                                    lambda = c(10^seq(3, -6, by = -.1), 0))$lambda.min}, 
                              FUN.VALUE = numeric(1))})}

opt_lambdas <- 
  vapply(lambdas, median, numeric(1), USE.NAMES = TRUE)

# bigmods <- 
#   cv.glmnet(y = training_data$raw_error_cm, 
#             x = as.matrix(training_data[, .(air_temperature_c, water_temperature_c, delta_at_c_min)]), 
#             alpha = 1, 
#             lambda = c(10^seq(3, -6, by = -.1), 0))
# 
# opt_lambda <- 
#   bigmod$lambda.min

mods <- 
  mclapply(dat, 
           function(x){
             
             i <- 1
             
             models <- 
               vector("list", 1000)
             
             wsn <- 
               unique(x$water_sn)
             
             bsn <- 
               unique(x$baro_sn)
             
             model_id <- 
               paste(wsn, bsn, sep = "_")
             
             opt_lambda <- 
               opt_lambdas[["bsn"]]
             
             nobs <- 
               nrow(x)
             
             while(i < length(models) + 1){
               
               samples <- 
                 sample.int(n = nobs, size = nobs, replace = TRUE)
               
               inst_err <- 
                 rnorm(nobs, 0, x$instrument_error_cm)
               
               at_err <- 
                 rnorm(nobs, 0, x$air_error_c)
               
               wt_err <- 
                 rnorm(nobs, 0, x$water_error_c)
               
               d_at_err <- 
                 rnorm(nobs, 0, x$delta_at_error_c_min)
               
               Y<- 
                 x$raw_error_cm[samples] + inst_err
               
               X <- 
                 matrix(c(air_temperature_c = x$air_temperature_c[samples] + at_err,
                          water_temperature_c = x$water_temperature_c[samples] + wt_err,
                          delta_at_c_min = x$delta_at_c_min[samples] + d_at_err),
                        nrow = nobs,
                        dimnames = list(NULL, 
                                        c("air_temperature_c",
                                          "water_temperature_c",
                                          "delta_at_c_min")))
               
               mod <- 
                 glmnet(x = X,
                        y = Y,
                        alpha = 1,
                        lambda = opt_lambda)
               
               adj_r2 <- 
                 1 - (1-cor(predict(mod, X, s = opt_lambda)[,1], Y)) * (nobs - 1) / (nobs - 3 - 1)
                 
               
               models[[i]] <- 
                 list(mod_summary = data.table(model_id, 
                                               rep = i, 
                                               tidy(mod)),
                      mod_fit = data.table(model_id, 
                                           rep = i, 
                                           sigma = sigma(mod), 
                                           adj_r2,
                                           nobs,
                                           deviance = deviance(mod)))
               
               i <- i + 1
             }
             models
           }
           )

RPushbullet::pbPost("note", title = "Bootstrapping Complete")

model_summaries <- 
  rbindlist(unlist(lapply(mods, function(x){lapply(x, `[[`, 1)}), recursive = FALSE))

model_summaries[, c("water_sn", "baro_sn") := tstrsplit(model_id, split = "_")]

# ggplot(data = model_summaries,
#        aes(x = estimate,
#            color = baro_sn,
#            group = model_id)) +
#   geom_density(aes(y = ..scaled..)) +
#   facet_wrap(~term, scales = "free")

# ggplot(data = model_summaries[term == "air_temperature_c"],
#        aes(x = estimate,
#            group = model_id)) +
#   geom_density(aes(y = ..scaled..)) +
#   facet_grid(water_sn ~ baro_sn)

fwrite(model_summaries,
       file = "output/tabular/bootstrap_model_summaries.csv")

RPushbullet::pbPost("note", title = "Model summaries saved")

model_fits <- 
  rbindlist(unlist(lapply(mods, function(x){lapply(x, `[[`, 2)}), recursive = FALSE))

model_fits[, c("water_sn", "baro_sn") := tstrsplit(model_id, split = "_")]

# ggplot(data = model_fits,
#        aes(x = sigma)) +
#   geom_density(aes(y = ..scaled..)) +
#   facet_grid(water_sn ~ baro_sn)

# drake::loadd(logger_metadata)
# logger_metadata[serial_number %in% c("2025928", "2013939")]

fwrite(model_fits,
       file = "output/tabular/bootstrap_model_fits.csv")

RPushbullet::pbPost("note", title = "Model fit info saved")

predictors <- 
  dcast(rbind(model_fits[, .(model_id, water_sn, baro_sn, rep, term = "sigma", estimate = sigma)], 
              model_summaries[, .(model_id, water_sn, baro_sn, rep, term, estimate)], 
              fill = TRUE),
        model_id + water_sn + baro_sn + rep ~ term,
        value.var = "estimate")

predictors[unique(model_fits[, .(model_id, nobs)]),
           nobs := nobs,
           on = "model_id"]

predictors[, `:=`(water_sn = as.character(water_sn),
                  baro_sn = as.character(baro_sn))]

setnames(predictors, 
         c("(Intercept)", "rep"),
         c("intercept", "replication"))

pred_mat <- 
  predictors[, .(intercept = list(matrix(rep(intercept, nobs),
                                         dimnames = list(NULL, "intercept"),
                                         ncol = 1)),
                 B = list(matrix(c(air_temperature_c,
                                   water_temperature_c,
                                   delta_at_c_min),
                                 dimnames = list(c("b_at",
                                                   "b_wt",
                                                   "b_dat"),
                                                 NULL),
                                 nrow = 3)),
                 sigma = list(matrix(rep(sigma, nobs),
                                     dimnames = list(NULL, "sigma"),
                                     ncol = 1))),
             by = .(water_sn, baro_sn, replication)]

mod_mat <- 
  training_data[, .(S = list(matrix(sample_time, 
                                    dimnames = list(NULL, "sample_time"),
                                    ncol = 1)),
                    L = list(matrix(raw_water_level_cm, 
                                    dimnames = list(NULL, "raw_water_level_cm"),
                                    ncol = 1)),
                    D = list(matrix(water_depth_cm, 
                                    dimnames = list(NULL, "water_depth_cm"),
                                    ncol = 1)),
                    E = list(matrix(raw_error_cm, 
                                    dimnames = list(NULL, "raw_error_cm"),
                                    ncol = 1)), 
                    X = list(matrix(c(air_temperature_c, 
                                      water_temperature_c, 
                                      delta_at_c_min), 
                                    dimnames = list(NULL, c("air_temperature_c", 
                                                            "water_temperature_c", 
                                                            "delta_at_c_min")),
                                    ncol = 3)),
                    nobs = .N),
                by = .(water_sn, baro_sn)]

dat <- 
  mod_mat[pred_mat, on = c("water_sn", "baro_sn")]

# This will provide prediction intervals for E_hat, removing sigma will provide
# confidence intervals for E_hat (when quantiles are taken)
# https://stackoverflow.com/questions/10584009/confidence-intervals-on-predictions-for-a-bayesian-linear-regression-model/10590067#10590067
dat[, E_hat := Map(function(x, i, b, s){y_hat <- i + x %*% b; matrix(rnorm(nrow(i), y_hat, s), dimnames = list(NULL, "predicted_error_cm"))}, X, intercept, B, sigma)]
dat[, L_rect := Map(`-`, L, E_hat)]
dat[, R := Map(`-`, E, E_hat)]
dat[, water_sn := Map(function(x, y){matrix(x,
                                            dimnames = list(NULL, "water_sn"),
                                            ncol = 1,
                                            nrow = y)},
                      water_sn,
                      nobs)]
dat[, baro_sn := Map(function(x, y){matrix(x,
                                            dimnames = list(NULL, "baro_sn"),
                                            ncol = 1,
                                           nrow = y)},
                      baro_sn,
                      nobs)]
# dat[, replication := Map(function(x, y){matrix(rep(x, y),
#                                            dimnames = list(NULL, "replication"),
#                                            ncol = 1)},
#                          replication,
#                          nobs)]

# Tried using Reduce, do.call, and rbindlist(as.data.table). rbind() showed the
# smallest increase as I increased the number of rows tested
# microbenchmark::microbenchmark(reduce = dat[1:10, lapply(.SD, Reduce, f = rbind),
#                                             .SDcols = c("water_sn", "baro_sn", "replication", "S", "L", "D", 
#                                                         "E_hat", "L_rect", "R", "sigma")],
#                                rbind = dat[1:10, lapply(.SD, function(x){do.call(rbind, x)}),
#                                            .SDcols = c("water_sn", "baro_sn", "replication", "S", "L", "D", 
#                                                        "E_hat", "L_rect", "R", "sigma")],
#                                DT = dat[1:10, lapply(.SD, function(x){rbindlist(lapply(x, as.data.table))}),
#                                         .SDcols = c("water_sn", "baro_sn", "replication", "S", "L", "D", 
#                                                     "E_hat", "L_rect", "R", "sigma")])


out <- 
  dat[, lapply(.SD, function(x){do.call(rbind, x)}),
      .SDcols = c("water_sn", "baro_sn", "S", "E_hat")]

setnames(out,
         names(out),
         sub("^.*\\.", "", names(out)))

out[, sample_time := as.POSIXct(sample_time, tz = "EST5EDT", origin = "1970-01-01")]

ints <- 
  out[, .(pred_lwr = quantile(predicted_error_cm, probs = 0.025),
          fit = quantile(predicted_error_cm, probs = 0.5), 
          pred_upr = quantile(predicted_error_cm, probs = 0.975)),
      by = .(water_sn, baro_sn, sample_time)]

ints[training_data,
     `:=`(experiment = experiment,
          water_depth_cm = water_depth_cm,
          raw_water_level_cm = raw_water_level_cm,
          raw_error_cm = raw_error_cm,
          instrument_error_cm = instrument_error_cm),
     on = c("water_sn", "baro_sn", "sample_time")]

ints[, rect_water_level_cm := raw_water_level_cm - fit]
ints[, `:=`(raw_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
            rect_water_level_cm = rect_water_level_cm - mean(rect_water_level_cm - water_depth_cm)),
     by = .(water_sn, baro_sn)]

library(ggplot2)
ggplot(data = ints[water_sn == "1033239"],
       aes(x = sample_time,
           y = rect_water_level_cm)) +
  geom_ribbon(aes(ymin = water_depth_cm - (fit - pred_lwr),
                  ymax = water_depth_cm + (pred_upr - fit)),
              fill = "gray20") + 
  geom_ribbon(aes(ymin = water_depth_cm - (qnorm(0.975) * instrument_error_cm),
                  ymax = water_depth_cm + qnorm(0.975) * instrument_error_cm),
              fill = "gray80") + 
  geom_line() +
  # geom_line(aes(y = raw_water_level_cm),
  #           linetype = "dotted") +
  geom_line(aes(y = water_depth_cm)) +
  facet_wrap(~baro_sn, scales = "free")

mod <- 
  lm(raw_error_cm ~ air_temperature_c + water_temperature_c + delta_at_c_min, 
     data = split(training_data, by = c("water_sn", "baro_sn"))[[1]])

ints[, .(pred_range = median(pred_upr - pred_lwr), 
         instrument_error_range = median(2 * qnorm(0.975) * instrument_error_cm),
         instrument_and_sigma_range = median(2 * qnorm(0.975) * sqrt(instrument_error_cm^2 + sigma(mod)^2))), 
     by = .(water_sn, baro_sn)]

ints[, .(pred_range = median(pred_upr - pred_lwr), 
         instrument_error_range = median(2 * sqrt((qnorm(0.975) * instrument_error_cm)^2 + sigma(mod)^2))), 
     by = .(water_sn, baro_sn)]


# Testing Data ------------------------------------------------------------

# library(data.table)
library(fst)
# library(disk.frame)
setup_disk.frame(4)
options(future.globals.maxSize = Inf)
setDTthreads(2)
threads_fst(2)

model_summaries <- 
  fread("output/tabular/bootstrap_model_summaries.csv")

model_fits <- 
  fread("output/tabular/bootstrap_model_fits.csv")

drake::loadd(compensated_data)

testing_data <- 
  compensated_data[experiment != "test-dat"]

testing_mat <- 
  testing_data[, .(S = list(matrix(sample_time, 
                                   dimnames = list(NULL, "sample_time"),
                                   ncol = 1)), 
                   X = list(matrix(c(air_temperature_c, 
                                     water_temperature_c, 
                                     delta_at_c_min), 
                                   dimnames = list(NULL, c("air_temperature_c", 
                                                           "water_temperature_c", 
                                                           "delta_at_c_min")),
                                   ncol = 3)),
                   nobs = .N),
               by = .(water_sn, baro_sn)]

# model_summaries[, estimate := fifelse(p.value > 0.05, 0, estimate)]

predictors <- 
  dcast(rbind(model_fits[, .(model_id, water_sn, baro_sn, rep, term = "sigma", estimate = sigma)], 
              model_summaries[, .(model_id, water_sn, baro_sn, rep, term, estimate)], 
              fill = TRUE),
        model_id + water_sn + baro_sn + rep ~ term,
        value.var = "estimate")

# predictors <- 
#   predictors[, lapply(.SD, rep, 4)]
# 
# predictors[, experiment := rep(unique(testing_data$experiment), each = nrow(predictors)/4)]

predictors[, `:=`(water_sn = as.character(water_sn),
                  baro_sn = as.character(baro_sn))]

predictors[testing_data[, .(nobs = .N), by = .(water_sn, baro_sn)],
           nobs := nobs,
           on = c("water_sn", "baro_sn")]

setnames(predictors, 
         c("(Intercept)", "rep"),
         c("intercept", "replication"))

pred_mat <- 
  predictors[, .(intercept = list(matrix(rep(intercept, nobs),
                                         dimnames = list(NULL, "intercept"),
                                         ncol = 1)),
                 B = list(matrix(c(air_temperature_c,
                                   water_temperature_c,
                                   delta_at_c_min),
                                 dimnames = list(c("b_at",
                                                   "b_wt",
                                                   "b_dat"),
                                                 NULL),
                                 nrow = 3)),
                 sigma = list(matrix(rep(sigma, nobs),
                                     dimnames = list(NULL, "sigma"),
                                     ncol = 1))),
             by = .(water_sn, baro_sn, replication)]

test_pred <- 
  testing_mat[pred_mat, on = c("water_sn", "baro_sn")]

# This will provide prediction intervals for E_hat, removing sigma will provide
# confidence intervals for E_hat (when quantiles are taken)
# https://stackoverflow.com/questions/10584009/confidence-intervals-on-predictions-for-a-bayesian-linear-regression-model/10590067#10590067
test_pred[, E_hat := Map(function(x, i, b, s){
                        y_hat <- i + x %*% b
                        matrix(rnorm(nrow(i), y_hat, s), 
                               dimnames = list(NULL, "predicted_error_cm"))
                        }, X, intercept, B, sigma)]

setkey(test_pred, "water_sn", "baro_sn")

predictions <-
  testing_data[, .(water_sn,
                   baro_sn,
                   sample_time,
                   predicted_error_cm = NA_real_,
                   lower_bound_cm = NA_real_,
                   upper_bound_cm = NA_real_)]

setkey(predictions, "water_sn", "baro_sn", "sample_time")

pb <- 
  txtProgressBar(min = 0, max = 36)

c <- 0

for(i in unique(test_pred$water_sn)){
  for(j in unique(test_pred$baro_sn)){
    post_pred <- 
      test_pred[CJ(i,j), 
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

# posterior_pred <-
#   disk.frame("output/tabular/bootstrap_predictions_test_data.df")
# 
# predictions <-
#   posterior_pred[, .(predicted_error_cm = quantile(predicted_error_cm, probs = 0.5),
#                      lower_bound_cm = quantile(predicted_error_cm, probs = 0.025),
#                      upper_bound_cm = quantile(predicted_error_cm, probs = 0.975)),
#                  by = .(water_sn, baro_sn, sample_time),
#                  keep = c("water_sn", "baro_sn", "sample_time", "predicted_error_cm")]

setkey(testing_data, "water_sn", "baro_sn", "sample_time")
setkey(predictions, "water_sn", "baro_sn", "sample_time")

predictions[testing_data,
            `:=`(experiment = i.experiment,
                 water_depth_cm = i.water_depth_cm,
                 raw_water_level_cm = i.raw_water_level_cm,
                 raw_error_cm = i.raw_error_cm,
                 instrument_error_cm = i.instrument_error_cm)]

# smoothScatter((predictions$predicted_error_cm - predictions$raw_error_cm) ~ predictions$raw_error_cm)
# smoothScatter(predictions$predicted_error_cm ~ predictions$raw_error_cm); abline(0, 1)
# cor(predictions$predicted_error_cm, predictions$raw_error_cm)
# fixest::feols(predicted_error_cm ~ raw_error_cm, data = predictions)

predictions[, rect_water_level_cm := raw_water_level_cm - predicted_error_cm]
predictions[, `:=`(raw_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
                   rect_water_level_cm = rect_water_level_cm - mean(rect_water_level_cm - water_depth_cm)),
            by = .(water_sn, baro_sn, experiment)]

mad_values <- 
  melt(predictions[, lapply(.SD, function(x) mad(x - water_depth_cm)), 
            by = .(water_sn, baro_sn, experiment),
            .SDcols = c("raw_water_level_cm", "rect_water_level_cm")], 
       id.vars = c("water_sn", "baro_sn", "experiment"),
       measure.vars = c("raw_water_level_cm", "rect_water_level_cm"),
       variable.name = "type", 
       value.name = "mad_cm")

rmse_values <- 
  melt(predictions[, lapply(.SD, function(x) sqrt(mean((x - water_depth_cm)^2))), 
                   by = .(water_sn, baro_sn, experiment),
                   .SDcols = c("raw_water_level_cm", "rect_water_level_cm")], 
       id.vars = c("water_sn", "baro_sn", "experiment"),
       measure.vars = c("raw_water_level_cm", "rect_water_level_cm"),
       variable.name = "type", 
       value.name = "rmse_cm")

library(ggplot2)
library(ggdist)
library(emmeans)

ggplot(rmse_values,
       aes(y = experiment,
           x = rmse_cm,
           fill = type)) + 
  stat_gradientinterval() +
  facet_wrap(~baro_sn)

rmse_mod <- 
  lm(rmse_cm ~ type*experiment, data = rmse_values)

emmeans(rmse_mod, pairwise ~ type | experiment, adjust = "bonferroni")

predictions[, `:=`(raw_ooir = !between(raw_water_level_cm, 
                                     water_depth_cm - qnorm(0.975) * instrument_error_cm,
                                     water_depth_cm + qnorm(0.975) * instrument_error_cm),
                   rect_ooir = !between(rect_water_level_cm, 
                                       water_depth_cm - qnorm(0.975) * instrument_error_cm,
                                       water_depth_cm + qnorm(0.975) * instrument_error_cm),
                   raw_oopr = !between(raw_water_level_cm, 
                                      water_depth_cm - (predicted_error_cm - lower_bound_cm),
                                      water_depth_cm + (upper_bound_cm - predicted_error_cm)),
                   rect_oopr = !between(rect_water_level_cm, 
                                       water_depth_cm - (predicted_error_cm - lower_bound_cm),
                                       water_depth_cm + (upper_bound_cm - predicted_error_cm)))]

ooir_mod <- 
  lm(value ~ variable*experiment,
     data = melt(predictions, id.vars = c("water_sn", "baro_sn", "sample_time", "experiment"),
                 measure.vars = c("raw_ooir", "rect_ooir")))

emmeans(ooir_mod, pairwise ~ variable | experiment, adjust = "bonferroni")

oor_summary <- 
  predictions[, lapply(.SD, sum),
              by = .(water_sn, baro_sn, experiment),
              .SDcols = patterns("_oo")]

oor_summary[, `:=`(delta_ooir = rect_ooir - raw_ooir,
                   delta_oopr = rect_oopr - raw_oopr)]

oor_summary[, .(reduced_ooir = sum(delta_ooir < 0),
                increased_ooir = sum(delta_ooir > 0),
                reduced_oopr = sum(delta_oopr < 0),
                increased_oopr = sum(delta_oopr > 0)), 
            by = .(experiment)]

ggplot(oor_summary[oor_summary$rect_ooir != 0 & oor_summary$raw_ooir != 0],
       aes(x = delta_ooir)) +
  geom_histogram(bins = 36) +
  facet_wrap(~experiment)

hist(oor_summary$delta_ooir[oor_summary$rect_ooir != 0 & oor_summary$raw_ooir != 0], nclass = 144)

ggplot(data = predictions[experiment == "var-dis"],
       aes(x = sample_time,
           y = rect_water_level_cm)) +
  # geom_ribbon(aes(ymin = water_depth_cm - (predicted_error_cm - lower_bound_cm),
  #                 ymax = water_depth_cm + (upper_bound_cm - predicted_error_cm)),
  #             fill = "gray20") + 
  geom_ribbon(aes(ymin = water_depth_cm - (qnorm(0.975) * instrument_error_cm),
                  ymax = water_depth_cm + qnorm(0.975) * instrument_error_cm),
              fill = "gray80") + 
  geom_line(aes(color = baro_sn)) +
  # geom_line(aes(y = raw_water_level_cm),
  #           linetype = "dotted") +
  geom_line(aes(y = water_depth_cm)) +
  facet_wrap(~ water_sn, scales = "free")

mod <- 
  lm(raw_error_cm ~ air_temperature_c + water_temperature_c + delta_at_c_min, 
     data = split(training_data, by = c("water_sn", "baro_sn"))[[1]])

ints[, .(pred_range = median(pred_upr - pred_lwr), 
         instrument_error_range = median(2 * qnorm(0.975) * instrument_error_cm),
         instrument_and_sigma_range = median(2 * qnorm(0.975) * sqrt(instrument_error_cm^2 + sigma(mod)^2))), 
     by = .(water_sn, baro_sn)]

# Test if uncertainty is higher in different settings
unc_dat <- 
  melt(predictions[, .(pred_range = median(upper_bound_cm - lower_bound_cm), 
                       instrument_error_range = median(2 * qnorm(0.975) * instrument_error_cm)), 
              by = .(water_sn, baro_sn, experiment)],
       measure.vars = c("pred_range", "instrument_error_range"))

pred_error_mod <- 
  lm(value ~ experiment, data = unc_dat[variable == "pred_range"])

emmeans(pred_error_mod, pairwise ~ experiment, adjust = "bonferroni")

error_type_mod <- 
  lm(value ~ variable*experiment, data = unc_dat)

emmeans(error_type_mod, pairwise ~ variable | experiment, adjust = "bonferroni")
