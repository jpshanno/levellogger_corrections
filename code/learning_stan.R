source("code/levellogger_packages.R"); source("code/levellogger_functions.R")
theme_set(ggthemes::theme_few())
brown_green_scale <- # Inspired by NOAA maps
  c('#97601c', '#a47138', '#b08353', '#bc956e', '#c6a78a', '#d0baa6', '#d8cdc3', '#e0e0e0', '#c6d1cf', '#adc2bf', '#93b3af', '#7aa49f', '#5f9690', '#438781', '#207972')

blue_orange_scale <-
  c('#3085c7', '#4a98d4', '#6babde', '#90bee3', '#b7cfe3', '#e0e0e0', '#eac79d', '#eaae66', '#e5943a', '#dd7a15', '#d55e00')

pale_pal <- 
  c("#7DA050", # green
    "#D19648", # orange
    "#329985", # teal
    "#5982A0", # blue
    "#966283", # purple
    "#9B5249") # red
loadd(combined_data)

# melt(compensated_data,
#      id.vars = c("baro_sn", "water_sn", "experiment"),
#      measure.vars = c("air_temperature_c", "water_temperature_c")) %>% 
#   rbind(., copy(.)[, experiment:= "combined"]) %>%
#   ggplot(aes(x = experiment, y = value)) + 
#   geom_boxplot() + 
#   facet_wrap(~variable, scales = "free_y")

training_data <- 
  combined_data$training

test_data <- 
  combined_data$testing

lm_mods <- 
  fit_correction_mod(training_data,
                     y = "raw_error_cm",
                     x = c("air_temperature_c", 
                           "water_temperature_c", 
                           "delta_at_01c_min"),
                     by = c("baro_sn", "water_sn"))

lm_mods[, raw_slope := lapply(slope, t)]
lm_mods[, slope_p := lapply(slope_p, t)]

library(glmnet)

loadd(lambdas)

opt_lambdas <- 
  lambdas$optimal
  
mods <- 
  training_data[, {
    Y <- 
      as.matrix(.SD[, .(raw_error_cm)])
    
    X <- 
      as.matrix(.SD[, .(air_temperature_c, 
                        water_temperature_c,
                        delta_at_01c_min)])
    
    opt_lambda <- 
      opt_lambdas[[paste(.BY[[1]], .BY[[2]], sep = "_")]]
    
    rmod <-
      # cv.glmnet(x = X,
      #           y = Y,
      #           alpha = 1,
      #           lambdas = c(10^seq(3, -6, by = -.1), 0))
      glmnet(x = X,
             y = Y,
             alpha = 1,
             lambda = opt_lambda)

    # opt_lambda <-
      # rmod$lambda.min
      # rmod$lambda[max(which(rmod$dev.ratio < 0.95))]
    
    # rmod <- 
    #   cv.glmnet(x = X,
    #             y = Y,
    #             alpha = 0,
    #             lambda = c(10^seq(2, -3, -0.1), 0))
    # 
    # opt_lambda <- 
    #   rmod$lambda.1se
    
    # if(opt_lambda == rcv$lambda.min){
    #   opt_lambda <-
    #     NA_real_
    #   message("Minimum lambda == 1 SE lambda. Check fit")
    # }

    coefs <- 
      coef(rmod, s = opt_lambda)
    
    list(mod = list(rmod),
         lambda = opt_lambda,
         intercept = list(matrix(coefs[1])),
         slope = list(matrix(coefs[-c(1)])))
  },
  by = .(water_sn, baro_sn)]


lm_mods[, is_sig := lapply(slope_p, function(x){as.numeric(x <= 0.05)})]
lm_mods[, slope := Map(`*`, raw_slope, is_sig)]

training_data[mods, 
          `:=`(intercept = i.intercept,
               slope = i.slope),
          on = c("baro_sn", "water_sn")]

training_data[, X := list(list(c(air_temperature_c, water_temperature_c, delta_at_01c_min))),
          by = .(water_sn, baro_sn, sample_time)]

training_data[, 
          correction_cm := mapply(function(i, b, x){i + x %*% b}, i = intercept, b = slope, x = X),
          by = .(water_sn, baro_sn, sample_time)]

training_data[, 
          `:=`(centered_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
               corrected_water_level_cm = (raw_water_level_cm - correction_cm) - mean((raw_water_level_cm - correction_cm) - water_depth_cm)),
          by = .(water_sn, baro_sn, experiment)]

training_data[, std_residuals := scale(correction_cm - raw_error_cm)]

# Normality of Residuals
ggplot(data = training_data, 
       aes(x = std_residuals)) + 
  geom_density() + 
  geom_density(data = NULL, 
               aes(x = rnorm(nrow(training_data))), 
               color = "red") + 
  facet_wrap(~water_sn + baro_sn, 
             scales = "free")

ggplot(data = training_data, 
       aes(x = water_temperature_c,
           y = std_residuals)) + 
  geom_point() + 
  facet_wrap(~water_sn + baro_sn, 
             scales = "free")

ggplot(data = training_data, 
       aes(x = air_temperature_c,
           y = std_residuals)) + 
  geom_point() + 
  facet_wrap(~water_sn + baro_sn, 
             scales = "free")

ggplot(data = training_data, 
       aes(x = delta_at_01c_min,
           y = std_residuals)) + 
  geom_point() + 
  facet_wrap(~water_sn + baro_sn, 
             scales = "free")

ggplot(data = training_data, 
       aes(x = delta_wt_c_min,
           y = std_residuals)) + 
  geom_point() + 
  facet_wrap(~water_sn + baro_sn, 
             scales = "free")

# Running plot(mods$mod[[1]]) shows a pattern in residuals ~ fitted, which looks
# like an artificat of the density correction. This is removed when the models
# are run without the density correction. However then the sqrt(std
# resid)~fitted plot shows that pattern and the plot of std_resid ~ air_temp
# shows the density error pattern really strongly. I think it is an issue of
# being an imperfect correction and so some error still shows through.

# ggplot(training_data,
#        aes(y = raw_water_level_cm - correction_cm - water_depth_cm,
#            x = air_temperature_c,
#            color = experiment)) +
#   geom_point() +
#   facet_wrap(~baro_sn + water_sn,
#              scales = "free")
# 
# ggplot(training_data,
#        aes(y = raw_water_level_cm - correction_cm - water_depth_cm,
#            x = water_temperature_c,
#            color = experiment)) +
#   geom_point() +
#   facet_wrap(~baro_sn + water_sn,
#              scales = "free")

ggplot(training_data[water_sn %in% c("1033239", "1062452", "2030899")],
       aes(x = sample_time,
           color = baro_sn)) +
  geom_ribbon(aes(ymin = water_depth_cm - instrument_error_cm * qnorm(0.995),
                  ymax = water_depth_cm + instrument_error_cm * qnorm(0.995)),
              fill = "#bfbfbf",
              color = NA) +
  geom_line(aes(y = centered_water_level_cm),
            linetype = "dotted") +
  geom_line(aes(y = corrected_water_level_cm)) +
  facet_wrap(~water_sn,
             scales = "free")

test_data[mods, 
          `:=`(intercept = i.intercept,
               slope = i.slope),
          on = c("baro_sn", "water_sn")]

test_data[, X := list(list(c(air_temperature_c, water_temperature_c, delta_at_01c_min))),
          by = .(water_sn, baro_sn, sample_time)]

test_data[, 
          correction_cm := mapply(function(i, b, x){i + x %*% b}, i = intercept, b = slope, x = X),
          by = .(water_sn, baro_sn, sample_time)]

test_data[, 
          `:=`(centered_water_level_cm = raw_water_level_cm - mean(raw_water_level_cm - water_depth_cm),
               corrected_water_level_cm = (raw_water_level_cm - correction_cm) - mean((raw_water_level_cm - correction_cm) - water_depth_cm)),
          by = .(water_sn, baro_sn, experiment)]

# ggplot(test_data,
#        aes(x = raw_error_cm,
#            y = correction_cm)) +
#   geom_point() +
#   facet_wrap(~baro_sn + water_sn,
#              scales = "free")

# Loggers that show a different slope sign for var-dis than for other
# experimental periods
# bad_loggers <- 
#   c("2064734", "1062528", "2030899")

ggplot(test_data[experiment == "var-dis"],
       aes(x = sample_time,
           color = baro_sn)) +
  geom_ribbon(aes(ymin = water_depth_cm - instrument_error_cm * qnorm(0.975),
                  ymax = water_depth_cm + instrument_error_cm * qnorm(0.975)),
              fill = "#bfbfbf",
              color = NA) +
  geom_line(aes(y = centered_water_level_cm),
            linetype = "dotted") +
  geom_line(aes(y = corrected_water_level_cm)) +
  facet_wrap(~water_sn + baro_sn,
             scales = "free")



figure_mad <- 
  test_data[,
            lapply(.SD, function(x) mad(x - water_depth_cm)),
            by = .(water_sn, baro_sn, experiment),
            .SDcols = patterns("ed_water_level_cm")] %>%
  melt(id.vars = c("water_sn", "baro_sn", "experiment")) %>%
  .[, id := paste(water_sn, baro_sn)] %>%
  .[, variable := fifelse(variable == "centered_water_level_cm", "Raw Error", "Temperature Corrected Error")] %>% 
  ggplot(data = .,
         aes(x = experiment,
             y = value,
             fill = variable)) +
  scale_fill_manual(name = NULL, values = pale_pal) +
  geom_boxplot() +
  theme(legend.position = c(0.2, 0.8))

ggsave("output/figures/boxplot_median_absolute_deviation_by_experiment.png",
       figure_mad)



library(brms)
options(mc.cores = 4)

# brm_dat <- 
#   compensated_data[water_sn %in% c("1066016", "2013939") & experiment == "var-dis"]

brm_form <- 
  bf(raw_error_cm | se(instrument_error_cm, sigma = TRUE) ~ 
       me(air_temperature_c, air_temp_error_c) + 
       me(water_temperature_c, water_temp_error_c) +
       delta_at_01c_min)

# Can't keep track of progress doing it this way
brm_dat <-
  split(training_data,
        by = c("water_sn", "baro_sn"))

names(brm_dat) <- 
  sub("\\.", "_", names(brm_dat))

# 
# brm_mods <- 
#   brm_multiple(formula = brm_form, 
#                data = brm_dat,
#                combine = FALSE)

base_brm <- 
  brm(brm_form, 
      data = training_data, 
      chains = 1,
      iter = 1)

priors <- 
  mods[, .(water_sn, baro_sn, rbindlist(lapply(slope, function(x){data.table(x)})))]

priors[,  c("air_temperature_c", "water_temperature_c", "delta_at_01c_min") :=
           lapply(.SD, 
                function(x){paste0("normal(", x, ",", round(sd(x), 2), ")")}),
       .SDcols = c("air_temperature_c", "water_temperature_c", "delta_at_01c_min")]

setkey(priors, "water_sn", "baro_sn")

brm_mods <-
  lapply(brm_dat,
         function(x){
           wsn <- 
             unique(x$water_sn)
           
           bsn <- 
             unique(x$baro_sn)
           
           model_id <- 
             paste0(wsn, "_", bsn)
           
           rds_file <- 
             paste0("output/models/brm/", 
                    model_id,
                    ".rds")
           
           stan_file <- 
             paste0("output/models/brm/", 
                    model_id,
                    ".stan")

           if(!file.exists(rds_file)){
             
             priors <- 
               c(set_prior(priors[CJ(wsn, bsn), air_temperature_c], 
                           class = "b", 
                           coef = "air_temperature_c"),
                 set_prior(priors[CJ(wsn, bsn), water_temperature_c],  
                           class = "b",
                           coef = "water_temperature_c"),
                 set_prior(priors[CJ(wsn, bsn), delta_at_01c_min],  
                           class = "b",
                           coef = "delta_at_01c_min"))
             
             fit_mod <- 
               update(base_brm,
                      newdata = x,
                      chains = 4,
                      iter = 2000,
                      prior = priors,
                      save_model = stan_file)
             
             saveRDS(fit_mod,
                     rds_file)
           }
           
           RPushbullet::pbPost("note",
                               title = paste0("BRM Model Fit ", 
                                              which(names(brm_dat) == model_id),
                                              "/", length(brm_dat)),
                               body = paste0("Water SN: ", wsn, "\nBaro SN: ", bsn))
         })

# fit_brm <- 
#   brm(brm_form,
#       data = vardis,
#       save_model = "full_model_brm")
# 
# saveRDS(fit_brm, "output/models/brm_hierachical_test_two_loggers.rds")

# Start sampling
# Warning messages:
#   1: There were 700 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 3300 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: The largest R-hat is 3.96, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#r-hat 
# 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 
# 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess 

fit_brm


library(rstan)

# rstan_options(auto_write = TRUE)


modcc <- 
  stanc("code/full_uncertainty.stan", 
        model_name = "uncertainty_mod")
 
mod_bin <- 
  stan_model("code/full_uncertainty.stan",
             model_name = "uncertainty_mod")

fit <- 
  sampling(mod_bin, 
           data = list(N = nrow(vardis), 
                       tau_x = vardis$air_temp_error_c, 
                       x_obs = vardis$air_temperature_c, 
                       y_obs = vardis$raw_error_cm,
                       tau_y = vardis$instrument_error_cm),
           chains = 3,
           pars = c("alpha", "beta", "sigma"))

saveRDS(fit, "output/models/stan_test.rds")

# fit <- 
#   stan("code/full_uncertainty.stan", 
#        data = list(N = nrow(dat), 
#                    tau = dat$air_temp_error_c, 
#                    x_meas = dat$air_temperature_c, 
#                    y = dat$raw_error_cm,
#                    tau_y = dat$instrument_error_cm))

fit <- readRDS("output/models/stan_test.rds")

print(fit, pars = c("alpha", "beta"))

# From:
# https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed

ex_fit <- 
  extract(fit)

alpha_post <- ex_fit$alpha
beta_post <- ex_fit$beta
sigma_post <- ex_fit$sigma

# Function for simulating y based on new x
predict_post <- function(x, alpha = 0.05) {
  set.seed(1234)
  alpha <- matrix(rep(sample(alpha_post, 100, replace = TRUE),
                      times = length(x)), 
                  nrow = length(x))
  beta <- matrix(rep(sample(beta_post, 100, replace = TRUE),
                     times = length(x)), 
                 nrow = length(x))
  mat <- 
    alpha + beta * x
  
  cbind(estimate = apply(mat, 1, median),
        lwr = apply(mat, 1, quantile, probs = 0.025),
        upr = apply(mat, 1, quantile, probs = 0.975))
}

trad_model <- 
  lm(raw_error_cm ~ air_temperature_c, data = vardis)


pred_trad <- 
  as.data.frame(cbind(air_temperature_c = seq(min(vardis$air_temperature_c),
                                        max(vardis$air_temperature_c),
                                        by = 0.5),
                predict(trad_model,
                        newdata = data.frame(air_temperature_c = seq(min(vardis$air_temperature_c),
                                                                     max(vardis$air_temperature_c),
                                                                     by = 0.5)),
                        interval = "prediction")))


# Run the function on x_test
plot(x = vardis$air_temperature_c,
     y = vardis$raw_error_cm,
     col = rgb(.25, .25, .25, .4))
points(x = pred_trad$air_temperature_c,
       y = pred_trad$fit, 
       type = "l")
points(x = vardis$air_temperature_c,
       y = pred_trad$lwr,
       type = "l", 
       lty = "dashed")
points(x = vardis$air_temperature_c,
       y = predict(trad_model, interval = "prediction")[,"upr"],
       type = "l", 
       col = "green",
       lty = "dashed")
points(x = vardis$air_temperature_c,
       y = predict_post(vardis$air_temperature_c)[, "estimate"], 
       col = "blue",
       pch = 20)
segments(x0 = vardis$air_temperature_c,
         y0 = predict_post(vardis$air_temperature_c)[, "lwr"],
         y1 = predict_post(vardis$air_temperature_c)[, "upr"],
         col = "red")


