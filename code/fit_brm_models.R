library(brms)
source("code/levellogger_packages.R")
source("code/levellogger_functions.R")
drake::loadd(compensated_data)
options(mc.cores = 4)

training_data <- 
  compensated_data[experiment == "test-dat"]

training_data[, delta_at_error_c_min := combine_errors(air_temp_error_c, air_temp_error_c)]

mods <- 
  fit_correction_mod(training_data,
                     y = "raw_error_cm",
                     x = c("air_temperature_c", 
                           "water_temperature_c", 
                           "delta_at_c_min"),
                     by = c("baro_sn", "water_sn"))

brm_form <- 
  bf(raw_error_cm | se(instrument_error_cm, sigma = TRUE) ~ 
       me(air_temperature_c, air_temp_error_c) + 
       me(water_temperature_c, water_temp_error_c) +
       me(delta_at_c_min, delta_at_error_c_min))

brm_dat <-
  split(training_data,
        by = c("water_sn", "baro_sn"))

names(brm_dat) <- 
  sub("\\.", "_", names(brm_dat))

priors <- 
  mods[, .(water_sn, baro_sn, mod)]

priors[, mod_summary := lapply(mod, tidy)]
priors[, mod_fit := lapply(mod, glance)]

priors <- 
  priors[, c(list(water_sn = water_sn,
                  baro_sn = baro_sn,
                  rbindlist(lapply(mod_summary,
                                   function(x){
                                     as.list(setNames(paste0("normal(",
                                                             x$estimate, ",",
                                                             x$std.error, ")"), 
                                                      x$term))})),
                  sigma = vapply(mod_fit, `[[`, numeric(1), "sigma")))]
setnames(priors, "X.Intercept.", "intercept")

priors[, sigma := paste0("normal(", sigma, ",", sd(sigma),")")]

setkey(priors, "water_sn", "baro_sn")

# Need to set priors and update with stanvars to avoid recompile
# https://github.com/paul-buerkner/brms/issues/219

brms_dir <- 
  "output/models/brms/"

base_brm <- 
  brm(brm_form, 
      data = brm_dat[[1]], 
      chains = 1,
      iter = 1,
      save_model = paste0(brms_dir, "base_model.stan"))


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
             paste0(brms_dir, 
                    "brm_",
                    model_id,
                    ".rds")
           
           stan_file <- 
             paste0(brms_dir,
                    "brm_",
                    model_id,
                    ".stan")
           
           if(!file.exists(rds_file)){
             
             new_priors <-
               c(set_prior(priors[CJ(wsn, bsn), intercept],
                           class = "b",
                           coef = "Intercept"),
                 set_prior(priors[CJ(wsn, bsn), air_temperature_c],
                           class = "b",
                           coef = "meair_temperature_cair_temp_error_c"),
                 set_prior(priors[CJ(wsn, bsn), water_temperature_c],
                           class = "b",
                           coef = "mewater_temperature_cwater_temp_error_c"),
                 set_prior(priors[CJ(wsn, bsn), delta_at_c_min],
                           class = "b",
                           coef = "medelta_at_c_mindelta_at_error_c_min"),
                 set_prior(priors[CJ(wsn, bsn), sigma],
                           class = "Intercept"))
             
             # new_priors <-
             #   stanvar(priors[CJ(wsn, bsn), air_temperature_c],
             #           name = "meair_temperature_cair_temp_error_c",
             #           block = "parameters") +
             #   stanvar(priors[CJ(wsn, bsn), water_temperature_c],
             #           name = "mewater_temperature_cwater_temp_error_c",
             #           block = "parameters") +
             #   stanvar(priors[CJ(wsn, bsn), delta_at_c_min],
             #           name = "medelta_at_c_mindelta_at_error_c_min",
             #           block = "parameters")
             
             fit_mod <- 
               update(base_brm,
                      newdata = x,
                      chains = 4,
                      iter = 2000,
                      prior = new_priors,
                      # stanvars = new_priors,
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


