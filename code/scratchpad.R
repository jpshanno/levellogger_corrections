# Supplemental Table 1

loadd(bootstrap_models)

summary_func <- function(x) {
  m <- mean(x)
  s <- sd(x)
  n <- length(x)
  sprintf("%.3f ± %.3f", m, s/sqrt(n))
}

model_summary <- bootstrap_models[
  experiment == "test-dat", 
  c(
    lapply(.SD, summary_func),
    rmse = round(mean(rmse), 3)
  ),
  by = .(water_sn, baro_sn),
  .SDcols = c("intercept", "air_temperature_c", 
              "water_temperature_c", "delta_at_01c_min")]


write.csv(model_summary, "output/tabular/Supplemental_Table_S1-Model_Fit_Summaries.csv", row.names = FALSE)


source("code/levellogger_packages.R")
source("code/levellogger_functions.R")
theme_set(ggthemes::theme_few())
brown_green_scale <- # Inspired by NOAA maps
  c('#97601c', '#a47138', '#b08353', '#bc956e', '#c6a78a', '#d0baa6', '#d8cdc3', '#e0e0e0', '#c6d1cf', '#adc2bf', '#93b3af', '#7aa49f', '#5f9690', '#438781', '#207972')

blue_orange_scale <-
  c('#3085c7', '#4a98d4', '#6babde', '#90bee3', '#b7cfe3', '#e0e0e0', '#eac79d', '#eaae66', '#e5943a', '#dd7a15', '#d55e00')

pale_pal <- 
  c(green = "#7DA050",
    orange = "#D19648",
    teal = "#329985",
    blue = "#5982A0",
    purple = "#966283",
    red = "#9B5249")

loadd(case_study)
loadd(water_balance)
loadd(daily_water_balance)

error_tests <- 
  water_balance[, .(sample_time, 
                    instrument_error_cm,
                    raw_compensated_level_cm,
                    corrected_compensated_level_cm)]

error_tests[, .N, 
                  by = .(error = ifelse(abs(raw_compensated_level_cm - corrected_compensated_level_cm) > instrument_error_cm,
                                          "significant",
                                          "insignificant"))]

error_tests[, error_instant := ifelse(abs(raw_compensated_level_cm - corrected_compensated_level_cm) > instrument_error_cm,
                                        "significant",
                                        "insignificant")]

error_tests[, error_hourly := ifelse(abs(frollmean(raw_compensated_level_cm - corrected_compensated_level_cm, 4, align = 'center')) > instrument_error_cm,
                                           "significant",
                                           "insignificant")]

error_tests[, error_daily := ifelse(abs(frollmean(raw_compensated_level_cm - corrected_compensated_level_cm, 96, align = 'center')) > instrument_error_cm,
                                           "significant",
                                           "insignificant")]

error_tests[, error_3day := ifelse(abs(frollmean(raw_compensated_level_cm - corrected_compensated_level_cm, 288, align = 'center')) > instrument_error_cm,
                                          "significant",
                                          "insignificant")]

error_tests[, error_weekly := ifelse(abs(frollmean(raw_compensated_level_cm - corrected_compensated_level_cm, 672, align = 'center')) > instrument_error_cm,
                                   "significant",
                                   "insignificant")]

error_tests[, error_bimonthly := ifelse(abs(frollmean(raw_compensated_level_cm - corrected_compensated_level_cm, 1440, align = 'center')) > instrument_error_cm,
                                     "significant",
                                     "insignificant")]

error_tests[, error_monthly := ifelse(abs(frollmean(raw_compensated_level_cm - corrected_compensated_level_cm, 2880, align = 'center')) > instrument_error_cm,
                                        "significant",
                                        "insignificant")]

table(error_tests$error_instant)
table(error_tests$error_hourly)
table(error_tests$error_daily)
table(error_tests$error_3day)
table(error_tests$error_weekly)
table(error_tests$error_bimonthly)
table(error_tests$error_monthly)

# Instantaneous 1
# Hourly 4
# Daily 96
# 3 day 288
# Weekly 672
# Bimonthly 1440
# Monthly 2880
# Seasonally (3 month) 8736
# Annually .N

# Compare to PET with New ESY ---------------------------------------------
daily_water_balance[, pet_cm_d := -pet_cm_d]

daily_water_balance[, `:=`(external_corrected_compensated_level_cm = corrected_compensated_level_cm,
                external_raw_compensated_level_cm = raw_compensated_level_cm)]

daily_water_balance[(dry_day)] %>% 
  {plot(raw_et_cm_d ~ pet_cm_d, data = ., col = 'red')
    points(corrected_et_cm_d ~ -pet_cm_d, data = ., col = 'blue')
    abline(0, 1)}

base_font <- 
  14

labels <- 
  data.frame(lab = c(glue("{`(Intercept)`} + {pet_cm_d}&times;PET<br>R<sup>2</sup> = {r2}", 
                          .envir = c(as.list(round(coef(corrected_mod), 2)), 
                                     list(r2 = round(corrected_summ$r.squared, 2)))),
                     glue("{`(Intercept)`} + {pet_cm_d}&times;PET<br>R<sup>2</sup> = {r2}", 
                          .envir = c(as.list(round(coef(raw_mod), 2)), 
                                     list(r2 = round(raw_summ$r.squared, 2))))),
             type = factor(c("Corrected", "Raw"),
                           levels = c("Raw", "Corrected"), 
                           labels = c("Raw", "Corrected"),
                           ordered = TRUE),
             x = c(0.21, 0.45),
             y = c(0.14, 0.28),
             angle = c(20, 7),
             lab11 = c("1:1", "1:1"))

melt(daily_water_balance[(dry_day)],
     id.vars = c("sample_date", "pet_cm_d", "total_input_cm_d"), 
     measure.vars = patterns("(_et_cm_d|compensated_level)"), 
     variable.name = "type", 
     value.name = "value") %>% 
  transform(variable = str_remove(type, "^(external_)?(corrected|raw)_"), 
            type = str_extract(type, "^(external_)?(corrected|raw)")) %>% 
  dcast(... ~ variable, value.var = "value") %>% 
  subset(et_cm_d > 0 & str_detect(type, "external", negate = TRUE)) %>% 
  transform(type = factor(type, 
                          levels = c("raw", "corrected"), 
                          labels = c("Raw", "Corrected"),
                          ordered = TRUE)) %>% 
  ggplot() +
  aes(x = pet_cm_d, 
      y = et_cm_d) +
  geom_point() + 
  geom_abline(linetype = 'dashed') +
  geom_richtext(data = labels,
                aes(x = -Inf, y = Inf, label = lab),
                hjust = 0,
                vjust = 1,
                fill = "#FFFFFF80",
                size = 0.9 * base_font * 5/14,
                label.margin = unit(c(1.5, 1), 'lines'),
                label.colour = NA) +
  geom_text(data = labels,
            aes(x = x,
                y = y, 
                label = lab11,
                angle = angle),
            size = 0.9 * base_font * 5/14) +
  geom_smooth(method = lm,
              formula = y ~ x,
              color = 'blue',
              se = FALSE) +
  labs(x = "PET (cm d<sup>-1</sup>)",
       y = "Evapotranspiration (cm d<sup>-1</sup>)") +
  facet_wrap(~type,
             scales = 'free') +
  theme_minimal(base_size = base_font) +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        strip.text = element_text(size = 0.9*base_font,
                                  face = 'bold'))

cor(daily_water_balance[(dry_day), .(pet_cm_d, raw_et_cm_d, corrected_et_cm_d)],
    use = "pairwise.complete",
    method = 'pearson')

# Tested glm with Gamma(identity) family. Did not improve model
raw_mod <- 
  lm(raw_et_cm_d ~ pet_cm_d, 
     data = daily_water_balance[(dry_day)]) 
raw_summ <- 
  summary(raw_mod)


corrected_mod <- 
  lm(corrected_et_cm_d ~ pet_cm_d, 
     data = daily_water_balance[(dry_day)]) 
corrected_summ <- 
  summary(corrected_mod)

paste(c(deparse(bquote(.(`(Intercept)`) + .(pet_cm_d)*PET,
                       where = as.list(round(coef(corrected_mod), 2)))),
        paste("R^2^ =", round(corrected_summ$r.squared, 2))), 
collapse = "\n")


plot(corrected_et_cm_d ~ raw_et_cm_d, 
     data = daily_water_balance[total_input_cm_d < 0.1 & shift(total_input_cm_d, 1) < 0.1 & shift(total_input_cm_d, 2) < 0.1])


melt(daily_water_balance[(dry_day)],
     id.vars = c("sample_date", "pet_cm_d"), 
     measure.vars = patterns("_et_cm_d"), 
     variable.name = "type", 
     value.name = "et_hat_cm_d") %>% 
  transform(type = str_extract(type, "^(external_)?[a-z]+")) %>% 
  subset(str_detect(type, "external", negate = TRUE)) %>% 
  ggplot() +
  aes(x = (et_hat_cm_d - pet_cm_d) / pet_cm_d, 
      y = type) +
  ggridges::geom_density_ridges(rel_min_height = 0.01,
                                color = 'white') +
  theme_minimal()

melt(daily_water_balance[(dry_day)],
     id.vars = c("sample_date", "pet_cm_d"), 
     measure.vars = patterns("_et_cm_d"), 
     variable.name = "type", 
     value.name = "et_hat_cm_d") %>% 
  transform(type = str_extract(type, "^(external_)?[a-z]+")) %>% 
  subset(str_detect(type, "external", negate = TRUE)) %>% 
  ggplot() +
  aes(x = type, 
      y = (et_hat_cm_d - pet_cm_d)) +
  ggdist::stat_halfeye(width = 0.4) +
  theme_minimal()

dry_data <- 
  melt(daily_water_balance[total_input_cm_d < 0.1 & shift(total_input_cm_d, 1) < 0.1],
       measure.vars = patterns("_et_cm_d"), 
       id.vars = c("sample_date", "pet_cm_d"),
       variable.name = 'type',
       value.name = 'et_hat_cm_d')[, .(dat = list(.SD)), by = .(type)]

dry_data[, mod := map(dat,
                      ~brm(pet_cm_d ~ et_hat_cm_d,
                           data = .x,
                           family = Gamma(identity),
                           prior = c(prior(normal(1, 1), class = b, lb = 0),
                                     prior(normal(0, 0.1), class = Intercept))))]

dry_data[, mod := map(mod,
                      add_criterion,
                      "loo")]

loo_compare(dry_data$mod[[1]], dry_data$mod[[2]])

loos <- 
  rbind(data.table(type = 'raw',
                   pet_cm_d = dry_data$mod[[1]]$data$pet_cm_d, 
                   loo = dry_data$mod[[1]]$criteria$loo$pointwise[, "elpd_loo"]),
        data.table(type = 'corrected',
                   pet_cm_d = dry_data$mod[[2]]$data$pet_cm_d,
                   loo = dry_data$mod[[2]]$criteria$loo$pointwise[, "elpd_loo"]))

ggplot(loos) +
  aes(x = pet_cm_d,
      y = loo,
      color = type) +
  geom_point()

# Posterior predictive plots
# Look into using tidybayes for doing these

{ggplot(data.table(obs = rep(dry_data$dat[[1]]$et_hat_cm_d, times = 100), sample_id = rep(1:100, each = 70), pred = as.numeric(posterior_predict(dry_data$mod[[1]], nsamples = 1000)))) +
  geom_density(aes(x = pred, group = sample_id), alpha = 0.001, col = 'gray80') +
  geom_density(aes(x = obs), color = 'red') +
  ggtitle("Raw") +
    theme_bw()} + 
  {ggplot(data.table(obs = rep(dry_data$dat[[2]]$et_hat_cm_d, times = 100), sample_id = rep(1:100, each = 70), pred = as.numeric(posterior_predict(dry_data$mod[[2]], nsamples = 1000)))) +
      geom_density(aes(x = pred, group = sample_id), alpha = 0.001, col = 'gray80') +
      geom_density(aes(x = obs), color = 'red') +
      ggtitle("Corrected") +
      theme_bw()}


# Esy Comparison
plot(raw_compensated_level_cm ~ ytd_water_balance, data = daily_water_balance[(daily_water_balance[1:which.min(raw_compensated_level_cm), which.max(raw_compensated_level_cm)]):which.min(raw_compensated_level_cm)], type = 'l')
lines(corrected_compensated_level_cm ~ ytd_water_balance, data = daily_water_balance[(daily_water_balance[1:which.min(corrected_compensated_level_cm), which.max(corrected_compensated_level_cm)]):which.min(corrected_compensated_level_cm)], col = 'red')




# Low R-squared not suprising:
# @zhu-2011 showed seasonal oscillations masking daily patterns
# @watras-2017 couldn't fit ESy function to pond-adjacent wetlands
# @lafluer-2005 direct measurement of ET via eddy flux & r^2 as low as 0.56 for ET ~ PET
# @wang-2014 different method but r^2 ~ 0.25
# @mclaughlin-2014 doesn't report ET ~ PET r^2, but is a noise relationship in their
# ET/PET, suggesting low r^2
# @soylu-2012 shows White method with r^2 <= 0.2 & proposed method with 0.33 <= r^2 <= 0.4 



# Temperature Ranges ------------------------------------------------------

daily_water_balance[water_balance[, map(.SD, value_range), 
                                  by = .(sample_date),
                                  .SDcols = patterns("temperature_c")],
                    on = "sample_date",
                    `:=`(water_temp_range = i.water_temperature_c,
                         air_temp_range = i.air_temperature_c)]

daily_water_balance[water_balance[, map(list("min", "mean", "max"),
                                        ~exec(.x, (air_temperature_c - water_temperature_c))), 
                                  by = .(sample_date)],
                    on = "sample_date",
                    `:=`(temp_diff_min = i.V1,
                         temp_diff_mean = i.V2,
                         temp_diff_max = i.V3)]

daily_water_balance[water_balance[, map(.SD, mean, na.rm = TRUE), 
                                  by = .(sample_date),
                                  .SDcols = patterns("temperature_c")],
                    on = "sample_date",
                    `:=`(water_temperature_c = i.water_temperature_c,
                         air_temperature_c = i.air_temperature_c)]

daily_water_balance[, daily_means_error_cm := raw_compensated_level_cm - corrected_compensated_level_cm]

ggplot(daily_water_balance) +
  aes(x = sample_date,
      y = air_temperature_c) +
  geom_point()