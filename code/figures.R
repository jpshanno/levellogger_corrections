source("code/levellogger_functions.R")
source("code/levellogger_packages.R")
library(ggplot2)
library(patchwork)
loadd(combined_data,
      bootstrap_models,
      predicted,
      fitted)

theme_set(ggthemes::theme_few())
brown_green_scale <- # Inspired by NOAA maps
  c('#97601c', '#a47138', '#b08353', '#bc956e', '#c6a78a', '#d0baa6', '#d8cdc3', '#e0e0e0', '#c6d1cf', '#adc2bf', '#93b3af', '#7aa49f', '#5f9690', '#438781', '#207972')

blue_orange_scale <-
  c('#3085c7', '#4a98d4', '#6babde', '#90bee3', '#b7cfe3', '#e0e0e0', '#eac79d', '#eaae66', '#e5943a', '#dd7a15', '#d55e00')

pale_pal <- 
  c(green = "#7DA050",
    orange = "#D19648",
    purple = "#966283",
    teal = "#329985",
    blue = "#5982A0",
    red = "#9B5249")

options(ggplot2.discrete.colour = unname(pale_pal),
        ggplot2.discrete.fill = unname(pale_pal))


# Figure AAA Three panels showing error ~ air temperature, error ~ water
# temperature, and air-temperature-corrected error ~ water temperature

dat_1 <- 
  combined_data$testing[experiment == "var-dis" & baro_sn == "1066019" & water_sn == "1062452"]

fig1a <- 
  ggplot(dat_1,
         aes(x = air_temperature_c,
             y = raw_error_cm)) + 
  geom_point(color = pale_pal["teal"]) +
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
  geom_point(color = pale_pal["teal"]) +
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
  geom_point(color = pale_pal["teal"]) +
  geom_smooth(method = "lm",
              formula = "y~x",
              se = FALSE,
              color = "black") +
  labs(x = expression(paste("Water Temperature, ", degree, "C")),
       y = "Raw Error minus Air Temperature Effect, cm")


fig1 <- 
  fig1a + fig1b + fig1c + plot_annotation(tag_levels = "A")

# Figure BBB: A line graph showing the predictions from about ~25 separate
# bootstrap models along with the final median prediction and uncertainty
# derived from all 1000 models

dat2 <- 
  fitted[experiment == "var-dis" & baro_sn == "1066019" & water_sn == "1062452"]

mods <- 
  read_fst("output/tabular/bootstrap_fit_values.df/1062452_1066019_var-dis.fst",
           as.data.table = TRUE)

r2_quartiles <- 
  quantile(bootstrap_models[water_sn == "1062452" &
                              baro_sn == "1066019" & 
                              experiment == "var-dis", adj_r2], 
           seq(0, 1, 0.5),
           type = 1)

quartile_reps <- 
  bootstrap_models[water_sn == "1062452" & 
                     baro_sn == "1066019" & 
                     experiment == "var-dis" & 
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
  labs(x = "Water Level", 
       y = "Sample Time") +
  facet_zoom(x = sample == "yes",
             zoom.size = 0.5, 
             horizontal = FALSE,
             ylim = range(mods2[sample == "yes", rect_water_level_cm])) + 
  theme(strip.background = element_rect(color = "black",
                                        size = rel(1.5)))


# Figure DDD: Plots showing the distribution of bootstrap model coefficients

dat4a <- 
  bootstrap_models[experiment == "test-dat", 
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
  labs(y = expression(beta[Air~Temperature]),
       x = "Barometric Pressure Serial Number") +
  scale_shape_manual(values = c(19, 21))

dat4b <- 
  bootstrap_models[experiment == "test-dat", 
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
                  show.legend = FALSE) +
  labs(y = expression(beta[Water~Temperature]),
       x = "Water Pressure Serial Number") +
  scale_shape_manual(values = c(19, 21)) +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1))

fig4a + fig4b

# Figure EEE: A panel of two plots. The first showing time series of corrected
# water level to highlight where the error occurs. The second would be a scatter
# plot showing residual error as a function of temperature

# Figure FFF: Point-range plots of published water level fluctuation rates with
# various uncertainty intervals from this study added as thresholds for
# comparison

fread(text = 
"source; signal_magnitude; standard_error; notes
Cuevas, 2010; ; ;
Carlson Mazur, 2014; ; ; thermally dissimilar transducers
diamond; ;
white; ;
mclaughlin; ;
kirchner; ;
Loheide, 2005; ;",
sep = ";")

data.table(source = c("carlson mazur", "diamond", "white", "mclaughlin", "kirchner", "Loheide, 2005"),
           signal_magnitude = ,
           standard_error = )
# 
# Mazur
# Diamond
# White
# McLaughlin
# Kirchner?
  
  

# Bootstrap vs OLS Coefs --------------------------------------------------

ols_models <- 
  fit_correction_mod(combined_data$training, 
                     x = c("air_temperature_c", "water_temperature_c", 
                           "delta_at_01c_min"),
                     y = "raw_error_cm", 
                     by = c("water_sn", "baro_sn"))

ols_coefs <- 
  ols_models[, map2(intercept, slope,
                    ~data.table(intercept = .x, .y))[[1]],
             by = .(water_sn, baro_sn)]


compare_coefs <- 
  melt(bootstrap_models[str_detect(model_id, "test-dat")], 
       id.vars = c("model_id", "water_sn", "baro_sn", "model_rep"), 
       measure.vars = c("intercept", "air_temperature_c", 
                        "water_temperature_c", "delta_at_01c_min"), 
       variable.name = "beta", 
       value.name = "boot_value")[melt(ols_coefs, 
                                       id.vars = c("water_sn", "baro_sn"), 
                                       variable.name = "beta", 
                                       value.name = "ols_value"), 
                                  on = c("water_sn", "baro_sn", "beta")]  

compare_labels <- 
  compare_coefs[beta == "intercept", 
                .(beta = unique(beta),
                  xend = quantile(boot_value, 0.75),
                  yend = 0.75,
                  x = quantile(boot_value, 0.9), 
                  y = 0.95, 
                  text = "Bootstrapped\nEstimates"),
                by = .(water_sn, baro_sn)]

ols_labels <- 
  compare_coefs[beta == "intercept", 
                .(beta = unique(beta),
                  xend = mean(ols_value),
                  yend = 0.75,
                  x = mean(ols_value) + 0.05 * diff(range(boot_value)), 
                  y = 0.75, 
                  text = "OLS\nEstimates"),
                by = .(water_sn, baro_sn)]

compare_coef_figures <- 
  split(compare_coefs,
        by = c("water_sn", "baro_sn")) %>% 
  map(~{ggplot(data = .x,
               aes(x = boot_value)) +
      geom_histogram(aes(y = ..ncount..), 
                     bins = 100,
                     fill = pale_pal[["teal"]],
                     alpha = 0.5) +
      geom_segment(aes(x = ols_value, 
                       xend = ols_value, 
                       y = 0, yend = 1),
                   color = pale_pal[["orange"]],
                   size = 1.5) +
      geom_text(data = compare_labels[.x[1, .(water_sn, baro_sn)], 
                                      on = c("water_sn", "baro_sn")],
                aes(x = x, y = y, label = text),
                vjust = 0.5,
                hjust = 0) +
      geom_curve(data = compare_labels[.x[1, .(water_sn, baro_sn)], 
                                       on = c("water_sn", "baro_sn")],
                 aes(x = x, xend = xend,
                     y = y, yend = yend),
                 arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
      geom_text(data = ols_labels[.x[1, .(water_sn, baro_sn)], 
                                  on = c("water_sn", "baro_sn")],
                aes(x = x, y = y, label = text),
                vjust = 0.5,
                hjust = 0) +
      geom_segment(data = ols_labels[.x[1, .(water_sn, baro_sn)], 
                                     on = c("water_sn", "baro_sn")],
                   aes(x = x, xend = xend,
                       y = y, yend = yend),
                   arrow = arrow(type = "closed", 
                                 length = unit(0.1, "inches"))) +
      geom_vline(data = .x[, .(boot_value = median(boot_value)), by = .(beta)],
                 aes(xintercept = boot_value)) +
      facet_wrap(~beta,
                 labeller = function(x){label_parsed(data.frame(str_replace_all(x[[1]], c("intercept" = "hat(beta)[0]",
                                                                                          "air_temperature_c" = "hat(beta)[1]",
                                                                                          "water_temperature_c" = "hat(beta)[2]",
                                                                                          "delta_at_01c_min" = "hat(beta)[3]"))))}, 
                 scales = "free_x") +
      labs(x = "Coefficient Estimate",
           y = "Proportion of Samples")})

compare_coef_figures[[2]]

# Figure GGG
