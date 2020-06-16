#' ---
#' title: "Research Proposal Defense"
#' subtitle: "PhD, Forest Science"
#' author: "Joe Shannon"
#' institute: "College of Forest Resources and Environmental Science,</br>Michigan Technological University"
#' date: "January 28, 2020"
#' output:
#'   xaringan::moon_reader:
#'   css: [shannon.css, xaringan-themer.css]
#' self_contained: false
#' lib_dir: libs
#' nature:
#'   highlightStyle: github
#' highlightLines: true
#' countIncrementalSlides: false
#' ---

library(data.table)
library(ingestr)
library(ggplot2)

water_files <- 
  list.files(pattern = "[0-9].xle")

baro_file <- 
  "solinst_temperature_test_baro.xle"

water <- 
  lapply(water_files,
         function(x){data.table(merge(ingest_xle(x), ingest_header(x), by = "input_source"))})

water <- 
  rbindlist(water, fill = TRUE)

water <- 
  water[, .(serial_number, 
            sample_time = as.POSIXct(paste(sample_date, sample_time), tz = "EST5EDT"),
            water_level_cm = ifelse(level_cm < 1000, level_cm + 960, level_cm),
            water_temp_c = pmax(temperature_c, temperature_deg_c, na.rm = TRUE))]

setkey(water, sample_time, serial_number)

baro <- 
  data.table(merge(ingest_xle(baro_file), ingest_header(baro_file), by = "input_source"))

baro <- 
  baro[, .(sample_time = as.POSIXct(paste(sample_date, sample_time), tz = "EST5EDT"),
           baro_level_cm = level_cm,
           air_temp_c = temperature_c)]

setkey(baro, sample_time)

combined <- 
  water[baro][sample_time <= as.POSIXct("2020-05-08 13:28", tz = "EST5EDT")]

# Adjust to equal temperatures
combined[, compensated_level_cm := water_level_cm - baro_level_cm]
combined[, adj_compensated_level_cm := water_level_cm - baro_level_cm * ((air_temp_c + 273.15) / (water_temp_c + 273.15))]
combined[, `:=`(delta_wl = c(NA_real_, diff(compensated_level_cm)),
                delta_wt = c(NA_real_, diff(water_temp_c))),
         by = .(serial_number)]

combined[, `:=`(delta_wl = nafill(delta_wl, "nocb"),
                delta_wt = nafill(delta_wt, "nocb")),
         by = .(serial_number)]

combined[, seconds_elapsed := as.numeric(sample_time - min(sample_time)),
         by = .(serial_number)]

ggplot(combined, 
       aes(x = sample_time, 
           y = adj_compensated_level_cm,
           color = serial_number)) +
  geom_line()

par(mar = c(5.1, 4.1, 2.1, 4.1))
plot(water_temp_c ~ sample_time,
     data = subset(combined, serial_number == "1065861"),
     main = "1065861", 
     type = "l",
     yaxt = "n",
     xaxt = "n",
     ylab = "",
     xlab = "",
     col = "gray60",
     lwd = 3)
axis(side = 4)
mtext("Water Temperature (C)", 
      side = 4,
      line = 3)
par(new = TRUE)
plot(adj_compensated_level_cm ~ sample_time,
     data = subset(combined, serial_number == "1065861"),
     type = "l",
     col = "#0072B2",
     ylab = "Compensated Water Level (cm)",
     xlab = "Sample Time")

plot(water_temp_c ~ sample_time,
     data = subset(combined, serial_number == "1066019"),
     main = "1066019",
     type = "l",
     yaxt = "n",
     xaxt = "n",
     ylab = "",
     xlab = "",
     col = "gray60",
     lwd = 3)
axis(side = 4)
mtext("Water Temperature (C)", 
      side = 4,
      line = 3)
par(new = TRUE)
plot(adj_compensated_level_cm ~ sample_time,
     data = subset(combined, serial_number == "1066019"),
     type = "l",
     col = "#0072B2",
     ylab = "Compensated Water Level (cm)",
     xlab = "Sample Time")


plot(water_temp_c ~ sample_time,
     data = subset(combined, serial_number == "2013939"),
     main = "2013939",
     type = "l",
     yaxt = "n",
     xaxt = "n",
     ylab = "",
     xlab = "",
     col = "gray60",
     lwd = 3)
axis(side = 4)
mtext("Water Temperature (C)", 
      side = 4,
      line = 3)
par(new = TRUE)
plot(adj_compensated_level_cm ~ sample_time,
     data = subset(combined, serial_number == "2013939"),
     type = "l",
     col = "#0072B2",
     ylab = "Compensated Water Level (cm)",
     xlab = "Sample Time")

ggplot(combined, 
       aes(x = water_temp_c, 
           y = adj_compensated_level_cm,
           color = delta_wt)) +
  geom_point() + 
  scale_color_gradientn(colors = blue_orange_scale) +
  # scale_color_distiller(type = "div") +
  # scale_color_binned(type = "viridis") +
  facet_wrap(~serial_number)

ggplot(combined, 
       aes(x = water_temp_c, 
           y = adj_compensated_level_cm,
           color = ifelse(abs(delta_wt) < 0.75, "stable", "changing"))) +
  geom_point() + 
  facet_wrap(~serial_number)

ggplot(combined, 
       aes(x = water_temp_c, 
           y = adj_compensated_level_cm,
           color = seconds_elapsed)) +
  geom_point() + 
  facet_wrap(~serial_number)


ggplot(combined, 
       aes(x = delta_wt, 
           y = delta_wl)) +
  geom_point() + 
  facet_wrap(~serial_number)
