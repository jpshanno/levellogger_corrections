library(data.table)
library(lme4)

round_time <- 
  function(x){
    x_lt <- 
      as.POSIXlt(x)
    
    x_lt$min <- 
      (1 + x_lt$min %/% 15) * 15
    
    as.POSIXct(x_lt)
  }

# Create air density adjustment
t <- seq(-40, 50, by = 0.1)
correction <- 
  (9.80665 * 1000 * 10) / ((t + 273.15) * 287.058)

# plot(rho ~ t)
lm(correction ~ t)

# Not using HOBO data because there is only a small range that has temperature
# and pressure
# hobo_baro_file <- 
#   "../Raw/Hobo/csv_met/__combined_barometric_pressure_data.csv"

mesowest_file <- 
  "data/__supplementary_barometric_data.csv"

solinst_file <- 
  "data/__combined_uncorrected_solinst_baro_data.csv"

dat <- 
  lapply(c(mesowest_file, solinst_file),
         fread,
         select = c(site = "character", 
                    sample_time = "character", 
                    air_temperature_c = "numeric", 
                    sea_level_pressure_m = "numeric", 
                    input_source = "character"))

dat <- 
  rbindlist(dat)

dat[, sample_time := as.POSIXct(sample_time, tz = "EST")]

serial_numbers <- 
  unique(fread(solinst_file, 
                 select = c(site = "character", serial_number = "character")))

# It looks like the Weather Underground uses more than just altimeter to adjust
# pressure, making it equivalent to "sea level pressure" not "altimeter pressure"

# wu_raw <-
#   fread("../Paired_Watersheds/Data_Weather_Underground.csv")[
#     station == "KMIBERGL2" & pressure > -100 & temp > -999]
# 
# wu_raw[, sample_time := as.POSIXct(date, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")]
# wu_raw[, sea_level_pressure_m := pressure * 0.34540396841]
# wu_raw[, air_temperature_c := 5 / 9 * (temp - 32)]
# 
# attributes(wu_raw$sample_time)$tzone <- "EST"
# 
# wu_raw[, sample_time := round_time(sample_time)]
# 
# wu <-
#   wu_raw[, .(site = "KMIBERGL2",
#              sea_level_pressure_m = mean(sea_level_pressure_m),
#              air_temperature_c = mean(air_temperature_c),
#              input_source = "Data_Weather_Underground.csv"),
#          by = .(sample_time)][, .(site, sample_time, air_temperature_c, sea_level_pressure_m, input_source)]
# 
# dat <-
#   rbindlist(list(dat, wu))

# Compare Datasets --------------------------------------------------------

# Get all matched pressure observations in a wide format

diff_dat <- 
  merge(dat[site %in% c("053", "113") & !is.na(sea_level_pressure_m), 
            .(sample_time, 
              site, 
              sea_level_pressure_m, 
              air_temperature_c)],
        dat[!(site %in% c("053", "113")) & !is.na(sea_level_pressure_m), 
            .(sample_time, 
              station = site, 
              ex_sea_level_pressure_m = sea_level_pressure_m, 
              ex_air_temperature_c = air_temperature_c)],
        by = "sample_time", 
        all= FALSE)

setkey(diff_dat, sample_time, site, station)

diff_dat[, pressure_diff := ex_sea_level_pressure_m - sea_level_pressure_m]
diff_dat[, station := factor(station)]

par(mfcol = c(3, 2),
    mar = c(4, 4, 1, 1) + 0.1)
diff_dat[order(site, station), 
         {
           smoothScatter(sea_level_pressure_m, 
                         ex_sea_level_pressure_m,
                         main = paste(.BY[[1]], .BY[[2]]))
           abline(0, 1)
         },
         by = .(site, station)]
par(mfrow = c(1, 1))

diff_dat[, .(corr_pearson = cor(sea_level_pressure_m, ex_sea_level_pressure_m),
             corr_spearman = cor(sea_level_pressure_m, ex_sea_level_pressure_m, method = "spearman")),
         by = .(site, station)]

# Air Density Correction --------------------------------------------------
# Solinst confirmed they use a constant water density of 1 kg/m^3, so the 
# temperature compensation addresses the transducer and electronics. This is 
# probably fine in a submerged setting in a well where there is little diurnal
# temperature flucation. But using the levellogger as a barologger means that
# it is necessary to correct for the diurnal temperature changes and the effect
# on air density, and on the 'water density' that the levellogger is actually
# measuring. I performed a multistep correction (commented code below) to convert recorded pressure to
# absolute pressure to head of air (m) to head of water (m).

# P = h_w * g * 1000
# rho_a = P / (T_a * R)
# h_a = P / (p_a * g)
# h_w,r = h_a * (rho_a / rho_w)

# P = absolute pressure in Pascals
# rho_w = density of water (kg/m3) as a function of temperature (C) [@tanaka-2001]
# rho_a = density of air as a function of temperature and absolute pressure 
# h_a = head of air (meters)
# h_w = head of water (meters)
# h_w,r = rectified head of water (meters)
# T_a = tempearture of air (K)
# R = ideal gas constant for air (287.058 J/kgT, T in Kelvins)
# g = gravitational acceleration (9.80665 m/s)

# The above equations can be combined and reduced into a single simplified 
# equation:

# h_w,r = h_w * (1000 / rho_w)

# Water density calculation is taken from @tanaka-2001

water_density <- 
  function(t){
    # From [@tanaka-2001]
    
    a1 <- -3.983035
    a2 <- 301.797
    a3 <- 522528.9
    a4 <- 69.34881
    a5 <- 999.97495
    
    a5 * (1 - ((t + a1)^2*(t+a2)) / (a3*(t+a4)))
  }

# curve(water_density(x), -40, 40)
# abline(h = water_density(-7))

# Calculate density of water
diff_dat[, rho_water_kg_m3 := water_density(air_temperature_c)]

# Something along this equation is necessary to avoid a break in the linear 
# temperature compensation between Solinst and external logger pressures
# -7 was just chosen visually
diff_dat[air_temperature_c < -10, rho_water_kg_m3 := water_density(-10)]

# Rectify barometric pressures (remove elevation adjustment to pressure becuase
# this was not measured by the transducer but was added to adjust recorded
# measurements to sea level pressure)
diff_dat[, elevation_offset := ifelse(site == "053", 473 / 826, 503 / 826)]

diff_dat[, rect_sea_level_pressure_m := (sea_level_pressure_m - elevation_offset) * (1000 / rho_water_kg_m3)]

# Add elevation adjustment back in 
diff_dat[, rect_sea_level_pressure_m := rect_sea_level_pressure_m + elevation_offset]

# # Full Steps for Rectification
# # Set constants
# g_m_s <- 9.80665
# R_J_Kg_K <- 287.058

# # Calculate absolute pressure (air pressure) in Pascals
# # Solinst uses 1000 kg/L constant density
# diff_dat[, abs_pressure_pa := (sea_level_pressure_m - elevation_offset) * g_m_s * 1000]
# 
# # Calculate Air Density
# diff_dat[, rho_air_kg_m3 := abs_pressure_pa / ((air_temperature_c + 273.15) * R_J_Kg_K)]
# 
# # Calculate pressure in meters of air
# diff_dat[, adj_pressure_m_air := abs_pressure_pa / (rho_air_kg_m3 * g_m_s) ]
# 
# # Convert Meters of Air to Meters of Water
# diff_dat[, rect_sea_level_pressure_m := adj_pressure_m_air / (rho_water_kg_m3 / rho_air_kg_m3)]
# 
# # Add elevation adjustment back in 
# diff_dat[, rect_sea_level_pressure_m := rect_sea_level_pressure_m + elevation_offset]

# Check Correlation betwen External and Rectified Pressure

par(mfcol = c(3, 2))
diff_dat[order(site, station), 
         {
           smoothScatter(rect_sea_level_pressure_m, 
                         ex_sea_level_pressure_m,
                         main = paste(.BY[[1]], .BY[[2]]))
           abline(0, 1)
         },
         by = .(site, station)]
par(mfrow = c(1, 1))

diff_dat[order(site, station), 
         .(correlation = cor(sea_level_pressure_m, ex_sea_level_pressure_m),
           dens_adj_correlation = cor(rect_sea_level_pressure_m, ex_sea_level_pressure_m)),
         by = .(site, station)]

# Temperature Compensation ------------------------------------------------
# If we assume that the external pressure is an authoritative source then we can
# calculate the error of Solinst measurements compared to external measurements
# to check for temperature artifacts [@moore-2016]

diff_dat[, pressure_error_m := ex_sea_level_pressure_m - rect_sea_level_pressure_m]

par(mfcol = c(3, 2))
diff_dat[order(site, station), 
         {
           smoothScatter(air_temperature_c, 
                         pressure_error_m,
                         main = paste(.BY[[1]], .BY[[2]]))
         },
         by = .(site, station)]
par(mfrow = c(1, 1))

# There is now a linear correction to be applied to the loggers (just a regular
# temperature compenstation). We can find this by comparing the pressure at the
# airport stations with the pressure measured at the reasearch sites. This is
# visible as the linear trend in the relationship between temperature difference
# and air temperature, with the x-intercept showing where the solinst loggers 
# record the 'correct' pressure. After that it is a matter of applying the 
# approaches in [@moore-2016]
# There is an issue below ~ -10 C, likely becuase the density of water should
# just be held constant there as the @tanaka-2001 equation only includes liquid
# water. Using the generic density of water ice results in a sharp drop in 
# density and severe discontinuties in the rectified data. Rectifying the data
# to express pressure in meters of water requires some extrapolation as liquid
# water would not exist below 0 C, and the Solinst temperature compensation only
# goes from 0 to 40 C.

# Fit linear model to each site-station combination to see if all sites show 
# similar relationships

temp_mods <- 
  diff_dat[, .(mod = list(lm(pressure_error_m ~ 
                               air_temperature_c, 
                             data = .SD))),
           by = .(site, station)]

temp_mods[, `:=`(intercept = vapply(mod, 
                                    function(x){coef(x)[1]},
                                    FUN.VALUE = numeric(1)),
                 slope = vapply(mod, 
                                function(x){coef(x)[2]}, 
                                FUN.VALUE = numeric(1)))]

temp_mods[, x_intercept := -intercept / slope]

temp_mods[order(site, station)]

# Remove KIWD from Datasets -----------------------------------------------

# It looks like KIWD has a different relationship to these sites than KLNL and
# E6920. This could be because KIWD shows more of an impact from Lake Superior
# as it is the closest station to the lake.

# library(sf)
# library(mapview)
# stations <- read_sf("data/__supplementary_baro_stations.gpkg")
# mapview(stations)

diff_dat <- 
  diff_dat[station != "KIWD"]

dat <- 
  dat[site != "KIWD"]

# Develop Genearlized Model Per Logger ------------------------------------

# Fitting a linear mixed model for each site (each logger since they are always
# deployed to the same place). The data are grouped by station as a random 
# effect to allow for the fixed effect model to incorporate the uncertainty from
# each site-station pair while allow prediction without the station term.

temp_mmods <- 
  diff_dat[, .(mod = list(lmer(pressure_error_m ~ 
                                 air_temperature_c + 
                                 (0 + air_temperature_c | station), 
                               data = .SD))),
           by = .(site)]

# Singular fit if intercept and slope can vary. Not so if only slope can vary.

# Extract the fixed effects terms for intercept and slope
temp_mmods[, `:=`(intercept = vapply(mod, 
                                     function(x){fixef(x)[1]},
                                     FUN.VALUE = numeric(1)),
                  slope = vapply(mod, 
                                 function(x){fixef(x)[2]}, 
                                 FUN.VALUE = numeric(1)))]

# Calculate the x-intercept, which represents the point where the Solinst and
# external pressures agree (error = 0)
temp_mmods[, x_intercept := -intercept / slope]


# Apply temperature compensation using equation 3 from @moore-2001
diff_dat[temp_mmods, 
         temp_rect_sea_level_pressure_m := rect_sea_level_pressure_m + slope*(air_temperature_c - x_intercept),
         on = c("site")]

diff_dat[, pressure_error_m := ex_sea_level_pressure_m - temp_rect_sea_level_pressure_m]

# Apply physics based correction
# par(mfcol = c(3, 2))
# diff_dat[, {smoothScatter(air_temperature_c, 
#                           ex_sea_level_pressure_m - (rect_sea_level_pressure_m + (9.80665 * 1000) / ((air_temperature_c + 273.15) * 287.058)), 
#                           ylab = "error (m)",
#                           main = paste(.BY[[1]], .BY[[2]])); 
#   abline(h = 0)},
#          by = .(site, station)]
# par(mfrow = c(1, 1))
# 
# diff_dat[, temp_rect_sea_level_pressure_m_physics := rect_sea_level_pressure_m + (1000) / ((air_temperature_c + 273.15) * 287.058)]
# diff_dat[, physics_error_m := ex_sea_level_pressure_m - temp_rect_sea_level_pressure_m_physics]
# diff_dat[, .(error = mean(physics_error_m), offset = mean(elevation_offset)), by = .(site, station)]

# par(mfrow = c(2, 2))
# diff_dat[, {smoothScatter(temp_rect_sea_level_pressure_m_physics, ex_sea_level_pressure_m); abline(0, 1)},
#          by = .(site, station)]
# par(mfrow = c(1, 1))
# 
# par(mfrow = c(2, 2))
# diff_dat[, {smoothScatter(air_temperature_c, physics_error_m); abline(h = 0)},
#          by = .(site, station)]
# par(mfrow = c(1, 1))

# Compare Solinst and external pressure correlation
diff_dat[, .(raw_correlation = cor(sea_level_pressure_m, ex_sea_level_pressure_m),
             rectified_correlation = cor(rect_sea_level_pressure_m, ex_sea_level_pressure_m),
             temp_rectified_correlation = cor(temp_rect_sea_level_pressure_m, ex_sea_level_pressure_m)),
         by = .(site, station)]

par(mfcol = c(2, 2))
diff_dat[order(site, station), 
         {
           smoothScatter(temp_rect_sea_level_pressure_m, 
                         ex_sea_level_pressure_m,
                         main = paste(.BY[[1]], .BY[[2]]))
           abline(0, 1)
         },
         by = .(site, station)]

diff_dat[order(site, station), 
         {
           smoothScatter(air_temperature_c, 
                         pressure_error_m,
                         main = paste(.BY[[1]], .BY[[2]]))
           abline(h = 0)
         },
         by = .(site, station)]

diff_dat[order(site, station), 
         {
           plot(density(pressure_error_m),
                main = paste(.BY[[1]], .BY[[2]]))
           abline(h = 0)
         },
         by = .(site, station)]

par(mfrow = c(1, 1))

# Evaluate Corrections ----------------------------------------------------

# Calculate total correction applied by rectifying for density and temperature
diff_dat[, correction_m := temp_rect_sea_level_pressure_m - sea_level_pressure_m]

# The correction matches the original difference in external and Solinst
# pressures
par(mfcol = c(3, 2))
diff_dat[, {
  smoothScatter(pressure_diff, 
                correction_m,
                main = paste(.BY[[1]], .BY[[2]]), 
                colramp = colorRampPalette(c("#FFFFFF", "#262626")))
  abline(0, 1, 
         col = "#82c729")},
  by = .(site, station)]
par(mfrow = c(1, 1))

# Check Models on Full Data -----------------------------------------------

dat[, elevation_offset_m := ifelse(site == "053", 473 / 826, 503 / 826)]
dat[!(site %in% c("053", "113")), elevation_offset_m := 0]

# Rectifiy original data by Density
dat[site %in% c("053", "113"),
    sea_level_pressure_m := elevation_offset_m + ((sea_level_pressure_m - elevation_offset_m) * (1000 / water_density(ifelse(air_temperature_c < -7, -7, air_temperature_c))))]

# Rectify original data by temperature
dat[temp_mmods,
    sea_level_pressure_m := sea_level_pressure_m + slope*(air_temperature_c - x_intercept),
    on = c("site")]

# Check retification
cor(dcast(dat, 
          sample_time ~ site, 
          value.var = "sea_level_pressure_m")[, -c(1)], 
    use = "pairwise.complete.obs")[1:2, 3:4]

# The correction matches the original difference in external and Solinst
# pressures
par(mfrow = c(2, 2))
mapply(function(x, y){
  df <- 
    dcast(dat,
          sample_time ~ site,
          value.var = "sea_level_pressure_m")
  
  assign(x, df[[x]])
  assign(y, df[[y]])
  
  smoothScatter(get(x), get(y), 
                main = paste(x, y),
                xlab = "Solinst Pressure (m)", ylab = "External Pressure (m)")
  abline(0, 1)},
  list("053", "113"),
  rep(grep("[A-Z]+", unique(dat$site), value = TRUE), 2))
par(mfrow = c(1, 1))

dat[, elevation_offset_m := NULL]


# Save Models -------------------------------------------------------------

# Add serial number to temperature models
temp_mmods[serial_numbers, baro_serial_number := serial_number, on = "site"]

saveRDS(temp_mmods, 
        "output/models/solinst_barologger_temperature_correction_models.rds")
