library(data.table)
library(ingestr)
library(lme4)
library(lubridate)

water_density <- 
  function(t){
    # From [@tanaka-2001]
    t <- ifelse(t < -7, -7, t)
    
    a1 <- -3.983035
    a2 <- 301.797
    a3 <- 522528.9
    a4 <- 69.34881
    a5 <- 999.97495
    
    a5 * (1 - ((t + a1)^2*(t+a2)) / (a3*(t+a4)))
  }

ex_baro <- 
  fread("external_barometric_data.csv")

ex_baro[,sample_time := as.POSIXct(sample_time, tz = "EST5EDT")]

setnames(ex_baro, 
         c("air_temperature_c", "air_pressure_m", "sea_level_pressure_m"),
         paste0("ex_", c("air_temperature_c", "air_pressure_m", "sea_level_pressure_m")))

baro_files <- 
  list.files(pattern = "20200601_baro")

baro <- 
  lapply(baro_files,
         function(x){
           dat <- 
             merge(as.data.table(ingest_xle(x)), 
                   as.data.table(ingest_header(x)),
                   by = "input_source")
           
           if("level_cm" %in% names(dat)){
             dat[, level_m := level_cm / 100]
           }
           
           if("altitude_m" %in% names(dat)){
             dat[, level_m := level_m - altitude_m * 0.00121]
           }
           
           if("temperature_deg_c" %in% names(dat)){
             dat[, temperature_c := temperature_deg_c]
           }
           
           if(dat$instrument_type[1] == "LT_Jr"){
             dat[, level_m := level_m + 9.5]
           }
           
           dat[, sea_level_pressure_m := level_m + 231 * 0.00121]
           
           dat[,sample_time := as.POSIXct(paste(sample_date, sample_time),
                                         tz = "EST5EDT")]
           
           dat[sample_time <= as.POSIXct("2020-06-01 08:05:00",
                                        tz = "EST5EDT"), 
               .(baro_serial_number = serial_number,
                   sample_time,
                   air_temperature_c = temperature_c,
                   pressure_m = level_m,
                   sea_level_pressure_m)]
           })

baro <- 
  rbindlist(baro)

# Compensate for density changes
baro[, `:=`(raw_pressure_m = pressure_m,
            pressure_m = pressure_m * 1000 / water_density(air_temperature_c),
            raw_sea_level_pressure_m = sea_level_pressure_m,
            sea_level_pressure_m = sea_level_pressure_m * 1000 / water_density(air_temperature_c))]

all_baro <- 
  ex_baro[baro,
          mult = "all",
          nomatch = NULL,
          on = "sample_time"]

all_baro[, `:=`(raw_pressure_error_m = ex_sea_level_pressure_m - raw_sea_level_pressure_m,
                pressure_error_m = ex_sea_level_pressure_m - sea_level_pressure_m)]

par(mfrow = c(4, 2),
    mar = c(4, 2, 2, 2) + 0.1)
all_baro[, smoothScatter(air_temperature_c, comp_pressure_diff),
         by = .(station, baro_serial_number)]
par(mfrow = c(1, 1))

temp_mmods <- 
  all_baro[, .(mod = list(lmer(pressure_error_m ~ 
                                 air_temperature_c + 
                                 (air_temperature_c | station), 
                               data = .SD))),
           by = .(baro_serial_number)]
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

par(mfrow = c(2, 2),
    mar = c(4, 2, 2, 2) + 0.1)
temp_mmods[, smoothScatter(mod[[1]]@frame$air_temperature_c, resid(mod)),
         by = .(baro_serial_number)]
par(mfrow = c(1, 1))

# Apply temperature compensation using equation 3 from @moore-2001
all_baro[temp_mmods, 
         rect_sea_level_pressure_m := sea_level_pressure_m - slope*(air_temperature_c - x_intercept),
         on = c("baro_serial_number")]

# Compare Solinst and external pressure correlation
all_baro[, .(raw_correlation = cor(raw_sea_level_pressure_m, ex_sea_level_pressure_m),
             rectified_correlation = cor(sea_level_pressure_m, ex_sea_level_pressure_m),
             temp_rectified_correlation = cor(rect_sea_level_pressure_m, ex_sea_level_pressure_m)),
         by = .(baro_serial_number, station)]

# Check Residuals

all_baro[temp_mmods, 
         pred_pressure_error_m := intercept + slope * air_temperature_c,
         on = c("baro_serial_number")]

all_baro[, pressure_error_resid_m := pressure_error_m - pred_pressure_error_m]

par(mfrow = c(2, 2))
all_baro[, {plot(air_temperature_c, pressure_error_resid_m)},
         by = .(baro_serial_number)]
par(mfrow = c(1, 1))

saveRDS(temp_mmods[baro_serial_number %in% c("2025928", "2013939")], "solinst_edge_barometric_models.rds")
