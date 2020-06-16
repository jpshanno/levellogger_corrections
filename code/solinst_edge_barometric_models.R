library(data.table)
library(ingestr)
library(lme4)

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
  list.files(pattern = "baro.xle")

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
           
           dat[, .(baro_serial_number = serial_number,
                   sample_time = as.POSIXct(paste(sample_date, sample_time),
                                            tz = "EST5EDT"),
                   air_temperature_c = temperature_c,
                   pressure_m = level_m,
                   sea_level_pressure_m)]
           })

baro <- 
  rbindlist(baro)

# Compensate for density changes
baro[, `:=`(pressure_m = pressure_m * 1000 / water_density(air_temperature_c),
            sea_level_pressure_m = sea_level_pressure_m * 1000 / water_density(air_temperature_c))]

all_baro <- 
  ex_baro[baro,
          mult = "all",
          nomatch = NULL,
          on = "sample_time"]

par(mfrow = c(3, 2),
    mar = c(4, 2, 2, 2) + 0.1)
all_baro[, plot(sea_level_pressure_m ~ ex_sea_level_pressure_m, 
                main = paste(.BY[[1]], .BY[[2]])),
         by = .(station, baro_serial_number)]
par(mfrow = c(1, 1))
