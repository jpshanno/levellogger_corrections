library(data.table)
library(httr)

url <- 
  paste("https://api.synopticdata.com/v2/stations/timeseries?stid=F3311,KCMX",
        paste0("token=", Sys.getenv("SYNOPTIC_KEY")),
        paste0("start=", 202005290000),
        paste0("end=", 202006012359),
        "obtimezone=UTC",
        "units=temp|C,pres|pa",
        "vars=air_temp,pressure,sea_level_pressure",
        sep = "&")

response <- 
  POST(url)

query_content <- 
  content(response, 
          type = "application/json")

dat <- 
  lapply(seq_along(query_content[["STATION"]]), 
         function(x){
           dat <- query_content[["STATION"]][[x]]
           
           cols <- 
             c("date_time", "pressure_set_1d", 
               "sea_level_pressure_set_1d", "air_temp_set_1")
           
           obs <- 
             as.data.table(lapply(dat$OBSERVATIONS[cols], unlist))
           
           info <- 
             data.table(station = dat$STID,
                        elevation_ft = as.numeric(dat$ELEV_DEM))
           
           cbind(info, obs)})

dat <- 
  rbindlist(dat)

dat <- 
  dat[, .(station,
          sample_time = as.POSIXct(date_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
          air_temperature_c = air_temp_set_1,
          air_pressure_m = (sea_level_pressure_set_1d / 9806.65) - 0.00121 * (elevation_ft / 3.2808),
          sea_level_pressure_m = sea_level_pressure_set_1d / 9806.65)]

attributes(dat$sample_time)$tzone <- "EST5EDT"

setkey(dat, station, sample_time)

write.table(dat, 
            "external_barometric_data.csv", 
            row.names = FALSE, 
            append = TRUE,
            sep = ",",
            dec = ".",
            qmethod = "double",
            col.names = FALSE)
