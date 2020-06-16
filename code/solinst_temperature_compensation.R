# Check manual measurements. Looks like two bad ones for stage

# There might be an issue with 2018 stage at 053 through ~early August. After I
# checked on and downloaded the loggers on 2018/8/8 the pattern of the signal
# changed and the large diurnal flucuations dropped off. Maybe there was an air
# bubble or something in the logger?

# library(ggplot2)
library(data.table)
library(ingestr)
library(tydygraphs)

dy_graph <- 
  function(...){
  dygraph(...) %>% 
  dyOptions(useDataTimezone = TRUE) %>% 
  dyRangeSelector()}

# Load & Prep Data --------------------------------------------------------
# All of these files have no altitude offset

flume_file <- 
  "flume_053_2018-09-05.xle"

well_file <- 
  "well_053_2019-04-19.xle"

baro_file <- 
  "baro_053_2018-09-05.xle"

# Flume levellogger
fl <-
  as.data.table(ingest_xle(flume_file))

fl <- 
  fl[, .(sample_time = as.POSIXct(paste(sample_date, sample_time),
                                       format = "%Y/%m/%d %H:%M:%S",
                                       tz = "EST"),
         stage_m = level_m,
         flume_temperature_c = temperature_c)]

setkey(fl, sample_time)

# Wetland level
wl <-
  as.data.table(ingest_xle(well_file))

wl <- 
  wl[, .(sample_time = as.POSIXct(paste(sample_date, sample_time),
                                       format = "%Y/%m/%d %H:%M:%S",
                                       tz = "EST"),
         level_m = level_m,
         well_temperature_c = temperature_c)]

setkey(wl, sample_time)

# Load Baro

baro <-
  as.data.table(ingest_xle(baro_file))

baro <-
  baro[,.(sample_time = as.POSIXct(paste(sample_date, sample_time),
                                        format = "%Y/%m/%d %H:%M:%S",
                                        tz = "EST"),
          baro_m = level_m + 9.5,
          air_temperature_c = temperature_deg_c,
          baro_input_source = input_source)]

setkey(baro, sample_time)

# Combine barometric, well level, and flume stage data
combined <-
  fl[wl][baro, nomatch = NULL]

# Raw Compensation

combined[, compensated_stage_m := stage_m - baro_m]

dy_graph(combined, compensated_stage_m)

# Density-Rectified Compensation

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

combined[, density_compensated_baro_m :=  baro_m * (1000 / water_density(ifelse(air_temperature_c < -7, -7, air_temperature_c)))]
combined[, density_compensated_stage_m := stage_m - density_compensated_baro_m]

dy_graph(combined, compensated_stage_m, density_compensated_stage_m)


# Add in Elevation



# Temperature-Corrected Compensation



# Full Compensation

mods <- 
  readRDS("../../output/solinst_barologger_tempeature_correction_models.rds")

combined[, site := "053"]

combined[mods, 
         `:=`(full_compensated_baro_m = (density_compensated_baro_m) + slope * (air_temperature_c - x_intercept)),
         on = "site"]

combined[, full_compensated_stage_m := stage_m - full_compensated_baro_m]

dy_graph(combined, compensated_stage_m, 
         density_compensated_stage_m,
         full_compensated_stage_m)

combined[, `:=`(compensated_level_m = level_m - baro_m,
                density_compensated_level_m = level_m - density_compensated_baro_m,
                full_compensated_level_m = level_m - full_compensated_baro_m)]

dy_graph(combined, 
         compensated_level_m, 
         density_compensated_level_m,
         full_compensated_level_m)



# Compensate Temperature Differences --------------------------------------

# Calculate the temperature difference adjustment factor (derived from 
# Gay-Lussac's Law, P1T1 = P2T2)
combined[, `:=`(level_temp_adjustment_factor = (air_temperature_c + 273.15) / (well_temperature_c + 273.15),
                stage_temp_adjustment_factor = (air_temperature_c + 273.15) / (flume_temperature_c + 273.15))]

# Calculate adjusted compensated level
combined[, `:=`(diff_adj_comp_level_cm = full_compensated_level_m * level_temp_adjustment_factor,
               diff_adj_comp_stage_cm = full_compensated_stage_m * stage_temp_adjustment_factor)]

# Compare base compensation with temperature adjusted compensation
# This does help reduce the strength of the diurnal signal
dy_graph(combined, 
        compensated_level_m, density_compensated_level_m, full_compensated_level_m, diff_adj_comp_level_cm, 
        main = "Well Level")

solinst %>% 
  jpshanno::spike_test(comp_level_cm) %>% 
  jpshanno::spike_test(adj_comp_level_cm) %>% 
  dygraph(comp_level_cm, adj_comp_level_cm, main = "Well Level") %>% 
  dyRangeSelector()

precip <- fread("../../../QAQC/Continuous_Sampling/precipitation_daily.csv")[, sample_date:=as.IDate(sample_date)]
precip[site == 53 & year(sample_date) == 2018, .(sample_date, precip_cm)][solinst, on = "sample_date"] %>% spike_test(comp_level_cm) %>% spike_test(adj_comp_level_cm) %>% dygraph(comp_level_cm, adj_comp_level_cm, precip_cm, time = "sample_time") %>% dyOptions(useDataTimezone = TRUE) %>% dySeries("precip_cm", axis = "y2") %>% dyRangeSelector()

dygraph(solinst, comp_stage_cm, diff_adj_comp_stage_cm, # flume_temperature_c, baro_temperature_c,
        main = "Flume Stage",
        ylab = "Compensated (uncorrected) Stage (cm)") %>% 
  # dySeries("flume_temperature_c",
  #          axis = "y2") %>% 
  # dySeries("baro_temperature_c",
  #          axis = "y2") %>% 
  dyRangeSelector()

# The diurnal signal is still reversed from anticipated. 
#I looked at compensating water and air to a reference temperature (tried 0 and
#25 C). It does shift the water levels (down for 0 C, up for 25 C), but doesn't
#seem to change much else
# solinst[, `:=`(water_adjustment_factor =  (0 + 273.15) / (temperature_c + 273.15),
#                baro_adjustment_factor =  (0 + 273.15) / (baro_temperature_c + 273.15))]
# 
# solinst[, dual_adj_comp_level_cm := 100 * (level_m * water_adjustment_factor - baro_level_m * baro_adjustment_factor)]
#  
# dygraph(solinst, comp_level_cm, adj_comp_level_cm, dual_adj_comp_level_cm) %>%
#   dyRangeSelector()

# Also tried to compensate to daily maximum water temp, which doesn't have a 
# big impact.
# solinst[, doy := yday(sample_time)]
# 
# solinst[solinst[, .(max_flume_temp = max(flume_temperature_c)), by = .(doy)],
#         on = "doy",
#         `:=`(max_flume_temp = max_flume_temp,
#              water_adjustment_factor =  (max_flume_temp + 273.15) / (flume_temperature_c + 273.15),
#              baro_adjustment_factor =  (max_flume_temp + 273.15) / (baro_temperature_c + 273.15))]
# 
# solinst[, dual_adj_comp_stage_cm := 100 * (stage_m * water_adjustment_factor - baro_m * baro_adjustment_factor)]
# 
# dygraph(solinst, comp_stage_cm, adj_comp_stage_cm, dual_adj_comp_stage_cm) %>% 
#   dyRangeSelector()



# Look at diurnal range for base and adjusted data
diurnal_ranges <- 
  solinst[, .(stage = diff(range(comp_stage_cm)),
              temp_adj_stage = diff(range(adj_comp_stage_cm)),
              level = diff(range(comp_level_cm)),
              temp_adj_level = diff(range(adj_comp_level_cm))), 
          by = .(doy = yday(sample_time))]

diurnal_ranges <- 
  melt(diurnal_ranges, 
       id.vars = "doy",
       variable.name = "measurement",
       value.name = "diurnal_range")

ggplot(data = diurnal_ranges, 
       aes(x = measurement, y = diurnal_range)) +
  geom_boxplot() +
  scale_y_log10()

# Load Manual Readings ----------------------------------------------------

manual_readings <- 
  fread("well_and_flume_measurements_053.csv")

manual_readings <- 
  manual_readings[comments == "",
                  .(sample_time = as.POSIXct(paste(sample_date, sample_time), 
                                             format = "%m/%d/%Y %H:%M",
                                             tz = "EST"),
                    location,
                    manual_level_cm)]

manual_readings <- 
  dcast(manual_readings, 
        sample_time ~ location, 
        value.var = "manual_level_cm")

setnames(manual_readings, 
         c("sample_time", "manual_stage_cm", "manual_level_cm"))

manual_readings[, manual_level_cm := 104.9 - manual_level_cm]

setkey(manual_readings, sample_time)

# Calculate & Compare Offsets ---------------------------------------------

manual <- 
  solinst[manual_readings, nomatch = NULL]

manual[manual_level_cm > 0,
       `:=`(diff_level_cm = comp_level_cm - manual_level_cm,
            adj_diff_level_cm = adj_comp_level_cm - manual_level_cm,
            diff_adj_diff_level_cm = diff_adj_comp_level_cm - manual_level_cm,
            diff_stage_cm = comp_stage_cm - manual_stage_cm,
            adj_diff_stage_cm = adj_comp_stage_cm - manual_stage_cm,
            diff_adj_diff_stage_cm = diff_adj_comp_stage_cm - manual_stage_cm)]

# The adjusted offsets show a lower coefficient of variation for water level, but
# higher for stage. Perhaps because the flume water levels vary with a similar
# signal as the air temperature

manual[, .(diff_stage_cm = sd(diff_stage_cm, na.rm = TRUE) / mean(diff_stage_cm, na.rm = TRUE),
           adj_diff_stage_cm = sd(adj_diff_stage_cm, na.rm = TRUE) / mean(adj_diff_stage_cm, na.rm = TRUE),
           diff_adj_diff_stage_cm = sd(diff_adj_diff_stage_cm, na.rm = TRUE) / mean(diff_adj_diff_stage_cm, na.rm = TRUE),
           diff_level_cm = sd(diff_level_cm, na.rm = TRUE) / mean(diff_level_cm, na.rm = TRUE),
           adj_diff_level_cm = sd(adj_diff_level_cm, na.rm = TRUE) / mean(adj_diff_level_cm, na.rm = TRUE),
           diff_adj_diff_level_cm = sd(diff_adj_diff_level_cm, na.rm = TRUE) / mean(diff_adj_diff_level_cm, na.rm = TRUE))]

offsets <- 
  manual[, 
         lapply(.SD, mean, na.rm = TRUE), 
         .SDcols = c("diff_stage_cm", "adj_diff_stage_cm", "diff_adj_diff_stage_cm", 
                     "diff_level_cm", "adj_diff_level_cm", "diff_adj_diff_level_cm")]

setnames(offsets, c("offset_stage", "offset_adj_stage", "offset_diff_adj_stage",
                    "offset_level", "offset_adj_level", "offset_diff_adj_level"))

# Correct Water Levels ----------------------------------------------------

corrected <- 
  cbind(solinst, offsets)

corrected[, `:=`(corrected_stage_cm = comp_stage_cm - offset_stage,
                 adj_corrected_stage_cm = adj_comp_stage_cm - offset_adj_stage,
                 diff_adj_corrected_stage_cm = diff_adj_comp_stage_cm - offset_diff_adj_stage,
                 corrected_level_cm = comp_level_cm - offset_level,
                 adj_corrected_level_cm = adj_comp_level_cm - offset_adj_level,
                 diff_adj_corrected_level_cm = diff_adj_comp_level_cm - offset_diff_adj_level)]

# corrected[corrected_stage_cm <0, corrected_stage_cm := 0]
# corrected[adj_corrected_stage_cm <0, adj_corrected_stage_cm := 0]

corrected[manual_readings, 
          `:=`(manual_level_cm = manual_level_cm,
               manual_stage_cm = manual_stage_cm)]

dygraph(corrected,
        corrected_level_cm, 
        diff_adj_corrected_level_cm,
        manual_level_cm,
        main = "Well Level") %>% 
  dySeries("manual_level_cm", 
           pointSize = 5) %>% 
  dyOptions(useDataTimezone = TRUE) %>% 
  dyRangeSelector()

dygraph(corrected,
        corrected_stage_cm, 
        diff_adj_corrected_stage_cm,
        manual_stage_cm,
        main = "Flume Stage") %>% 
  dySeries("manual_stage_cm", 
           pointSize = 5) %>% 
  dyOptions(useDataTimezone = TRUE) %>% 
  dyRangeSelector()
     
# Temp and baro difference are inversely proportional
solinst[, plot(adj_comp_level_cm ~ temperature_c)]
, 
   .(sample_timestamp, 
     temp_diff = as.numeric(scale(temperature_c - baro_temperature_c)), 
     level_diff = as.numeric(scale(level_m - baro_level_m)))] %>%
  dygraph(temp_diff, level_diff)

wl[baro, 
   plot(as.numeric(scale(temperature_c - baro_temperature_c)), 
        as.numeric(scale(level_m - baro_level_m)))]

# Adjusting baro pressure accourding to a combination of Boyle's and Charles'
# Laws so that P1V1/T1 = P2V2/T2. V1 and V2 are set equal as the atmosphere is
# not confined to a given volume (see USGS document, p.2). Now let's substitute
# in the subscripts with 'b' for barometer reading and 'w' for water reading,
# and 'a' for barometer reading at water temperature: Pb/Tb = Pa/Tw. Solving for
# Pa gives us Pa = PbTw/Tb

wl[baro, .(sample_timestamp, comp_level = 100*((level_m - baro_level_m)-9.5))] %>% 
  dygraph(comp_level) %>% 
  dyRangeSelector()


wl[baro,
   .(sample_timestamp,
     adj_baro = baro_level_m * (temperature_c - 273.15) / (baro_temperature_c - 273.15),
     baro_level_m,
     level_m = 100 * (level_m - 9.5),
     base_comp_level = 100 * (level_m - baro_level_m - 9.5),
     comp_level = 100*(level_m - baro_level_m * (temperature_c - 273.15) / (baro_temperature_c - 273.15)-9.5))] %>%
  dygraph(comp_level, base_comp_level) %>% 
  dyRangeSelector()

# Trying to use correction factors for gas meters
test <- 
  wl[baro]

test[, `:=`(temp_correction = ((25-273.15) / (temperature_c - 273.15)),
            baro_temp_correction = ((25-273.15) / (baro_temperature_c - 273.15)))]

test[, `:=`(temp_adj_level_m = level_m * temp_correction,
            temp_adj_baro_level_m = baro_level_m * baro_temp_correction)]
     
test[, `:=`(comp_level = level_m - baro_level_m - 9.5,
            temp_adj_comp_level = temp_adj_baro_level_m - temp_adj_level_m + 9.5)]

# test[, .(sample_timestamp, 
#          temp_adj_level_m = temp_adj_level_m - 9.5, 
#          temp_adj_baro_level_m)][, 
#   .(sample_timestamp, temp_adj_level_m, temp_adj_baro_level_m, 
#     final_level = temp_adj_baro_level_m- temp_adj_level_m)]


wl[baro][, .(sample_timestamp, comp_level_cm = 100 * (level_m - baro_level_m - 9.5), adj_comp_level_cm = 100 * (level_m - baro_level_m * (temperature_c - 273.15) / (baro_temperature_c - 273.15) - 9.5))][manual,.(sample_timestamp, manual_level_cm = 100 * manual_level_m, comp_level_cm, diff = comp_level_cm - 100*manual_level_m, adj_diff = adj_comp_level_cm - 100 * manual_level_m)][, .(mean(diff, na.rm = TRUE), mean(adj_diff, na.rm = TRUE))]
wl[baro][, .(sample_timestamp, comp_level_cm = 100 * (level_m - baro_level_m - 9.5), adj_comp_level_cm = 100 * (level_m - baro_level_m * (temperature_c - 273.15) / (baro_temperature_c - 273.15) - 9.5))][manual,.(sample_timestamp, manual_level_cm = 100 * manual_level_m, comp_level_cm, diff = comp_level_cm - 100*manual_level_m, adj_diff = adj_comp_level_cm - 100 * manual_level_m)][, .(mean(diff, na.rm = TRUE), mean(adj_diff, na.rm = TRUE))]


manual <- 
  fread("well_and_flume_measurements.csv")[,
    .(sample_timestamp = as.POSIXct(paste(sample_date, sample_time), 
                                   format = "%m/%d/%Y %H:%M",
                                   tz = "EST"),
      manual_level_m = manual_level / 3.2808)]

setkey(manual, sample_timestamp)

test[manual,
     .(sample_timestamp,
       manual_level_m,
       temp_adj_comp_level,
       diff = temp_adj_comp_level - manual_level_m)][, lm(manual_level_m ~ temp_adj_comp_level)]

manual[test,
       .(sample_timestamp,
         manual_level_m,
         temp_adj_comp_level,
         corrected_level = -0.04144 + 0.15153 * temp_adj_comp_level)][,
        .(sample_timestamp, 
          manual_level_m,
          corrected_level = ifelse(corrected_level < 0, 0, corrected_level))] %>% 
  dygraph(corrected_level,
          manual_level_m)


solinst[manual, diff = comp_level_cm - 100*manual_level_m, adj_diff = adj_comp_level_cm - 100 * manual_level_m)][, .(mean(diff, na.rm = TRUE), mean(adj_diff, na.rm = TRUE))]

test[, 
     .(sample_timestamp,
       temp_adj_level_m,
       temp_adj_baro_level_m,
       final_level = temp_adj_comp_level + 9.5)] %>% 
  dygraph(final_level) %>% 
  dyRangeSelector()

# The level below is adjusted using the density of air equation from wikipedia
# and solving for pressure -> it may not be the right conversion, need to 
# confirm it with some manual measurements. Also should look to see if this is
# the standard formula for the relationship between density and pressure of any
# liquid



dat <- 
  lapply(dat.list,
         function(x){x$density <- gsub("(test|\\.xle)", "", x$input_source); x})

dat <- 

raw <- 
  as.data.table(ingest_xle("Desktop/test.xle", header.info = FALSE, collapse.timestamp = TRUE))

raw <- 
  as.data.table(ingest_xle("Desktop/test.xle", header.info = FALSE, collapse.timestamp = TRUE))

full_join(raw, adj, by = "sample_timestamp", suffix = c(".raw", ".adj")) %>% mutate(level_ratio = level_cm.adj / level_cm.raw, level_diff = level_cm.adj - level_cm.raw) %>% lm(level_cm.adj ~ level_cm.raw, data = .)

dat <- bind_rows(raw, adj)

# Exploring Density Adjustment from Solinst -------------------------------

files <-
  setNames(nm = list.files(pattern = "test[0-9]"))

dat.list <- 
  lapply(files,
         ingest_xle,
         header.info = FALSE,
         collapse.timestamp = TRUE)

dat.list <- 
  lapply(dat.list,
         as.data.table)

baseline <- 
  dat.list[["test1.000.xle"]]

models.list <- 
  lapply(dat.list,
         function(x){data.table(t(coef(lm(x$level_cm ~ baseline$level_cm))))})

models <- 
  rbindlist(models.list, 
            idcol = "input_source")

setnames(models, c("input_source", "intercept", "slope"))

models[, density := as.numeric(gsub("(test|\\.xle)", "", input_source))]

models[, slope_check := slope - 1/density]


plot(intercept ~ density, data = models)
points(models$density, 
       fitted(lm(intercept ~ density, data = models)),
       col = "red",
       type = "l")
plot(slope ~ density, data = models)
points(models$density, 
       fitted(lm(slope ~ density, data = models)), 
       col = "red",
       type = "l")

dat <- 
  rbindlist(dat.list)

dat[, density := as.numeric(gsub("(test|\\.xle)", "", input_source))]
