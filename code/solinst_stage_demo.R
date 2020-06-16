# Plot corrected water levels for different baro sources

# Email solinst
# - Are xle pressures already temperture compenstated? Seems to be a strong
#   correlation between temperature and pressure
#   - Should diurnal tmeperature signal be adjusted? Doesnt' seem to help
# - Is there a process I can use to correct my levellogger barometric readings
#   to more accurately reflect barologger readings?

# Dual adjustment of water and air pressure to a daily reference temperature
# does not change results much from just adjusting air pressure. That makes me
# think that the internal temperature compensation is mostly working. But there
# is still a diurnal signal that is reverse of anticipated (wetland water levels
# and streamflow increasing through the day and then dropping at night). It 
# really should be the other way around

library(ggplot2)
library(data.table)
library(ingestr)
library(tydygraphs)


# Load & Prep Data --------------------------------------------------------

# Flume levellogger
fl <-
  as.data.table(ingest_xle("flume_053_2018-09-05.xle"))
  # as.data.table(ingest_xle("flume_053_2016-03-10.xle"))

fl <- 
  fl[, .(sample_time = as.POSIXct(paste(sample_date, sample_time),
                                       format = "%Y/%m/%d %H:%M:%S",
                                       tz = "EST"),
         stage_m = level_m - 9.5,
         flume_temperature_c = temperature_c)]
# For 2016 flume file  
# fl[, .(sample_time = as.POSIXct(paste(sample_date, sample_time), 
  #                                 format = "%Y/%m/%d %H:%M:%S",
  #                                 tz = "EST"),
  #        stage_m = level_ft / 3.28 - 9.5, # Altitude offset
  #        flume_temperature_c = temperature_c)]

setkey(fl, sample_time)

# Wetland level
wl <-
  as.data.table(ingest_xle("well_053_2019-04-19.xle"))
  # as.data.table(ingest_xle("well_053_2015-12-08.xle"))

wl <- 
  wl[, .(sample_time = as.POSIXct(paste(sample_date, sample_time),
                                       format = "%Y/%m/%d %H:%M:%S",
                                       tz = "EST"),
         level_m = level_m - 9.5,
         well_temperature_c = temperature_c)]
  # For 2016 well file  
  # wl[, .(sample_time = as.POSIXct(paste(sample_date, sample_time), 
  #                                 format = "%Y/%m/%d %H:%M:%S",
  #                                 tz = "EST"),
  #        level_m = 0.01*level_cm - 488 / 826,
  #        well_temperature_c = temperature_deg_c)]

setkey(wl, sample_time)

# Load external baro pressure
ex_baro <-
  fread("supplementary_baro_data.csv")

ex_baro[,baro_sample_time := as.POSIXct(baro_sample_time, 
                                        format = "%Y-%m-%d %H:%M:%S", 
                                        tz = "UTC")]

ex_baro[, baro_m := baro_level_m - 9.5]

setkey(ex_baro, baro_sample_time, baro_site)

# Load Solisnt Baro
baro <- 
  as.data.table(ingest_xle("baro_053_2018-09-05.xle"))
  # as.data.table(ingest_xle("baro_053_2016-03-10.xle"))

baro <- 
  baro[,.(baro_site = "053",
          baro_sample_time = as.POSIXct(paste(sample_date, sample_time),
                                        format = "%Y/%m/%d %H:%M:%S",
                                        tz = "EST"),
          baro_m = level_m,
          baro_temperature_c = temperature_deg_c,
          baro_input_source = input_source)]
  # For 2016 well file  
  # baro[,.(baro_site = "053", 
  #         baro_sample_time = as.POSIXct(paste(sample_date, sample_time), 
  #                                       format = "%Y/%m/%d %H:%M:%S",
  #                                       tz = "EST"),
  #         baro_m = level_ft / 3.28 - 1600/3.28/826,
  #         baro_temperature_c = temperature_deg_c,
  #         baro_input_source = input_source)]

setkey(baro, baro_sample_time, baro_site)

# Compare Levellogger Barometric Pressure and External Sources ------------

# Combine external and Solinst baro
all_baro <- 
  rbind(ex_baro, 
        baro,
        fill = TRUE)

setkey(all_baro, baro_sample_time)

# Look at barometric pressure across sites
# Levellogger data matches pretty well with 2 of three RAWS data sites
dcast(all_baro[year(baro_sample_time) == 2018], 
      baro_sample_time ~ baro_site, 
      value.var = "baro_m") %>% 
  dygraph(`053`, E6920, KD25, KIWD)

# Combine barometric, well level, and flume stage data
combined <-
  fl[wl][all_baro, nomatch = NULL, on = c("sample_time" = "baro_sample_time")]

# Compensate water levels based on raw barometric pressure
combined[, `:=`(comp_stage_cm = 100 * (stage_m - baro_m),
                comp_level_cm = 100 * (level_m - baro_m))]

# Look at water levels compensated with each barometric pressure
# Again, levellogger data matches pretty well with 2 of three RAWS data sites
dcast(combined,
      sample_time ~ baro_site,
      value.var = "comp_level_cm") %>% 
  dygraph(`053`, E6920, KD25, KIWD) %>% 
  dyRangeSelector()

dcast(combined,
      sample_time ~ baro_site,
      value.var = "comp_stage_cm") %>% 
  dygraph(`053`, E6920, KD25) %>% 
  dyRangeSelector()

# Compensate Temperature Differences --------------------------------------

# Pull out just Solinst data
solinst <- 
  combined[baro_site == "053", 
           .(sample_time, 
             level_m, 
             stage_m,
             flume_temperature_c,
             well_temperature_c, 
             baro_m, 
             baro_temperature_c,
             comp_stage_cm,
             comp_level_cm)]

# Calculate the temperature difference adjustment factor (derived from 
# Gay-Lussac's Law, P1T1 = P2T2)
solinst[, `:=`(level_temp_adjustment_factor = (baro_temperature_c + 273.15) / (well_temperature_c + 273.15),
               stage_temp_adjustment_factor = (baro_temperature_c + 273.15) / (flume_temperature_c + 273.15))]

# Calculate adjusted compensated level
solinst[, `:=`(adj_comp_level_cm = 100 * (level_m - baro_m * level_temp_adjustment_factor),
               adj_comp_stage_cm = 100 * (stage_m - baro_m * stage_temp_adjustment_factor))]

# Compare base compensation with temperature adjusted compensation
# This does help reduce the strength of the diurnal signal
dygraph(solinst, comp_level_cm, adj_comp_level_cm, main = "Well Level") %>% 
  dyRangeSelector()

dygraph(solinst, comp_stage_cm, adj_comp_stage_cm, # flume_temperature_c, baro_temperature_c,
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

manual[, `:=`(diff_level_cm = comp_level_cm - manual_level_cm,
              adj_diff_level_cm = adj_comp_level_cm - manual_level_cm,
              diff_stage_cm = comp_stage_cm - manual_stage_cm,
              adj_diff_stage_cm = adj_comp_stage_cm - manual_stage_cm)]

# The adjusted offsets show a lower coefficient of variation for water level, but
# higher for stage. Perhaps because the flume water levels vary with a similar
# signal as the air temperature

manual[, .(diff_stage_cm = sd(diff_stage_cm, na.rm = TRUE) / mean(diff_stage_cm, na.rm = TRUE),
           adj_diff_stage_cm = sd(adj_diff_stage_cm, na.rm = TRUE) / mean(adj_diff_stage_cm, na.rm = TRUE),
           diff_level_cm = sd(diff_level_cm, na.rm = TRUE) / mean(diff_level_cm, na.rm = TRUE),
           adj_diff_level_cm = sd(adj_diff_level_cm, na.rm = TRUE) / mean(adj_diff_level_cm, na.rm = TRUE))]

offsets <- 
  manual[, 
         lapply(.SD, mean, na.rm = TRUE), 
         .SDcols = c("diff_stage_cm", "adj_diff_stage_cm", 
                     "diff_level_cm", "adj_diff_level_cm")]

setnames(offsets, c("offset_stage", "offset_adj_stage", 
                    "offset_level", "offset_adj_level"))

# Correct Water Levels ----------------------------------------------------

corrected <- 
  cbind(solinst, offsets)

corrected[, `:=`(corrected_stage_cm = comp_stage_cm - offset_stage,
                 adj_corrected_stage_cm = adj_comp_stage_cm - offset_adj_stage,
                 corrected_level_cm = comp_level_cm - offset_level,
                 adj_corrected_level_cm = adj_comp_level_cm - offset_adj_level)]

# corrected[corrected_stage_cm <0, corrected_stage_cm := 0]
# corrected[adj_corrected_stage_cm <0, adj_corrected_stage_cm := 0]

corrected[manual_readings, 
          `:=`(manual_level_cm = manual_level_cm,
               manual_stage_cm = manual_stage_cm)]

dygraph(corrected,
        corrected_level_cm, 
        adj_corrected_level_cm,
        manual_level_cm,
        main = "Well Level") %>% 
  dySeries("manual_level_cm", 
           pointSize = 5) %>% 
  dyRangeSelector()

dygraph(corrected,
        corrected_stage_cm, 
        adj_corrected_stage_cm,
        manual_stage_cm,
        main = "Flume Stage") %>% 
  dySeries("manual_stage_cm", 
           pointSize = 5) %>% 
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
