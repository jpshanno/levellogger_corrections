source("code/levellogger_functions.R")
source("code/levellogger_packages.R")
library(ggplot2)

fits <-
  readRDS(glue::glue(".drake/data/{readLines('.drake/keys/objects/fitted')}.rds"))
models <-
  readRDS(glue::glue(".drake/data/{readLines('.drake/keys/objects/bootstrap_models')}.rds"))
exp <- "var-dis"
water.sn <- "1062452"
baro.sn <- "1066019"

    dat2 <- 
      fits[experiment == exp & baro_sn == baro.sn & water_sn == water.sn]
    
    model_file <- 
      glue("output/tabular/bootstrap_fit_values.df/{water.sn}_{baro.sn}_{exp}.fst")
    
    mods <- 
      read_fst(model_file,
               as.data.table = TRUE)
    
    r2_quartiles <- 
      quantile(models[water_sn == water.sn &
                        baro_sn == baro.sn & 
                        experiment == exp, 
                      adj_r2], 
               seq(0, 1, 0.5),
               type = 1)
    
    quartile_reps <- 
      models[water_sn == water.sn &
               baro_sn == baro.sn & 
               experiment == exp & 
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
    
    empty_plot <- 
      ggplot(
        mods2,
        aes(x = sample_time)
      ) +
      coord_cartesian(expand = FALSE) +
      labs(y = "Water Level", 
           x = "Sample Time") +
      theme_minimal(base_size = 12)
      
    
    base <- empty_plot +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  alpha = 0.3,
                  fill = palette.colors()[["skyblue"]]) +
      geom_line(aes(y = centered_water_level_cm),
                color = NA) +
      geom_line(aes(y = water_depth_cm),
                linetype = "dashed",
                color = palette.colors()[["blue"]])
    
    ggsave(plot = base,
           filename = "~/phd/defense_figures/correction_timeseries_01.png",
           dpi = 300,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
    raw_water_level <- empty_plot +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  alpha = 0.3,
                  fill = palette.colors()[["skyblue"]]) +
      geom_line(aes(y = centered_water_level_cm),
                color = NA) +
      geom_line(aes(y = centered_water_level_cm),
                color = "gray75") +
      geom_line(aes(y = water_depth_cm),
                linetype = "dashed",
                color = palette.colors()[["blue"]])
    
    ggsave(plot = raw_water_level,
           filename = "~/phd/defense_figures/correction_timeseries_02.png",
           dpi = 300,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
    corrected_water_level <-
      empty_plot +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  alpha = 0.3,
                  fill = palette.colors()[["skyblue"]]) +
      geom_line(aes(y = centered_water_level_cm),
                color = NA) +
      geom_line(aes(y = centered_water_level_cm),
                color = "gray75") +
      geom_line(aes(y = orig_rect_water_level_cm),
                color = "gray5") +
      geom_line(aes(y = water_depth_cm),
                linetype = "dashed",
                color = palette.colors()[["blue"]])
    
    ggsave(plot = corrected_water_level,
           filename = "~/phd/defense_figures/correction_timeseries_03.png",
           dpi = 300,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
    bootstrap_plot <- empty_plot +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  alpha = 0.3,
                  fill = "gray80") +
      geom_line(aes(y = centered_water_level_cm),
                color = "gray75") +
      geom_line(aes(y = rect_water_level_cm,
                    color = as.factor(model_rep)),
                show.legend = FALSE,
                size = rel(0.2)) +
      geom_line(aes(y = water_depth_cm),
                linetype = "dashed",
                color = palette.colors()[["blue"]]) +
      geom_line(aes(y = orig_rect_water_level_cm),
                color = "gray5")
    
    ggsave(plot = bootstrap_plot,
           filename = "~/phd/defense_figures/correction_timeseries_04.png",
           dpi = 300,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
    
    subplot <- empty_plot +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  alpha = 0.3,
                  fill = "gray80") +
      geom_line(aes(y = rect_water_level_cm,
                    color = as.factor(model_rep)),
                show.legend = FALSE,
                size = rel(0.2)) +
      geom_line(aes(y = orig_rect_water_level_cm),
                color = "gray5") +
      geom_line(aes(y = water_depth_cm),
                linetype = "dashed",
                color = palette.colors()[["blue"]]) +
      scale_x_datetime(
        limits = c(as.POSIXct("2020-06-10 18:00", tz = "EST5EDT"), 
                   as.POSIXct("2020-06-10 20:00", tz = "EST5EDT")),
        breaks = c(as.POSIXct("2020-06-10 18:30", tz = "EST5EDT"),
                   as.POSIXct("2020-06-10 19:30", tz = "EST5EDT")),
        date_labels = "%b %d %H:%M")

    ggsave(plot = subplot,
           filename = "~/phd/defense_figures/correction_timeseries_05.png",
           dpi = 300,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
    full_error <-
      empty_plot +
        geom_ribbon(aes(ymin = propagated_lower,
                        ymax = propagated_upper),
                    alpha = 0.3,
                    fill = "gray30") +
        geom_ribbon(aes(ymin = instrument_lower,
                        ymax = instrument_upper),
                    fill = "gray90") +
        geom_line(aes(y = rect_water_level_cm,
                      color = as.factor(model_rep)),
                  show.legend = FALSE,
                  size = rel(0.3)) +
      geom_line(aes(y = orig_rect_water_level_cm),
                color = "gray5") +
      geom_line(aes(y = water_depth_cm),
                linetype = "dashed",
                color = palette.colors()[["blue"]]) +
        scale_x_datetime(
          limits = c(as.POSIXct("2020-06-10 18:00", tz = "EST5EDT"), 
                     as.POSIXct("2020-06-10 20:00", tz = "EST5EDT")),
          breaks = c(as.POSIXct("2020-06-10 18:30", tz = "EST5EDT"),
                     as.POSIXct("2020-06-10 19:30", tz = "EST5EDT")),
          date_labels = "%b %d %H:%M")
        
      
    ggsave(plot = full_error,
           filename = "~/phd/defense_figures/correction_timeseries_06.png",
           dpi = 300,
           height = 3.75,
           width = 8.4,
           units = 'in')
      
    
    
    
    child <-
      ggplot(mods2[sample == "yes"], 
             aes(x = sample_time)) +
      geom_ribbon(aes(ymin = propagated_lower,
                      ymax = propagated_upper),
                  fill = "gray65") +
      geom_ribbon(aes(ymin = instrument_lower,
                      ymax = instrument_upper),
                  fill = "gray85") +
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
      labs(y = "Water Level (cm)", 
           x = "Sample Time") +
      scale_y_continuous(breaks = c(27.5, 28.5, 29.5),
                         position = 'right') +
      scale_x_datetime(breaks = c(as.POSIXct("2020-06-10 18:30", tz = "EST5EDT"),
                                  as.POSIXct("2020-06-10 19:30", tz = "EST5EDT")),
                       date_labels = "%b %d %H:%M",
                       position = 'top') +
      theme_minimal(base_size = 12) +
      theme(axis.title = element_blank(),
            # axis.text.y = element_blank(),
            plot.background = element_rect(fill = 'white',
                                           color = 'black'))
    
    fig2 <- 
      main + inset_element(child, left = 0.025, bottom = 0.05, right = 0.65, top = 0.5)
    
    ggsave(plot = fig2,
           filename = paste0(out.path, ".pdf"),
           device = cairo_pdf,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
    ggsave(plot = fig2,
           filename = paste0(out.path, ".tiff"),
           type = "cairo",
           dpi = tiff_dpi,
           height = 3.75,
           width = 8.4,
           units = 'in')
    
  }