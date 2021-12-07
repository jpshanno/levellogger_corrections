
plan <- 
  drake_plan(

# Original Analysis -------------------------------------------------------

    # Read in manual measurements of cable length and water depth
    levellogger_measurements = 
      read_manual_measurements(file_in("data/levellogger_measurements.csv")),
    
    # Load xle data from water transducers
    raw_water_data =
      read_xle_data(levellogger_measurements, 
                    "water"),
    
    # Load xle data from barometric transducers
    raw_barometric_data =
      read_xle_data(levellogger_measurements, 
                    "air"),
    
    # Extract and save metadata for loggers
    instrument_metadata =
      get_metadata(raw_barometric_data, 
                   raw_water_data),
     
    # Combine logger data
    combined_data = 
      combine_datasets(raw_water_data, 
                       raw_barometric_data, 
                       levellogger_measurements, 
                       file_in("data/experimental_design.csv")),
                       
    # Fit 2-stage bootstrap models
    bootstrap_models = 
      fit_models(combined_data),
     
    # Calculated fitted values and bounds for each experimental period
    fitted = 
      calculate_fitted_values(bootstrap_models, combined_data),
    
    # Calculated predicted values for testing data with models built on training
    # data
    predicted =
      calculate_predicted_values(bootstrap_models, combined_data),
    
    # Export models for use
    final_models = 
      export_correction_models(models = bootstrap_models,
                               instrument.metadata = instrument_metadata,
                               file.out = file_out("output/tabular/final_correction_model_coefficients.csv")),


# Case Study --------------------------------------------------------------
# Case study is site 152
    # Load and compensate case study data
    case_study =
      load_case_study(data = file_in("data/case_study/example_data.csv")) %>% 
      correct_data(models = final_models) %>% 
      compensate_data(calibration.data = file_in("data/case_study/calibration_data.csv")) %>%
      smooth_data(n = 13, "raw_compensated_level_cm", "corrected_compensated_level_cm") %>%
      subset_year(year = c(2012, 2018)) %>%
      add_met_data(file_in("data/case_study/water_budget.fst")) %>% 
      # calculate_sy() now handles doing the final subset of data, using 2012
      # for Esy and returning Esy for the case study
      calculate_sy(),

    # Calculate water balance components
    water_balance =
      case_study %>% 
      calculate_detrended_g() %>%
      calculate_delta_s() %>%
      calculate_et() %>% 
      adjust_water_balance(),

    daily_water_balance =
      water_balance %>% 
      expand_instantaneous_rates(reg.exp = "(in|dwt|et)_cm_s") %>% 
      summarize_by_day(mean = c("corrected_compensated_level_cm",
                                "raw_compensated_level_cm",
                                "air_temperature_c",
                                "water_temperature_c"),
                       sum = c("raw_et_cm_15m",
                               "corrected_et_cm_15m",
                               "external_raw_et_cm_15m",
                               "external_corrected_et_cm_15m"),
                       first = c("instrument_error_cm", 
                                 "total_input_cm_d",
                                 "pet_cm_d", 
                                 "tmin_c",
                                 "tmax_c",
                                 "melt_cm",
                                 "rain_cm",
                                 "ytd_water_balance")) %>% 
      transform(dry_day = frollapply(total_input_cm_d, 
                                     n = 2, 
                                     FUN = max,
                                     align = 'right', 
                                     na.rm = TRUE) < 0.1,
                error_cm = raw_compensated_level_cm - corrected_compensated_level_cm),

# Figures -----------------------------------------------------------------

    fig_drivers_panel = 
      create_drivers_panel(combined_data,
                           out.path = "output/figures/Figure_1-Sources_of_Error_Scatter_Plot_Panel",
                           tiff_dpi = 600),
    
    fig_ols_boostrap_comparison = 
      create_bootstrap_ols_comparison(raw.data = combined_data$training, 
                                      models = bootstrap_models,
                                      baro.sn = "1066019",
                                      water.sn = "2069158",
                                      out.path = "output/figures/Figure_3-Density_Plot_Bootstrap_and_OLS_Coefficients",
                                      tiff_dpi = 600),

    fig_bootstrap_timeseries = 
      create_bootstrap_timeseries(fits = fitted,
                                  models = bootstrap_models,
                                  exp = "var-dis", 
                                  water.sn = "1062452", 
                                  baro.sn = "1066019",
                                  out.path = "output/figures/Figure_4-Bootstrap_Correction_Timeseries_Linegraph_and_Ribbon_Panel",
                                  tiff_dpi = 600),
    
    fig_coefficients_panel = 
      create_coefficients_panel(bootstrap_models,
                                out.path = "output/figures/Figure_5-Bootstrap_Coefficients_Pointrange_Panel",
                                tiff_dpi = 600),

    case_study_panel = 
      create_case_study_panel(case_study, 
                              out.path = "output/figures/Figure_6-Case_Study_Panel",
                              tiff_dpi = 600),

    et_to_pet_panel =
      create_et_to_pet_panel(daily_water_balance,
                             out.path = "output/figures/Figure_7-Scatterplot_and_Smooth_ET_to_PET_Panel",
                             tiff_dpi = 600),

    supplemental_figure_s1 = 
      create_figure_s1(file_name = file_out("output/figures/Supplemental_Figure_S1_Ecosystem_Specific_Yield_Functions_Line_and_Smooth.tiff")),

    supplemental_figure_s2 = 
      create_figure_s2(file_name = file_out("output/figures/Supplemental_Figure_S2-Error_Drivers_by_Experimental_Period.tiff")),
    
    supplemental_figure_s3 = 
      create_figure_s3(file_name = file_out("output/figures/Supplemental_Figure_S3-Water_Temperature_Error_Reversal.tiff")),

# Tables ------------------------------------------------------------------


  )
