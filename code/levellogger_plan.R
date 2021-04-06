
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

    # Load and compensate case study data
    case_study =
      load_case_study(data = file_in("data/case_study/example_data.csv")) %>% 
      correct_data(models = final_models) %>% 
      compensate_data(calibration.data = file_in("data/case_study/calibration_data.csv")) %>%
      smooth_data(n = 13, "raw_compensated_level_cm", "corrected_compensated_level_cm") %>%
      subset_year(year = 2018) %>%
      add_met_data(file_in("../Climate_Change_Impacts/_targets/objects/water_budget")) %>% 
      calculate_sy(),

    # Calculate water balance components
    water_balance =
      case_study %>% 
      calculate_detrended_g() %>%
      calculate_delta_s() %>%
      calculate_et() %>% 
      adjust_water_balance() %>% 
      drop_trailing_data("2018-08-23"),

# Figures -----------------------------------------------------------------

    fig_drivers_panel = 
      create_drivers_panel(combined_data),
    
    fig_bootstrap_timeseries = 
      create_bootstrap_timeseries(fitted,
                                  bootstrap_models,
                                  exp = "var-dis", 
                                  water.sn = "1062452", 
                                  baro.sn = "1066019"),
    
    fig_coefficients_panel = 
      create_coefficients_panel(bootstrap_models),

    case_study_panel = 
      create_case_study_panel(case_study, 
                              file_out("output/figures/case_study_panel.pdf"))


# Tables ------------------------------------------------------------------


  )
