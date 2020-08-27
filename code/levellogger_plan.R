
plan <- 
  drake_plan(

# Analysis ----------------------------------------------------------------

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
                               file.out = file_out("output/tabular/final_correction_model_coefficients.csv")),

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
      create_coefficients_panel(bootstrap_models)


# Tables ------------------------------------------------------------------


  )
