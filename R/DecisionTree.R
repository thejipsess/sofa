train_decision_tree <- function(x, rec = NULL, depth = 10, min_n = 1,
                                complexity = 0.1, save = FALSE,
                                save_name = "Decision_Tree"){
  # Apply recipe
  if(!is.null(rec)){
    x <- prep(rec, new_data = x) %>%
      juice()
  }
  
  # Define the logistic regression model specification
  model_spec <- decision_tree(cost_complexity = complexity,
                              min_n = min_n,
                              tree_depth = depth) %>%
    set_engine("rpart") %>% # Set which engine to use
    set_mode("classification") %>% # Set model to classification mode
    translate() # Convert generic parameters into engine specific parameters.
  
  # Train the model to predict y
  model <- model_spec %>%
    fit(y ~ ., x)
  
  # Locally save model if requested
  if(save){
    saveRDS(model, file = paste("Objects/", save_name, ".Rds", sep = ""))
  }
  
  return(model)
}


tune_decision_tree = function(x, dt_rec, depth = NULL, min_n = NULL,
                              complexity = NULL, k = 10, seed = 123, rep = 1,
                              parallel_comp = TRUE,
                              verbose = TRUE, save = FALSE,
                              save_name = "Decision_Tree_tune_res"){
  tryCatch({
    # Remove duplicate entries
    depth <- unique(depth)
    min_n <- unique(min_n)
    complexity <- unique(complexity)
    
    # Create the grid with possible values
    grid <- expand_grid(tree_depth = depth, min_n = min_n,
                        cost_complexity = complexity, .name_repair = "universal")
    
    # Define the logistic regression model
    # Define the logistic regression model specification
    model_spec <- decision_tree(cost_complexity = tune::tune(),
                                min_n = tune::tune(),
                                tree_depth = tune::tune()) %>%
      set_engine("rpart") %>% # Set which engine to use
      set_mode("classification") %>% # Set model to classification mode
      translate() # Convert generic parameters into engine specific parameters.
    
    # Create the workflow for hyperparameter tuning
    tune_wf <- workflow() %>%
      add_recipe(dt_rec) %>%
      add_model(model_spec)
    
    # Setup k-fold crossvalidation
    k_fold <- vfold_cv(x, v = k, repeats = rep, strata = "y")
    
    # Start parallesisation cluster if enabled
    if(parallel_comp){
      startParallelisation_old()
      registerDoRNG(seed = seed) # Set parallel seed for reproducibilty
    }
    
    # If verbose, print the estimated parallelised computation time.
    if(verbose){
      tic()
      invisible(train_decision_tree(juice(prep(dt_rec, x))))
      dtime <- toc(quiet = TRUE)
      dtime <- dtime$toc - dtime$tic
      etime <- dtime * nrow(grid) / parallel::detectCores(logical = TRUE) + 24
      message(sprintf("Estimated parallelised computation time: %1.0f seconds",
                      etime))
    }
    
    # Perform the hyperparameter tuning
    tic()
    tune_res <- tune::tune_grid(
      tune_wf,
      resamples = k_fold,
      grid = grid,
      metrics = metric_set(yardstick::bal_accuracy, yardstick::accuracy,
                           yardstick::roc_auc, yardstick::precision,
                           yardstick::recall, yardstick::sensitivity,
                           yardstick::specificity, yardstick::mn_log_loss,
                           gain_capture, yardstick::f_meas, yardstick::mcc))
    
    # Generate final message
    dtime <- toc(quiet = TRUE)
    dtime <- (dtime$toc - dtime$tic)/60
    message(sprintf("Hyperparameter tuning finished - took %1.0f seconds", dtime))
    
    # Record status
    status <- list("status" = "Succesful",
                   "dtime" = dtime,
                   "error" = NA)
    
  }, error = function(err){
    # If there was an error during the gridsearch, record this in the status
    print(paste("Gridsearch crashed with error:",err))
    status <- list("status" = "Failed",
                   "dtime" = NA,
                   "error" = err)
    
    # Stop parallel cluster
    if(parallel_comp)
      parallel::stopCluster(cl)
    
    # Return status
    return(list("tune_res" = NA,
                "status" = status))
  })
  
  
  
  
  # Stop the parallel cluster
  if(parallel_comp)
    stopCluster(cl)
  
  
  # Put tuning results and status in 
  return_list <- list("tune_res" = tune_res,
                      "status" = status)
  
  # Locally save tune results if requested
  if(save){
    saveRDS(return_list, file = paste("Objects/", save_name, ".Rds", sep = ""))
  }
  
  # Return the results
  return(return_list)
}
