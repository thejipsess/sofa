#' Train a logistic regression model
#'
#' @param x dataframe where rows represent sampels and columns represent
#' variables. The dependent variable should be in a column names 'y'.
#' @param rec optional recipe object to perform pre-processing steps.
#' @param penalty float specifying the total amount of regularisation in the
#' model. 
#' @param mixture float specifying a value inbetween 0 and 1 that represents the
#' ratio between l1 lasso, and l1 ridge regularisation. Here, mixture = 1 means
#' pure l1 regularisation, and mixture = 1, means pure l2 regularisation.
#' @param save boolean indicating whether to locally save the model or not.
#' @param save_name String specifying the name of the file to save the model to.
#' Only used if save=true.
#'
#' @return trained logistic regression model
#' @export
#'
#' @examples
logistic_regression <- function(x, rec = NULL, penalty = 0, mixture = 0, save = FALSE,
                                save_name = "LR"){
  # Apply recipe
  if(!is.null(rec)){
    x <- prep(rec, new_data = x) %>%
      juice()
  }
  
  # Define the logistic regression model specification
  model_spec <- logistic_reg(penalty = penalty,
                             mixture = mixture) %>%
    set_engine("glmnet") %>% # Set which engine to use
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


#' Hyper parameter tuning for a logistic regression model
#'
#' @param x dataframe where rows represent sampels and columns represent
#' variables. The dependent variable should be in a column names 'y'.
#' @param lr_rec recipe that specifies the preprocessing steps for x.
#' @param penalty float specifying the total amount of regularisation in the
#' model. 
#' @param mixture float specifying a value inbetween 0 and 1 that represents the
#' ratio between L1 lasso, and L2 ridge regularisation. Here, mixture = 1 means
#' pure L1 regularisation, and mixture = 1, means pure L2 regularisation.
#' @param k integer of how many folds to create in k-fold cross-validation
#' @param seed integer that specifies the random seed to use for reproducibility
#' @param rep integer specifying the number of repeats in the repeated k-fol
#' cross-validation
#' @param parallel boolean indicating whether or not to perform gridsearch in 
#' parallel.
#' @param verbose boolean indicating whether or not to print additional
#' information like the expected computation time.
#' @param save boolean indicating whether to locally save the results or not.
#' @param save_name String specifying the name of the file to save the results
#' to. Only used if save=true.
#'
#' @return tuning results as output of the tune::tune_grid function.
#' @export
#'
#' @examples
tune_logistic_regression = function(x, lr_rec, penalty = 10,
                                    mixture = c(0, 0.5, 1), k = 10,
                                    seed = 123, rep = 1, parallel_comp = TRUE,
                                    verbose = TRUE, save = FALSE,
                                    save_name = "LR_tune_res"){
  tryCatch({
    # Remove duplicate entries
    penalty <- unique(penalty)
    mixture <- unique(mixture)
    
    # Create the grid with possible values
    grid <- expand_grid(penalty, mixture, .name_repair = "universal")
    
    # Define the logistic regression model
    model_spec <- logistic_reg(mode = "classification",
                               penalty = tune::tune(),
                               mixture = tune::tune()) %>%
      set_engine("glmnet") %>% 
      set_mode("classification") %>% 
      translate() # Convert generic parameters into engine specific parameters.
    
    # Create the workflow for hyperparameter tuning
    tune_wf <- workflow() %>%
      add_recipe(lr_rec) %>%
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
      invisible(logistic_regression(x, penalty = 1, mixture = 0.5))
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
