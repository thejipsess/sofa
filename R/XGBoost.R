
#' Hyperparameter optimisation for XGBoost
#'
#' @param x dataframe containing the data to train the model on. Rows are
#' samples and columns are variables. Should also contain a column named 'y'
#' that contains the dependent variable.
#' @param xg_rec the recipe describing the prepocessing steps for the data, like
#' handling with imblanced data. If you are unfamiliar with recipes than please
#' check out the tidyverse package.
#' @param trees_val positive integer or array of integers specifying the number of
#' trees to grow.
#' @param mtry_val positive integer or array of integers specifying the number of
#' predictors that will be randomly sampled at each split when creating the
#' tree models.
#' @param min_n_val non-negative integer or array of integers specifying the minimum
#' number of data points in a node that is required for the node to be split further.
#' @param tree_depth_val positive integer or array of integers specifying the
#' maximum depth of the tree.
#' @param learn_rate_val positive float or array of floats specifying the rate at
#' which the boosting algorithm adapts from iteration-to-iteration.
#' @param loss_reduction_val positive float or array of floats specifying the
#' reduction in the loss function required to split further.
#' @param entropy_grid Boolean indicating whether to automatically set the
#' search space by maximising the entropy. Setting this parameter to true will
#' prevent any predefined sets of hyperparameter values to be used. See
#' ?dials::grid_max_entropy documentation for more information.
#' @param k positive integer speciying the number of folds for k-fold cross
#' validation.
#' @param rep integer specifying how often to repeat the cross validation.
#' @param seed integer specifying a random seed for reproducibility. If you want
#' the grid search to be completely random and not reproducible, use seed = NULL.
#' @param parallel_comp boolean indicating whether to paralellise the grid search.
#' @param verbose Boolean indicating whether additional information should be 
#' printed.
#' @param save boolean indicating whether to locally save the results or not.
#' @param save_name String specifying the name of the file to save the results
#' to. Only used if save=true.

#'
#' @return the result of the tuning process.
#' @export
#'
#' @examples
tune_xgboost <- function(x, xg_rec = recipe(y ~., x), trees_val = 1000, mtry_val = 5, min_n_val = 0,
                         tree_depth_val= 10, learn_rate_val = 1e-6,loss_reduction_val = 1e-6,
                         entropy_grid = FALSE, k = 10, rep = 1, seed = NULL,
                         parallel_comp = TRUE, verbose = TRUE,
                         save = FALSE, save_name = "XGBoost_tune_res"){
  tryCatch({
    # XGBoost model specification
    xgboost_model <- 
      parsnip::boost_tree(
        mode = "classification",
        trees = tune(),
        mtry = tune(),
        min_n = tune(),
        tree_depth = tune(),
        learn_rate = tune(),
        loss_reduction = tune()
      ) %>%
      set_engine("xgboost", objective = "binary:logistic",
                 eval_metric = "logloss")
    
    # grid specification
    xgboost_params <- 
      dials::parameters(
        trees(),
        finalize(mtry(), x),
        min_n(),
        tree_depth(),
        learn_rate(),
        loss_reduction()
      )
    
    # Define hyperparameter values to evaluate
    if(entropy_grid){
      set.seed(seed)
      xgboost_grid <- 
        dials::grid_max_entropy(
          xgboost_params, 
          size = 250)
    } else{
      # Remove duplicate entries
      trees_val <- unique(trees_val)
      mtry_val <- unique(mtry_val)
      min_n_val <- unique(min_n_val)
      tree_depth_val <- unique(tree_depth_val)
      learn_rate_val <- unique(learn_rate_val)
      loss_reduction_val <- unique(loss_reduction_val)
      
      # Create the grid with possible values
      xgboost_grid <- expand_grid(trees_val, mtry_val, min_n_val, tree_depth_val,
                                  learn_rate_val, loss_reduction_val,
                                  .name_repair = "universal") %>%
        dplyr::rename(c("trees" = "trees_val", # Change paramter names so they will be recognised
                        "mtry" = "mtry_val",
                        "min_n" = "min_n_val",
                        "tree_depth" = "tree_depth_val",
                        "learn_rate" = 'learn_rate_val',
                        "loss_reduction" = "loss_reduction_val"))
    }
    
    # Create workflow
    xgboost_wf <- 
      workflows::workflow() %>%
      add_recipe(xg_rec) %>%
      add_model(xgboost_model)
    
    # Enable paralellisation
    if(parallel_comp == T){
      startParallelisation_old()
    }
    
    # Set resamples
    set.seed(seed) # Set seed for reproducibility
    resamples <- vfold_cv(x, v = k, repeats = rep, strata = "y")
    
    # If verbose, print the estimated parallelised computation time.
    if(verbose){
      # Create temporary XGBoost model to train
      xgboost_model_temp <- 
        parsnip::boost_tree(
          mode = "classification",
          trees = max(trees_val),
          mtry = max(mtry_val),
          min_n = max(min_n_val),
          tree_depth = max(tree_depth_val),
          learn_rate = max(learn_rate_val),
          loss_reduction = max(loss_reduction_val)
        ) %>%
        set_engine("xgboost", objective = "binary:logistic",
                   eval_metric = "logloss")
      tic()
      invisible(train_xgboost(x, xg_rec))
      dtime <- toc(quiet = TRUE)
      dtime <- dtime$toc - dtime$tic + 3
      etime <- (dtime * nrow(xgboost_grid) / parallel::detectCores(logical = FALSE))/60*rep
      message(sprintf("Parallelised computation time shouldn't be more than %1.0f minutes",
                      etime))
    }
    
    # hyperparameter tuning
  
    tic()
    xgboost_tuned <- tune::tune_grid(
      object = xgboost_wf,
      resamples = resamples,
      grid = xgboost_grid,
      metrics = yardstick::metric_set(bal_accuracy, accuracy,
                                      yardstick::roc_auc, precision, recall,
                                      sensitivity, specificity, mn_log_loss,
                                      gain_capture, f_meas, mcc),
      control = tune::control_grid(verbose = F, parallel_over = "everything"))
    
    # Generate final message
    dtime <- toc(quiet = TRUE)
    dtime <- (dtime$toc - dtime$tic)/60
    message(sprintf("Hyperparameter tuning finished - took %1.0f minutes", dtime))
    
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
  
  # stop parallel cluster
  if(parallel_comp == T){
    stopCluster(cl)
  }
  
  return_list <- list("tune_res" = xgboost_tuned,
                      "status" = status)
  
  # Locally save tune results if requested
  if(save){
    saveRDS(return_list, file = paste("Objects/", save_name, ".Rds", sep = ""))
  }
  
  return(return_list)
}

#' Train an XGBoost model
#'
#' @param x dataframe containing the data to train the model on. Rows are
#' samples and columns are variables. Should also contain a column named 'y'
#' that contains the dependent variable.
#' @param rec optional recipe object to perform pre-processing steps.
#' @param xg_rec the recipe describing the prepocessing steps for the data, like
#' handling with imblanced data. If you are unfamiliar with recipes than please
#' check out the tidyverse package.
#' @param save boolean indicating whether to locally save the model or not.
#' @param save_name String specifying the name of the file to save the model to.
#' Only used if save=true.
#'
#' @return
#' @export
#'
#' @examples
train_xgboost <- function(x, rec = NULL, trees_val = 1000, mtry_val = 5, min_n_val = 0,
                          tree_depth_val= 10, learn_rate_val = 1e-6,
                          loss_reduction_val = 1e-6, save = FALSE,
                          save_name = "XGBoost"){
  # Apply recipe
  if(!is.null(rec)){
    x <- prep(rec, new_data = x) %>%
      juice()
  }
  
  # Initiate xgboost model
  # XGBoost model specification
  model <- 
    parsnip::boost_tree(
      mode = "classification",
      mtry = mtry_val,
      trees = trees_val,
      min_n = min_n_val,
      tree_depth = tree_depth_val,
      learn_rate = learn_rate_val,
      loss_reduction = loss_reduction_val
    ) %>%
    set_engine("xgboost", objective = "binary:logistic",
               eval_metric = "logloss")
  
  model_trained <- model %>%
    fit(formula = y ~ .,
        data = x)
  
  # Locally save model if requested
  if(save){
    saveRDS(model_trained, file = paste("Objects/", save_name, ".Rds", sep = ""))
  }
  
  return(model_trained)
}
