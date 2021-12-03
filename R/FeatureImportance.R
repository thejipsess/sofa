#' Creates a shapley summary plot
#'
#' @param train dataframe containing the data on which to train the XGBoost model.
#' @param label Vector containing labels of the dependent variable. Does not
#' have to be supplied if dependant variable is present in 'train'.
#' @param params optional dataframe containing hyperparameter values to supply
#' to the XGBoost model. Should be of type dataframe with the following columns:
#' min_n, tree_depth learn_rate and loss_reduction.
#' @param n_features integer specifying the number of top n features to return.
#' @param save Boolean indicating whether to save the plot locally
#' @param save_name String specifying the name of the file to save the plot to
#' if save=TRUE.
#' @param mod If the model has already been trained you can supply it here
#'
#' @return list containing the shapley values, the shapley information in long
#' format, and the ggplot object.
#' @export
#'
#' @examples
plot_shaply_summary <- function(train, mod = NULL, label = NULL, params = NULL,
                                n_features = 10, save = FALSE,
                                save_name = "Shapley summary plot"){
  # Determine first and second label of dependent variable
  label_0 <- names(table(train$y))[1]
  label_1 <- names(table(train$y))[2]
  
  
  if(is.null(params)){
    param_list <- list(objective = "binary:logistic",
                       eval_metric = "logloss")
  } else {
    param_list <- list(objective = "binary:logistic",
                       min_child_weight  = params$min_n,
                       max_depth = params$tree_depth,
                       eta = params$learn_rate,
                       gamma = params$loss_reduction,
                       eval_metric = "logloss")
  }
  
  # If dependent variable is present in train, then remove it
  if('y' %in% colnames(train)){
    label <- train$y
    train <- select(train, -y)
  }
  
  # Make sure train is a matrix
  train <- as.matrix(train)
  
  # Make sure label consists of zeros and ones
  if(class(label) != "numeric")
    label <- as.numeric(label)
  step <- 0
  while(max(label) > 1){
    label <- label - 1
    step <- step + 1
    if(step>100) stop("Something seems to be wrong with the supplied label")
  }
  
  if(is.null(mod)){
    print("Model not supplied so training it now")
    mod <- xgboost::xgboost(data = train,
                          label = as.matrix(label),
                          params = param_list, nrounds = 100,
                          verbose = FALSE, nthread = parallel::detectCores())
    
    # To return the SHAP values and ranked features by mean|SHAP|
    shap_values <- shap.values(xgb_model = mod, X_train = train)
    
    # To prepare the long-format data:
    shap_long <- shap.prep(xgb_model = mod, top_n = n_features,
                           X_train = train)
  } else {
    # To return the SHAP values and ranked features by mean|SHAP|
    shap_values <- shap.values(xgb_model = mod$fit, X_train = train)
    
    # To prepare the long-format data:
    shap_long <- shap.prep(xgb_model = mod$fit, top_n = n_features,
                           X_train = train)
    }
  
  
  # **SHAP summary plot**
  plt <- shap.plot.summary(shap_long) +
    ylab(bquote("Higher predicted probability towards" ~ .(label_0) %<->% ~ "Higher predicted probability towards"~ .(label_1)))
  
  if(save)
    dev.print(svg, paste("Figures/", save_name, ".svg", sep = ""),
              width = 14, height = 7)
  
  return(list(shap_values = shap_values,
              shap_long = shap_long,
              plt = plt))
}

# Get shapley information of top predictors
#
#
#
#   simplify: Whether to simplify the output such that only the mean results
#     per feature are returned rather then results for each sample.
get_shaply_info <- function(train, label = NULL, params = NULL, n_features = 10,
                            simplify = FALSE){
  if(is.null(params)){
    param_list <- list(objective = "binary:logistic",
                       eval_metric = "logloss")
  } else {
    param_list <- list(objective = "binary:logistic",
                       min_child_weight  = params$min_n,
                       max_depth = params$tree_depth,
                       eta = params$learn_rate,
                       gamma = params$loss_reduction,
                       eval_metric = "logloss")
  }
  
  
  # If dependent variable is present in train, then remove it
  if('y' %in% colnames(train)){
    label <- train$y
    train <- select(train, -y)
  }
  
  # Make sure train is a matrix
  train <- as.matrix(train)
  
  # Make sure label consists of zeros and ones
  if(class(label) != "numeric")
    label <- as.numeric(label)
  step <- 0
  while(max(label) > 1){
    label <- label - 1
    step <- step + 1
    if(step>100) stop("Something seems to be wrong with the supplied label")
  }
  
  mod <- xgboost::xgboost(data = train, 
                          label = as.matrix(label), 
                          params = param_list, nrounds = 100,
                          verbose = FALSE, nthread = parallel::detectCores())
  
  mod <- train_xgboost(train, ml_rec, save = F,
                       mtry = show_best(tune_res_xg$tune_res,
                                        metric = metric, n=1)$mtry,
                       trees = show_best(tune_res_xg$tune_res,
                                         metric = metric, n=1)$trees,
                       min_n = show_best(tune_res_xg$tune_res,
                                         metric = metric, n=1)$min_n,
                       tree_depth = show_best(tune_res_xg$tune_res,
                                              metric = metric,
                                              n=1)$tree_depth,
                       learn_rate = show_best(tune_res_xg$tune_res,
                                              metric = metric,
                                              n=1)$learn_rate,
                       loss_reduction = show_best(tune_res_xg$tune_res,
                                                  metric = metric,
                                                  n=1)$loss_reduction)
  
  # To return the SHAP values and ranked features by mean|SHAP|
  shap_values <- shap.values(xgb_model = mod, X_train = train)
  
  # To prepare the long-format data:
  shap_long <- shap.prep(xgb_model = mod, top_n = n_features,
                         X_train = train)
  
  # Simplify result to one row per feature if desired
  if(simplify == T){
    shap_long <- select(shap_long, c(variable, mean_value))
    shap_long <- subset(shap_long, subset = !duplicated(shap_long))
  }
  
  return(shap_long)
}
