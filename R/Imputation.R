#' Random Forest imputation
#'
#' @param df 
#' @param seed 
#' @param max_missing maximum percentage of missing data per sample, samples
#' with more missings will be dropped
#' @param variable_max_missing maximum percentage of missing data per variable,
#' variables with more missings will be dropped
#'
#' @return
#' @export
#'
#' @examples
RF_impute <- function(df, seed = 1, max_missing = 20, variable_max_missing= 100){
  # Temporarily remove Sample ID column
  ID <- df$Record.Id
  df <- select(df, -Record.Id)
  
  # Drop samples with too many missing data points
  miss_rate <- apply(df, 1, function(x) sum(is.na(x))/length(x)*100)
  df <- df[miss_rate <= max_missing,]
  ID <- ID[miss_rate <= max_missing]
  
  # Drop variables with too many missing data points
  miss_rate <- apply(df, 2, function(x) sum(is.na(x))/length(x)*100)
  if(any(miss_rate > variable_max_missing))
  message(paste("The following variables will be removed due to excessive missings:",
              names(df[,miss_rate > variable_max_missing])))
  df <- df[,miss_rate <= variable_max_missing]
  
  
  # Determine best paralellisation method and number of threads to use
  if (parallel::detectCores() > ncol(df))
    parallel_method <- "no" else if (parallel::detectCores() > 10)
      parallel_method <- "variable" else
        parallel_method <- "forests"
  if (parallel::detectCores() > ncol(df))
    threads <- ncol(df) else
      threads <- parallel::detectCores()
        
  # Mark for each samples whether it has got missings/to-be-imputed values
  imp_index <- rep(0, nrow(df))
  imp_index[rowSums(is.na(df)) > 0] <- 1
  imp_index <- factor(imp_index, levels = c(0, 1),
                       labels = c("not_imputed", "imputed"))
  
  startParallelisation_old(n_threads = threads) # Enable paralellisation
  registerDoRNG(seed = seed) # Set parallel seed for reproducibility
  df_res <- missForest(df, maxiter = 20, variablewise = T,
                       parallelize = parallel_method)
  df <- df_res$ximp
  df$Record.Id <- ID
  return(list("df" = df, "imputation_index" = imp_index))
}
