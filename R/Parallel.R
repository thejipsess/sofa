# This script contains function to help with parallelisation of various tasks
# and through so speed up the computation time.

# Start a parallel cluster
startParallelisation_old = function(n_threads = NULL){
  library(parallel)
  library(doParallel)
  # Setup parallel computing
  if(is.null(n_threads)){
    # Detect number of available threads
    n_threads <- parallel::detectCores(logical = TRUE)
  } 
  # Generate one R process for each core
  cl <<- parallel::makePSOCKcluster(n_threads) # <<- ensure that the variable is global
  # Register the paralel backedn for the cluster of processes
  registerDoParallel(cl)
  
  message("Parallel comuting enabled")
}

# Start a parallel cluster
startParallelisation_oldfuture = function(n_threads = NULL){
  # Setup parallel computing
  if(is.null(n_threads)){
    # Detect number of available threads
    n_threads <- parallel::detectCores(logical = FALSE)
  } 
  # Generate one R process for each core
  cl <<- parallel::makePSOCKcluster(n_threads) # <<- ensure that the variable is global
  
  # Initate parallelisation through a future plan
  plan(cluster, workers = cl)
  
  
  message("Parallel comuting enabled")
}

# Start a parallel cluster
startParallelisation1 = function(n_threads = NULL){
  # Setup parallel computing
  
  all_cores <- parallel::detectCores(logical = TRUE)
  
  library(doFuture)
  registerDoFuture()
  cl <- makeCluster(all_cores)
  plan(cluster, workers = cl)
}


# Start a parallel cluster
startParallelisation = function(n_threads = NULL){
  if (is.null(n_threads))
    n_threads <- future::availableCores()
  
  registerDoFuture()
  cl <<- future::makeClusterPSOCK(n_threads)
  plan(cluster, workers = cl)
}
