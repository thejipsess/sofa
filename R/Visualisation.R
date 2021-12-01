# This script contains visualisation function

# Ensure there's a folder to save plots in
if(!dir.exists("Figures")){
  message("Making 'Figures' directory to store plots")
  dir.create(file.path("Figures"))}

# Initiate colorblind friendly color pallete
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp2 <- c("#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Set some general theme options for all plots from this file
windowsFonts("Helvetica" = windowsFont("Helvetica"))

set_custom_theme <- function(){
  theme_set(theme_light() +
            theme(plot.title = element_text(size = 20, family = "Helvetica",
                                            face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, family = "Helvetica",
                                               hjust = 0.5),
                  text = element_text(size = 18, family = "Helvetica"),
                  axis.title = element_text(family = "Helvetica"),
                  axis.text.x=element_text(size = 16, family = "Helvetica"),
                  strip.background = element_rect(fill=NA),
                  strip.text = element_text(colour = 'black', face = 'bold',
                                            size = 18, family = "Helvetica"),
                  panel.background = element_rect(fill = "transparent")))
}

set_custom_theme()




#' Create distribution plot of all numeric variables in a dataframe
#'
#' @param x dataframe with rows as samples and columns as variables.
#' @param log_scale boolean indicating whether to plot x-axis on log scale
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_distribution <- function(x, log_scale= FALSE){
  
  # Convert dataframe to long format of all numeric variables
  df_numeric <- x %>%
    select_if(is.numeric) %>%
    melt() %>%
    drop_na()
  
  # Generate ordered boxplot
  p <- ggplot(df_numeric,aes(x=reorder(variable, value), y = value))+
    geom_boxplot() +
    coord_flip() +
    xlab("Variables") + 
    ylab("Values")
  
  if(log_scale)
    p <- p + scale_y_log10()
  
  # Print the plot
  print(p)
  
  return(p)
}




#' Plots independent variables and their correlation with a variable of intertest.
#'
#' @param df dataframe of dimension n*m with n samples and m variables.
#' @param method string indicating which correlation to compute. Possible values
#' are: "pearson", "kendall" and "spearman".
#' @param threshold optional float value indicating the p-value threshold for
#' which to include variables in the plot. Variables with p <= threshold will
#' be used.
#' @param save Boolean indicating whether the plot should be saved locally.
#' @param var string indicating the variable of interest with which the
#' correllation will be calcualated.
#' @param title String indicating the title of plot, which if save==TRUE will
#' also be used to name the file.
#'
#' @return ggplot object
#' @export
#'
#' @examples
plot_correlations <- function(df, method = "spearman", threshold = NULL,
                              save = F, var = 'y',
                              title = "Spearman correlation with asthma"){
  # If the variable of interest is not the dependent variable, then replace
  # the dependent variable with the variable of interest and store the
  # dependent variable in the 'dep_var' column.
  if (!var == "y"){
    df$dep_var <- df$y
    df$y <- df[,var]
    df <- select(df, -all_of(var))
  }
  
  # store df in x
  x <- df
  # Append '999' to all variables except the variable of interest. This will
  # later be used to remove everything in fron of the '999' to prevent
  # redundant variable name repetition after one hot encoding
  for(i in seq(1, ncol(x))){
    if(class(x[,i]) == "factor" & colnames(x)[i] != "y"){
      colnames(x)[i] <- str_c(colnames(x)[i], "999")
    }
  }
  
  # Create and juice recipe to perform one hot encoding.
  x <- recipe(y ~ ., data = x) %>%
    step_dummy(all_predictors(), one_hot = F) %>%
    prep() %>%
    juice() %>%
    as.data.frame() # convert to dataframe
  
  # Remove all character in front of and inclduing "999", to create clear\
  # variable names without redundancy.
  colnames(x) <- gsub(".*999_","",colnames(x))
  
  # Loop over all columns and calculate their correlation with the variable of
  # interest.
  res <- data.frame(row.names = colnames(x))
  for(i in seq(1, ncol(x))){
    # Get correlation coefficient
    res[i,1] <- round(as.numeric(cor.test(as.numeric(x[,i]), as.numeric(x$y),
                                          method = method, exact = F)[4]), 4)*-1
    # Get the p-value
    res[i,2] <- round(as.numeric(cor.test(as.numeric(x[,i]), as.numeric(x$y),
                                          method = method, exact = F)[3]), 2)
  }
  
  # Apply threshold filter if requested
  if(!is.null(threshold)){
    res <- res[res[,2] <= threshold,]
  }
  
  colnames(res) <- c("corr", "p") # Set column names
  res$variable <- row.names(res) # add column specifying the variable
  res <- res[res$variable != 'y',] # Remove the variable of interest
  res$p[res$p > 0.1] <- "0.1+" # Simplify p-values bigger than 0.1 to "0.1+"
  # Simplify p-values between 0.05 and 0.1 to "0.05 ~ 0.1"
  res$p[res$p >= 0.05 & res$p <= 0.1] <- "0.05 ~ 0.1"
  
  # Sort rows on the correlation coefficient
  res <- dplyr::arrange(res, corr) %>%
    mutate(variable = factor(variable, levels=variable))
  
  # Plot the variable correlations
  p <- ggplot(data = res, aes(x = variable, y = corr, fill = p)) +
    geom_bar(stat = 'identity', alpha = 0.8, width = 0.8,
             colour = rgb(0,0,0, 0.5)) +
    scale_fill_brewer(direction = -1, labels = c("< 0.01", "< 0.02", "< 0.03",
                                                 "< 0.04", "< 0.05", "0.05 ~ 0.1",
                                                 "> 0.1")) +
    coord_flip() +
    labs(title = title,
         #subtitle = expr(!!names(table(df$y))[1] %<->% !!names(table(df$y))[2]),
         fill = "p-value") +
    ylab("Correlation coefficient") +
    xlab(NULL)
  
  # Locally save the plot if requested
  if(save) ggsave(plot = p, paste("Figures/", title, ".svg", sep=""),
                  height = 7, width = 14)
  
  return(p)
}

plot_correlations_heatmap <- function(df, method = "spearman",
                                      save = F, var = 'y',
                                      save_name = "Spearman correlation heatmap"){
  
  # store df in x
  x <- df
  
  # Swap dependent variable factor level order
  if("y" %in% colnames(x)){
    x$y <- factor(x$y, levels = c(levels(x$y)[2],
                                  levels(x$y)[1]))
  }
  if("y_asthma_atopic_parents" %in% colnames(x)){
    x$y_asthma_atopic_parents <- factor(x$y_asthma_atopic_parents,
                                       levels = c(levels(x$y_asthma_atopic_parents)[2],
                                                  levels(x$y_asthma_atopic_parents)[1]))
  }
  if("y_eczema_first_2yr" %in% colnames(x)){
    x$y_eczema_first_2yr  <- factor(x$y_eczema_first_2yr ,
                                       levels = c(levels(x$y_eczema_first_2yr )[2],
                                                  levels(x$y_eczema_first_2yr )[1]))
  }
  if("y_eczema_6_7yr" %in% colnames(x)){
    x$y_eczema_6_7yr  <- factor(x$y_eczema_6_7yr ,
                                    levels = c(levels(x$y_eczema_6_7yr )[2],
                                               levels(x$y_eczema_6_7yr )[1]))
  }
  
  
  # Append '999' to all variables except the variable of interest. This will
  # later be used to remove everything in fron of the '999' to prevent
  # redundant variable name repetition after one hot encoding
  for(i in seq(1, ncol(x))){
    if(class(x[,i]) == "factor" & colnames(x)[i] != "y"){
      colnames(x)[i] <- str_c(colnames(x)[i], "999")
    }
  }
  
  # Create and juice recipe to perform one hot encoding.
  x <- recipe( ~ ., data = x) %>%
    step_dummy(all_predictors(), one_hot = F) %>%
    prep() %>%
    juice() %>%
    as.data.frame() # convert to dataframe
  
  # Remove all character in front of and inclduing "999", to create clear\
  # variable names without redundancy.
  colnames(x) <- gsub(".*999_","",colnames(x))
  
  # Convert x to data frame
  x <- as.data.frame(x)
  
  # Compute correlations
  cormat <- round(cor(x, method = method, use = "pairwise.complete.obs"), 2)
  
  # Transform correlation matrix to long format
  melted_cormat <- melt(cormat)
  
  # Create correlation heatmap plot
  p <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    labs(x = "", y = "", fill = paste(method, "correlation")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plot if requested
  if(save) ggsave(plot = p, paste("Figures/", save_name, ".svg", sep=""),
                  height = 14)
    
  
  print(p)
  return(p)
}


# This function plots the hyperparameter tuning results
#
# Arguments:
#   tune_res: A list with information about the hypermarater tuning. This is the
#     output of the tune::tune_grid() function.
#
#   Metric: The metric type you'd like to evaluate the hyperparameters on
#' plots the hyperparameter tuning results
#'
#' @param tune_res tune result object as is generated with the tune package
#' @param metric string indicating the evaluation metric to use
#' @param save boolean indicating whether plot should be saved locally
#' @param save_name String specifying how to name the save file if save=TRUE
#' @param summarise boolean indicating whether to summarise the folds into a
#' mean score. summarise=FALSE will plot each individually calculated metric.
#'
#' @return
#' @export
#'
#' @examples
plot_tune_res <- function(tune_res, metric = "roc_auc", save = FALSE,
                          save_name = "plot_tune_res", summarise = TRUE,
                          title = "hyperparameter tuning results"){
  # Plot model's performance with the different hyperparameter values
  if(summarise == TRUE){
    tune_plot_data <- tune_res %>%
      tune::collect_metrics(summarize = summarise) %>%  # Collect tuning results
      filter(.metric == metric) %>%               # Only keep desired metric
      select(-.metric, -.estimator, -n,           # remove unwanted columns
             -std_err, -.config) %>%
      pivot_longer(-mean,                    # Convert wide to long format
                   values_to = "value",
                   names_to = "parameter") %>%
      mutate(value = factor(value))         # Convert value to type factor
  } else {
    tune_plot_data <- tune_res %>%
      collect_metrics(summarize = summarise) %>%     # Collect tuning results
      filter(.metric == metric) %>%               # Only keep desired metric
      select(-id, -id2, -.metric, -.estimator,           # remove unwanted columns
            -.config) %>%
      pivot_longer(-.estimate,                    # Convert wide to long format
                   values_to = "value",
                   names_to = "parameter") %>%
      rename("mean" = ".estimate") %>%  # Rename .estimate to mean to be consistent
      mutate(value = factor(value))          # Convert value to type factor
  }
  
    
    
    tune_plot <- ggplot(tune_plot_data, aes(value, mean, color = value)) +   # Plot the 
      geom_boxplot(show.legend = FALSE,
                   outlier.shape = 4,
                   outlier.stroke = 2,
                   outlier.color = 'red',
                   outlier.size = 2) +
      geom_point(size = 3, alpha = 0.25, show.legend = FALSE) +
      ggtitle(title) +
      facet_wrap(~parameter, scales = "free_x") +
      labs(x = "Hyperparameter value", y = metric) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Show plot
  print(tune_plot)
  
  
  # Save plot if requested
  if(save) ggsave(plot = tune_plot, paste("Figures/", save_name, ".svg", sep=""),
                  height = 7)
}


#' Title
#'
#' @param pca prcomp object as returned by the prcomp function
#' @param y optional array of labels to color the score plot with
#' @param PC_x integer indicating principal component to use on x-axis
#' @param PC_y integer indicating principal component to use on y-axis
#' @param PC_z integer indicating principal component to use on z-axis
#' @param title string containing title of the plot. Also used to name file when
#' save=TRUE
#' @param save boolean indicating whether to locally save the plots
#' @param include_3D boolean indicating whether to include a 3D score plot
#' @param alpha float specifying the degree of transparancy of the dots in the 
#' score plot.
#' @param legend_title string specifying the title to put above the legend of
#' the score plot.
#'
#' @return
#' @export
#'
#' @examples
plot_pca <- function(pca, y = NULL, PC_x = 1, PC_y = 2, PC_z=3, title = "",
                     save=FALSE, include_3D = FALSE, alpha = 1, legend_title = 'legend'){
  # This function take the pca variable generetad by prcomp() and visualises it.
  # By default it plots a histogram of the variance percentage of the PCs
  # Also, it plots a PCA plot of the specified PCs.
  
  # Compute the variance (percentage) of the PCs
  pca_var <- pca$sdev^2
  pca_var_per <- data.frame(round(pca_var/sum(pca_var)*100, 1)) %>%
    rownames_to_column()
  colnames(pca_var_per) <- c('PC', 'variance')
  
  # Plot histograms of the PCs and their variance percentage
  ggplot(pca_var_per, aes(x = reorder(PC, -variance), y = variance)) +
    geom_bar(stat="identity") +
    xlab("Principal component")
  
  # Extrapolate the PCs of interest
  if(!is.null(y)){
    pca.data <- data.frame(Sample=y, X = pca$x[,PC_x], Y = pca$x[,PC_y])
  } else{
    pca.data <- data.frame(X = pca$x[,PC_x], Y = pca$x[,PC_y])
  }
    
  # Plot the PCs of interest
  if(!is.null(y)){
    PCA_score_plot <- ggplot(data = pca.data, aes(x=X, y=Y, label = Sample,
                                                colour = factor(Sample)))
  } else{
    PCA_score_plot <- ggplot(data = pca.data, aes(x=X, y=Y))
  }
  PCA_score_plot <- PCA_score_plot +
    geom_point(alpha = alpha) +
    xlab(paste(colnames(pca$x)[PC_x], " - ", pca_var_per[PC_x, 2], "%", sep="")) +
    ylab(paste(colnames(pca$x)[PC_y], " - ", pca_var_per[PC_y, 2], "%", sep="")) +
    ggtitle(paste("PCA score plot: ", title)) +
    theme(axis.title = element_text(face = "plain", family = "Helvetica")) +
    coord_equal() # Make sure x and y axis are on the same scale
  
  # Add legend title
  PCA_score_plot <- PCA_score_plot + scale_colour_discrete(name = legend_title)
  
  print(PCA_score_plot)
  
  if(save)
    ggsave(paste("Figures/PCA_score_plot_", title, ".svg"),
           plot = PCA_score_plot, height = 7)
  
  if (include_3D){
    pca3d(pca, group=pca.data$Sample, components = c(PC_x, PC_y, PC_z))
  }
}

# MAKE sure no independant variables start with the letter "y"!!
pca_projection <- function(train, test, recipe, train_imputed = NULL,
                           test_imputed = NULL, save = FALSE, transp = FALSE,
                           title = "PCA score plot projection"){
  # Check if boith sets have identical columns
  if(!identical(colnames(train), colnames(test))){
    stop("There is a missmatch between the columns of the inputted train and 
         test set! Both sets require identical columns.")
  }
  
  if(is.null(train_imputed)){
    train$impute <- "unknown"
    test$impute <- "unknown"
  } else{
    train$impute <- train_imputed
    test$impute <- test_imputed
    
  }
  # Create one variable containing imputation information for train and test set
  imputed <- rbind(matrix(train_imputed), matrix(test_imputed))
  
  
  # This function performs pca on a training set, plots it, and projects the set set on top of the PCA
  # Prepare the recipe
  recipe <- prep(recipe, train)
  
  test_norm <- bake(recipe, test) %>%
    select(!starts_with(c("y", "impute")))
  train_norm <- bake(recipe, train) %>%
    select(!starts_with(c("y", "impute")))
  
  # Perform PCA on training set
  pca <- prcomp(train_norm, center = FALSE, scale=FALSE)
  train_scores <- data.frame(pca$x)
  train_scores$type <- rep('train', nrow(pca$x))
  
  # Project test set onto the PCA of the training set
  test_proj <- data.frame(as.matrix(test_norm) %*% pca$rotation)
  test_proj$type <- rep('test', nrow(test_proj))
  
  # Combine PcA and projection in one variable
  pca_proj <- rbind(train_scores, test_proj)
  
  # Add column with imputation information
  
  
  p <- ggplot(data = data.frame(Sample=pca_proj$type, X = pca_proj[,1],
                                Y = pca_proj[,2], impute = imputed),
              aes(x=X, y=Y,fill = factor(Sample), colour = factor(Sample),
                  stroke = 1, shape = impute)) +
    geom_point(alpha = 0.5, size = 4) +
    labs(x = "PC1", y = "PC2", fill = "data split") +
    scale_fill_discrete() +
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    scale_shape_manual(name = "imputation", values = c(24, 21)) +
    guides(colour=FALSE) +
    ggtitle(title) +
    coord_equal() # Make sure x and y axis are on the same scale
  
  # Make background transparent if requested
  if(transp) p <- p + theme(plot.background = element_rect(fill = "transparent",
                                                           color = NA))
  
  #show plot
  print(p)
  
  # Save plot if requested
  if(save){
    ggsave(p, filename = paste("Figures/", title, ".svg", sep = ""),
           bg = "transparent", height = 7)
  }
}


plot_barplot <- function(data, x = 'index', y='data',
                         x_name = "name", y_name="value",
                         title = "barplot", save = FALSE){
  # Convert to dataframe
  if(class(data) != "data.frame")
    data <- data.frame(data)
  
  # Set index as column
  data$index <-rownames(data)
  
  # Sort the dataframe on the y column
  data <- data %>% arrange(!!sym(y)) # Here y is a astring and arrange() automatically quotes inputes meaning that you cannot input variables as they will be quoted. Here I convert the string value of y into symbols with sym() and unquote it with '!!'
  
  # Plot barplot of data
  p <- ggplot(data=data,
                aes(x = reorder(!!sym(x),desc(!!sym(y))), y = !!sym(y),
                    label = round(!!sym(y), digits = 1))) +
                geom_bar(stat="identity", alpha = 0.8, width = 0.8,
                         colour = rgb(0,0,0, 0.5), size = 1) +
                xlab(x_name) +
                ylab(y_name) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                scale_y_continuous(expand = expansion(add = c(0, .5))) +
                ggtitle(title)
  
  if(save){
    ggsave(p, filename = paste("Figures/", title, ".png", sep = ""))
  }
  return(p)
}


plot_missings <- function(data, variables, title = "Percentage of missings",
                          save = FALSE){
  data <- data %>%
    select(variables)
  
  calc_missings <- function(x){
    return((sum(is.na(x)))/length(x)*100)
  }
  
  # Calculate percentage of missings per column
  col_mis <- apply(data, 2, FUN = calc_missings)
  # Calculate percentage of missings per row
  row_mis <- apply(data, 1, FUN = calc_missings)
  
  # PLot a barplot of the percentages of missings per column
  p1 <- plot_barplot(col_mis, x_name = "Variable", y_name="Percentage of missing values",
                     title = paste(title," (per variable)")) +
                      geom_text(angle = 0, nudge_y = .5) +
                      xlab(NULL)
  
  # PLot a barplot of the percentages of missings per row
  p2 <- plot_barplot(row_mis, x_name = "Sample", y_name="Percentage of missing values",
                     title = paste(title," (per sample)")) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        xlab(NULL)
  
  print(p1)
  print(p2)
  
  if(save){
    ggsave(p1, filename = paste("Figures/", title," (per variable)", ".svg",
                                sep = ""))
    ggsave(p2, filename = paste("Figures/", title," (per sample)", ".svg",
                                sep = ""))
  }
  
  return(list("variable" = p1,
              "sample" = p2))
}


plot_venn_diagram <- function(x, type = "percent", save = FALSE,
                              vars = c("mom_asthma",
                                       "mom_mite",
                                       "mom_pet",
                                       "mom_hay" ), title = "venn diagram"){

    if(length(vars) > 4) stop("Cannot create venn diagram for more than 4 variables")

  
  # Get number of variables
  n_vars = length(vars)
  
  x_mat <- x %>%
    select(all_of(vars)) %>%
    mutate_all(as.numeric) %>%
    drop_na()
  
  # Check intersections of the first variable (v1)
  v1 <- x_mat[,1]
  v1[v1 == 2] <- "v1"
  v1[x_mat[,1] == 2 & x_mat[,2] == 2] <- "v1v2"
  if(n_vars >= 3) v1[x_mat[,1] == 2 & x_mat[,3] == 2] <- "v1v3"
  if(n_vars == 4) v1[x_mat[,1] == 2 & x_mat[,4] == 2] <- "v1v4"
  if(n_vars >= 3) v1[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2] <- "v1v2v3"
  if(n_vars == 4) v1[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,4] == 2] <- "v1v2v4"
  if(n_vars == 4) v1[x_mat[,1] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v3v4"
  if(n_vars == 4) v1[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v2v3v4"
  v1 <- v1[v1 != 1]
  v1 <- make.unique(v1, sep = '_') # make duplicated entries unique
  
  # Check intersections of the second variable (v2)
  v2 <- x_mat[,2]
  v2[v2 == 2] <- "v2"
  v2[x_mat[,1] == 2 & x_mat[,2] == 2] <- "v1v2"
  if(n_vars >= 3) v2[x_mat[,2] == 2 & x_mat[,3] == 2] <- "v2v3"
  if(n_vars == 4) v2[x_mat[,2] == 2 & x_mat[,4] == 2] <- "v2v4"
  if(n_vars >= 3) v2[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2] <- "v1v2v3"
  if(n_vars == 4) v2[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,4] == 2] <- "v1v2v4"
  if(n_vars == 4) v2[x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v2v3v4"
  if(n_vars == 4) v2[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v2v3v4"
  v2 <- v2[v2 != 1]
  v2 <- make.unique(v2, sep = '_') # make duplicated entries unique
  
  # Check intersections of the third variable (v3)
  if(n_vars >= 3){
    v3 <- x_mat[,3]
    v3[v3 == 2] <- "v3"
    v3[x_mat[,1] == 2 & x_mat[,3] == 2] <- "v1v3"
    if(n_vars >= 3) v3[x_mat[,2] == 2 & x_mat[,3] == 2] <- "v2v3"
    if(n_vars == 4) v3[x_mat[,3] == 2 & x_mat[,4] == 2] <- "v3v4"
    if(n_vars >= 3) v3[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2] <- "v1v2v3"
    if(n_vars == 4) v3[x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v2v3v4"
    if(n_vars == 4) v3[x_mat[,1] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v3v4"
    if(n_vars == 4) v3[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v2v3v4"
    v3 <- v3[v3 != 1]
    v3 <- make.unique(v3, sep = '_') # make duplicated entries unique
  }
  
  # Check intersections of the third variable (v3)
  if(n_vars >= 4){
    v4 <- x_mat[,4]
    v4[v4 == 2] <- "v4"
    v4[x_mat[,1] == 2 & x_mat[,4] == 2] <- "v1v4"
    if(n_vars >= 3) v4[x_mat[,2] == 2 & x_mat[,4] == 2] <- "v2v4"
    if(n_vars == 4) v4[x_mat[,3] == 2 & x_mat[,4] == 2] <- "v3v4"
    if(n_vars >= 3) v4[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,4] == 2] <- "v1v2v4"
    if(n_vars == 4) v4[x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v2v3v4"
    if(n_vars == 4) v4[x_mat[,1] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v3v4"
    if(n_vars == 4) v4[x_mat[,1] == 2 & x_mat[,2] == 2 & x_mat[,3] == 2 & x_mat[,4] == 2] <- "v1v2v3v4"
    v4 <- v4[v4 != 1]
    v4 <- make.unique(v4, sep = '_') # make duplicated entries unique
  }
  
  
  if (n_vars == 2) venn_data <- list(v1 = v1, v2 = v2)
  if (n_vars == 3) venn_data <- list(v1 = v1, v2 = v2, v3 = v3)
  if (n_vars == 4) venn_data <- list(v1 = v1, v2 = v2, v3 = v3, v4 = v4)
  
  p <- ggVennDiagram(venn_data, label_alpha = 0, label = type,
                category.names = vars)
  if(save){
    ggsave(p, filename = paste("Figures/", title, ".png", sep = ""))
  }
  
  return(p)
  
}


plot_stacked_variable_description <- function(x, save = FALSE,
                                              title = "stacked_variable_description"){
  
  for (i in seq(1, ncol(x))){
    mat <- as.data.frame(table(x[,i]))
    mat[nrow(mat)+1, 2] <- nrow(x) - sum(mat[,2])
    
    mat$Perc <- mat$Freq/sum(mat$Freq)
    mat$Pos <- cumsum(mat$Perc)
    
    p <- ggplot(mat, aes(fill=mat[,1], y=mat[,2], x=colnames(x[i]))) + 
                  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
                  geom_text(aes(y=mat[,4], label=mat[,1]), vjust=1, 
                            color="white", size=3.5)
    
    if(save){
      ggsave(p, filename = paste("Figures/", title, ".png", sep = ""))
    }
    
    print(p)
    
    return(p)
  }
}

# Create diversity plot of a phyloseq object
#
# Parameters:
#   x = phyloseq object
#
#   agg = an aggregation level if you want to to aggregate the data to a
#     specific taxonomic level (e.g agg = "genus")
#
#   ... = If you want to subset the data, you can specify your criteria here
#        (e.g. sex=="male")
plot_diversity <- function(x, ..., label = "y_asthma", div_type = "shannon",
                           agg = NULL, plot_type = "violin", save = FALSE,
                           title = "alpha diversity",
                           save_name = "alpha diversity"){
  
  # Aggregate taxonomy if requested
  if(!is.null(agg)){
    x <- x %>%
      aggregate_taxa(agg)
  }
  
  # Subset the data on user defined criteria if available
  x <- x %>%
    subset_samples(...)

  # Calculate diversity
  x_div <- diversity(x, div_type)
  x_div$label <- apply(sample_data(x)[,label], 1, function(x) as.character(x))

  # Create the plot
  div_plot <- ggplot(x_div, aes_string(x = 'label', y = div_type,
                                       colour = 'label', fill = "label"))

  if(plot_type == "violin"){
    div_plot <- div_plot +
      geom_violin(size = 2, alpha =0.1)
  }
  else if (plot_type == "boxplot"){
    div_plot <- div_plot + geom_boxplot()
  }
  
  # Improve some visual aspects of the plot
  div_plot <- div_plot +
    scale_y_continuous(limits = c(0, NA), expand = c(0,0))+
    labs(x = "", colour = "") +
    ylab("Shannon diveristy") +
    guides(fill = FALSE) +
    theme(text = element_text(size = 24, family = "Helvetica"))
    

  # show graph
  print(div_plot)
  
  if(save){
    ggsave(div_plot, filename = paste("Figures/", save_name, ".svg", sep = ""),
           width = 10, height = 7)
  }
}


#' Create plot of taxonomic composition for all samples
#'
#' @param x 
#' @param tax_level 
#' @param n_taxa 
#' @param group_by 
#' @param save 
#' @param bar_outline_colour 
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
plot_taxonomic_composition <- function(x, tax_level = "genus", n_taxa = 10,
                                       group_by = "y", save = FALSE,
                                       bar_outline_colour = 'black',
                                       title = "Taxonomic composition",
                                       save_name = "Taxonomic composition"){
  
  composition_plot_list <- x %>%
    comp_barplot(
      tax_level = tax_level,
      n_taxa = n_taxa, group_by = group_by, bar_outline_colour = bar_outline_colour
    )
  
  composition_plot <- patchwork::wrap_plots(
    composition_plot_list,
    ncol = 1, guides = "collect"
  ) &
    labs(x = NULL, y = NULL) & theme(axis.text.x = element_blank(),
                                     legend.position="bottom",
                                     legend.direction='vertical')
  
  
  
  if(save){
    ggsave(composition_plot, width = 14, height = 7,
           filename = paste("Figures/", save_name, ".svg", sep = ""))
  }
  
  print(composition_plot)
}


#' Title
#'
#' @param x phyloseq object
#' @param tax_level string describing taxonomic level to use
#' @param n_taxa number of distinct taxa to show
#' @param save boolean indicating whether to save plot
#' @param title string of the title of the plot
#' @param clean_name 
#' @param save_name 
#' @param anno_colour name of sample_data variable to use for colouring
#' geom_segment annotation ring
#'
#' @return
#' @export
#'
#' @examples
plot_taxonomic_composition_iris <- function(x, tax_level = "genus", n_taxa = 10,
                                            save = FALSE, clean_name = FALSE,
                                            save_name = "iris plot",
                                            title = "taxonomic composition",
                                            anno_colour = NULL){
  tax_renamer <- function(x) identity(x)
  if(clean_name){
    tax_renamer <- function(tax) {
      gsub(".*__","",tax)
    }
  }
  
  # Perform ordination
  ordination_var <- ord_calc(x, method = "PCA")
  
  # Set some theme parameters
  theme_set(theme_light() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8)
    ))
  
  # Create iris plot (circular  stacked bar plot)
  composition_plot <- microViz::ord_plot_iris(ordination_var,
                                              taxon_renamer = tax_renamer,
                                              tax_level = tax_level,
                                              n_taxa = n_taxa,
                                              anno_colour = anno_colour) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank()) +
    labs(fill = tax_level)
  
  # Show plot
  print(composition_plot)
  
  
  
  if(save){
    ggsave(composition_plot, width = 10, height = 8,
           filename = paste("Figures/", save_name, ".svg", sep = ""))
  }
  
  # Reset custom theme for other plots
  set_custom_theme()
  
  return(composition_plot)
}


#' Generates ordination score plots of microbiome data after specifiable
#' pre-processing steps.
#'
#' @param x phyloseq object containing the microbiome data
#' @param pc1 integer specifying which first principal component to use.
#' @param pc2 integer specifying which second principal component to use.
#' @param pc3 integer specifying which third principal component to use. Only
#' used if threeD=TRUE.
#' @param min_prevalence number or proportion of samples that a taxon must be
#' present in
#' @param save boolean indicating whether to locally save the plot
#' @param title string specifying the title of the plot. If save=TRUE the title
#' will also be used to generate the file name.
#' @param tax_level string specifying to which taxonomic level to aggregate.
#' @param transform string specifying how to transform the data. Can take any
#' valid taxa transformation from microbiome::transform.
#' @param label string specifying with which variable to label the score plot.
#' @param threeD boolean indicating whether to also generate a 3D score plot.
#' @param clean_names boolean indicating whether to remove prefixes from the
#' taxonomic names.
#'
#' @return
#' @export
#'
#' @examples
plot_microbiome_ordination <- function(x, pc1 = 1, pc2 = 2, pc3 = 3,
                                       min_prevalence = 0.1, save = FALSE,
                                       title = "microbiome_ordination",
                                       clean_name = FALSE,
                                       tax_level = "genus", transform = "clr",
                                       label = "y", threeD = FALSE,
                                       save_name = "Microbiome ordination"){
  tax_renamer <- function(x) identity(x)
  if(clean_name){
    tax_renamer <- function(tax) {
      gsub(".*__","",tax)
    }
  }
  
  # Create RDA plot
  micro_genus_pca <-
    x %>%
    tax_filter(min_prevalence = min_prevalence, tax_level = tax_level) %>% # Filter out rare taxa
    tax_agg(tax_level) %>% # Aggregate taxa onto genera
    tax_transform(transform) # Transform with requested method
  
  
  micro_genus_pca <- ord_calc(micro_genus_pca, method = "PCA") # Perform PCA
  
  # create plot
  pca_plot <- micro_genus_pca %>%
    ord_plot(
      plot_taxa = 1:6,
      colour = label, 
      tax_vec_length = 0.2,
      auto_caption = F,
      axes = c(pc1, pc2),
      taxon_renamer = tax_renamer,
      tax_lab_style = aes(size = 5),
      size = 4, shape = 21, fill = label, alpha = 0.5, stroke = 1)
  
  # customise plot
  customised_plot <- pca_plot +
    # stat_ellipse(aes_string(colour = label, fill = label),
    #              geom="polygon", alpha=0.05, lwd = 1, segments = 100) +
    labs(fill = "") +
    guides(colour = FALSE) +
    theme(plot.title = element_text(size = 20, family = "Helvetica",
                                    face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 14, family = "Helvetica",
                                       hjust = 0.5),
          text = element_text(size = 18, family = "Helvetica"),
          axis.title = element_text(family = "Helvetica"),
          axis.text.x=element_text(size = 16, family = "Helvetica"),
          strip.background = element_rect(fill=NA),
          strip.text = element_text(colour = 'black', face = 'bold',
                                    size = 18, family = "Helvetica"),
          panel.border = element_rect(colour = rgb(0, 0, 0, 0.5), fill = NA,
                                      size = 0.2)) +
    coord_equal() # Make sure x and y axis are on the same scale
  
  if(save){
    ggsave(customised_plot, height = 7,
           filename = paste("Figures/", save_name, ".svg", sep = ""))
  }
  
  # show plot
  print(customised_plot)
  
  if(threeD){
    micro_genus_pca_scores <- vegan::scores(micro_genus_pca$ord,
                                            choices = c(pc1, pc2, pc3))$sites
    pca3d(micro_genus_pca_scores[,1:3], group = df_imputed[, label])
  }
  
  return(vegan::scores(micro_genus_pca$ord))
}

#' Title
#'
#' @param x 
#' @param label 
#' @param scaler 
#' @param agg 
#' @param prevalence 
#' @param title 
#' @param save 
#'
#' @return
#' @export
#'
#' @examples
plot_taxa_comparison <- function(x, label = "y",
                                 agg = "genus", prevalence = 0.1,
                                 title = "taxa_comparison", save = FALSE){
  
  # Aggregate the phyloseq object
  micro_scale <- x %>%
    aggregate_taxa(agg)
  
  
  df_sample <- data.frame(sample_data(micro_scale))
  df_sample$label <- data.frame(sample_data(micro_scale)[,label])
  
  sample_data(micro_scale) <- df_sample
  
  # Extrapolate all possible values
  
  values <- data.frame(sample_data(micro_scale)$label[,1])
  
  # Get the different unique group labels
  groups <- unique(sample_data(micro_scale)$label[,1])
  
  # Initaite empty list to store results
  list_core <- c()
  
  for (group in groups){ # for each variable group in DiseaseState
    group <<- group
    # Subset current group
    micro_sub <- subset_samples(micro_scale, label == group) 
    
    # Determine the core taxa for the current group
    core_m <- tax_filter(micro_sub, min_prevalence = prevalence,
                         prev_detection_threshold = 0.001) %>%
      taxa_names()
    # print core taxa identified in current group
    print(paste0("No. of core taxa in ", group, " : ", length(core_m)))
    # Store results in the list
    print("group is:")
    print(group)
    list_core[[group]] <- core_m
  }
  
  # show the resulting list of core taxa
  print(list_core)
  
  # Create venn diagram
  p <- ggVennDiagram(list_core) +
    scale_colour_brewer(palette = "Set1")
    
  
  if(save){
    ggsave(p, filename = paste("Figures/", title, ".png", sep = ""))
  }
  
  return(p)
}


#' Title
#'
#' @param x Object that can take two types. Either you input trained model or
#' you input a dataframe. If you input a dataframe, it should contain the
#' predictons of your model where rows represent samples. The columns that should
#' be present are "pred", "proba_1", and "true". The "pred" column should contain
#' the factors of the class that the model predicted for the corresponding
#' samples, "proba_1" should contain the predicted probability of each sample
#' belonging to the first class, and the "true" column should contains factors
#' of the true classes of the samples.
#' @param test dataframe of the test set, including the dependant variable in
#' the "y" column. This parameter should only be set when you inputted a trained
#' model for x.
#' @param rec preprocessing recipe
#' @param title string specifying the title of the plot.
#'
#' @return
#' @export
#'
#' @examples
plot_proba_truth <- function(x, test = NULL, rec = NULL, type = "scatter",
                             title = "predicted probabilities", save = FALSE,
                             save_name = "predicted probabilities"){
  if(!is.data.frame(x)){
    # Set class_one to the class name of class one
    class_one = levels(test$y)[1]
    if(!is.null(rec)){
      test <- bake(prep(rec), new_data = test)
      
    } else{
      # Set class_one to generic string as class one is unknown
      class_one = "class one"
    }
    
    # Predict test set
    test_pred <- predict(x, test)
    test_pred_proba <- predict(x, test, type = "prob")
    
    # Store true and predicted class
    x <- bind_cols(test_pred, test_pred_proba,
                            test$y)
    colnames(x) <- c("pred", "proba_1", "proba_2", "true")
  }
  
  # Order samples on their predicted probability
  x <- x[order(x$proba_1),]
  x$sample <- seq(1, nrow(x))
  
  # check which factor level corresponds to class one
  class_one <- levels(df$y)[1]
  
  # Create plot
  if(type == "scatter"){
    p <-ggplot(x, aes(x = sample, y = proba_1, color = true, alpha = 0.6)) +
      geom_point(size = 5) +
      ylab(paste("Predicted probability of", class_one)) +
      xlab("True class") +
      labs(fill = "True class") +
      theme(legend.title=element_blank()) +
      ggtitle(title)
  } else if(type == "boxplot"){
    p <- ggplot(x, aes(x = true, y = proba_1, group = true, color = true)) +
      geom_boxplot(size = 2) +
      ylab(paste("Predicted probability of", class_one)) +
      xlab("True class") +
      theme(legend.title=element_blank()) +
      labs(fill = "True class") +
      ggtitle(title)
  }else if(type == "violin"){
    p <- ggplot(x, aes(x = true, y = proba_1, group = true, color = true)) +
      geom_violin(size = 2) +
      theme(legend.title=element_blank()) +
      ylab(paste("Predicted probability of", class_one)) +
      xlab("True class") +
      stat_summary(fun.data = mean_sdl, size = 2,
                   geom = "pointrange", alpha = 0.3, fun.args = list(mult = 1))+
      ggtitle(title)
  }
  
  
  # print plot
  print(p)
  
  if(save)
    ggsave(p, filename = paste("Figures/", save_name, ".png", sep = ""))
  
  return(p)
}


#DISCLAIMER: CURRENTLY ONLY WORKS WHEN INPUTTING MULTIPLE MODELS, SINGLE MODEL
# FUNCTIONAILITY IS CURRENTLY BROKEN.

#' Title
#'
#' @param model a list or model. If you enter a single model, it can be any model
#' that can be used by the predict function, e.g. a ranger model object. If you
#' want to plot multiple models simultaneously, you have to enter a list
#' with multiple models.
#' @param test a list or dataframe. if you input one model, this should be one
#' dataframe containing the test data for the model. If you input multiple
#' models, this parameter should be a list, where each ith item is a dataframe
#' of the test set corresponding to the ith model in model.
#' @param dep_var A factor array containing the dependant variable for the test
#' set.
#' @param model_name 
#' @param color 
#' @param save 
#' @param vlcex Font size magnification for vlabels.
#' @param caxislabels Character vector for center axis labels, overwriting
#' values specified in axistype option. If NULL, the values specified by
#' axistype option are used.
#' @param title 
#' @param alpha 
#' @param ... 
#' @param legend_offset 
#' @param legend_horizontal 
#' @param chance Boolean indicating whether to add the 0.5 threshold which
#' indicates how a purely random model would perform.
#' @param rec preprocessing recipe that IS PREPARED on the training set.
#' @param save_name 
#'
#' @return
#' @export
#'
#' @examples
plot_radar <- function(model, test, dep_var, rec = NULL, model_name = "model",
                       color = cbp2,
                       save = FALSE, vlcex = 0.7, caxislabels = NULL,
                       title = "Radar chart", alpha = 0.5, chance = TRUE,
                       legend_offset = c(0,0), legend_horizontal = FALSE,
                       save_name = "Radar chart", ...){
  
  df <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(df) <- c("accuracy", "balanced accuracy", "AUROC", "F1", "precision",
                    "sensitivity", "specificity", "1 - log loss")
  df["max",] <- c(1, 1, 1, 1, 1, 1, 1, 0)
  df["min",] <- c(0, 0, 0, 0, 0, 0, 0, 1)
  if(chance)
    df["chance threshold (0.5)",] <- t(rep(0.5, 8))
  
  # Create variable to count how many models are added
  count <- 1
  
  add_metrics <- function(model, model_name, test, df, i){
    # Apply preprocessing recipe to test set
    if(!is.null(rec))
      test <- bake(rec, new_data = test)
    
    # Predict the test set
    if(class(model)[1] == "mixo_splsda"){ # If model is sPLS-DA prediction works a bit different
      prediction <- predict(model, select(test, -y), dist = "mahalanobis.dist")
      df_test_pred <- as.data.frame(prediction$MajorityVote) %>%
        last() %>% # Pick last column
        factor(levels = levels(dep_var)) # Convert to factor
      df_test_pred_proba <- as.data.frame(prediction$predict) %>%
        nth(-2) # Pick column before the last column
    } else{ # If not an sPLS-DA model, use the regular framework
      df_test_pred <- predict_class.model_fit(model, test)
      df_test_pred_proba <- predict_classprob.model_fit(model, test)[,1]
    }
    
    
    # Store true and predicted class
    truth_pred <- bind_cols(dep_var, df_test_pred, df_test_pred_proba)
    colnames(truth_pred) <- c("truth", "pred", "proba")
    
    df[model_name, "accuracy"] <- yardstick::accuracy(truth_pred, truth, pred)$.estimate
    df[model_name, "balanced accuracy"] <- yardstick::bal_accuracy(truth_pred, truth, pred)$.estimate
    df[model_name, "AUROC"] <- yardstick::roc_auc(truth_pred, truth, proba)$.estimate
    df[model_name, "F1"] <- yardstick::f_meas(truth_pred, truth, pred)$.estimate
    df[model_name, "precision"] <- yardstick::precision(truth_pred, truth, pred)$.estimate
    df[model_name, "sensitivity"] <- yardstick::sensitivity(truth_pred, truth, pred)$.estimate
    df[model_name, "specificity"] <- yardstick::specificity(truth_pred, truth, pred)$.estimate
    df[model_name, "1 - log loss"] <- yardstick::mn_log_loss(truth_pred, truth, proba)$.estimate
    
    # If any log loss score is bigger than one, set it to one
    # This is only done so it can always be plotted with the other metrics which
    # are always withing the 0 - 1 range.
    df$`1 - log loss`[df$`1 - log loss` > 1] <- 1
    
    # Keep track of how many models are added
    count <- count + 1
    
    return(df)
  }
  
  if(class(model)[1] == "list"){
    for(i in seq(1, length(model))){
      model_i <- model[[i]]
      model_name_i <- names(model)[i]
      test_i <- test[[i]]
      df <- add_metrics(model_i, model_name_i, test_i, df, i)
    }
  } else {
    # Apply preprocessing recipe to test set
    if(!is.null(rec))
      test <- bake(rec, new_data = test)
    
    # Calculate metrics for model
    df <- add_metrics(model, test, df, 1)
  }
  
  
  # Create radar chart
  p <- radarchart(
    df, axistype = 1,
    # Customize styling of the polgyon
    pcol = color, pfcol = scales::alpha(color, alpha), plwd = 2, plty = 1,
    # Customize styling of the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Set variable labels
    vlcex = vlcex,
    caxislabels = caxislabels, title = title,
    # Add other arguments if available
    ...
  )
  
  # Add legend
  legend(x = "topright", legend = rownames(df[-c(1,2),]), bty = "n", col = color,
         text.col = "black", cex = 1.5, pt.cex = 3, pch = 20,
         inset = legend_offset, horiz = legend_horizontal, xpd = TRUE)
  
  # save plot if requested
  if(save)
    dev.print(svg, paste("Figures/", save_name, ".svg", sep = ""), width = 15,
              height = 10)
  
  print(p)
  
  return(p)
}
