library(lme4)
library(rms)
library(pROC)
library(data.table)
library(tidyverse)
library(tidymodels)
library(xgboost)
library(themis)
library(tictoc)
library(ggplot2)
library(SHAPforxgboost)
library(fmsb)
library(svglite)

# Load local scripts
source("R/XGBoost.r")
source("R/Parallel.r")
source("R/Visualisation.r")
source("R/FeatureImportance.r")

#=======================#
#                       #
#   Data preparation    #
#                       #
#=======================#

## Databestand gecreeerd met sofa_data.R inlezen
load("../../Data/sofa_data.Rda")


# Set on what day you want to evaluate the samples
day <- 7

## Persoonslevel regressie om dag 1 en t-dagen verandering sofa score met
## minimale meetfout te bepalen
d <- subset(d, d$dag <= day)


SOFA_variables <- c("PF_low", "vent_mode", "trombocytes", "bilirubine",
                    "vasopressors", "nor", "MAP_low", "GCS", "sedation",
                    "GCS_admission", "dialysis", "creatinine", "urine_output")


df <- d %>%
  select(Record.Id, dag, all_of(SOFA_variables)) %>%
  arrange(Record.Id, dag) %>%
  filter(dag %in% c(1, 2, 3, 4, 5, 6, 7)) %>% 
  group_by(Record.Id) %>% filter(n() == 7) %>% # Only keep sampels with data on both days
  pivot_wider(id_cols = Record.Id, names_from = dag, # Encapsulate all sample information in one row
              values_from = SOFA_variables) %>%
  merge(distinct(d[,c("Record.Id","ICU_mortality")]), by = "Record.Id") %>% # Include mortality variable
  rename("y" = "ICU_mortality") %>% # rename mortality variable as y
  column_to_rownames("Record.Id") %>% # Turn sample ID variable into rownames
  mutate_at(vars(starts_with("vent_mode")), factor, levels = c(1,2,3,4,5,6,7,8),
            labels = c("PC", "PS", "VC", "CPAP", "NAVA", "Intellivent", "Other",
                      "None")) %>%
  mutate_at(vars(starts_with("vasopressors")), factor, levels = c(0, 1),
            labels = c("no", "yes")) %>%
  mutate_at(vars(starts_with("sedation")), factor, levels = c(0, 1),
            labels = c("no", "yes")) %>%
  mutate_at(vars(starts_with("dialysis")), factor, levels = c(0, 1),
            labels = c("no", "yes")) %>%
  mutate_at(vars(starts_with("vasopressors")), factor, levels = c(0, 1),
            labels = c("no", "yes")) %>%
  mutate(y = factor(y))

  # mutate(vent_mode_1 = factor(vent_mode_1,levels = c(1,2,3,4,5,6,7,8),
  #                             labels = c("PC", "PS", "VC", "CPAP", "NAVA",
  #                                        "Intellivent", "Other", "None")),
  #        vent_mode_7 = factor(vent_mode_7,levels = c(1,2,3,4,5,6,7,8),
  #                             labels = c("PC", "PS", "VC", "CPAP", "NAVA",
  #                                        "Intellivent", "Other", "None")),
  #        vasopressors_1 = factor(vasopressors_1, labels = c("no", "yes")),
  #        vasopressors_7 = factor(vasopressors_7, labels = c("no", "yes")),
  #        sedation_1 = factor(sedation_1, labels = c("no", "yes")),
  #        sedation_7 = factor(sedation_7, labels = c("no", "yes")),
  #        dialysis_1 = factor(dialysis_1, labels = c("no", "yes")),
  #        dialysis_7 = factor(dialysis_7, labels = c("no", "yes")),
  #        y = factor(y)) # Change columns into appropriate data types



#=======================#
#                       #
#       XGBOOST         #
#                       #
#=======================#
#!!! still need to implement train/test split beforehand

# Set a random seed
seed <- 7649

# Create necessary directories if they do not yet exist
dir.create("Figures", showWarnings = FALSE)
dir.create("Objects", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)

# Create preprocessing recipe for tuning and training
ml_rec <- recipe(y ~ ., data = df) %>%
  #step_range(all_numeric()) %>% # Min-max normalisation
  step_dummy(all_predictors() & where(is.factor)) %>% # Convert to dummy variables
  themis::step_upsample(y)

# Create preprocessing recipe for testing/predicting
test_rec <- recipe(y ~ ., data = df) %>%
  step_dummy(all_predictors() & where(is.factor))# Convert to dummy variables

# Set some naming to base the exported file names on
tune_name <- "_tuning_sofa_seperated_all_days"
train_name <- "_training_sofa_seperated_all_days"
plot_name <- "_plot_sofa_seperated_all_days"

# Set hyperparameter values to evaluate
trees_val <- c(10, 100,1000, 2000)
mtry_val <- c(1, 10, 20)
min_n_val <- c(0, 1, 2, 5, 10, 20, 40)
tree_depth_val <- c(1, 2, 5, 10, 20)
learn_rate_val <- c(1e-8, 1e-4, 1e-2, 1e-1, 0.5, 1)
loss_reduction_val <- c(1e-10, 1e-5, 0.1, 10)

# Tune XGBoost hyperparamters
tune_res_xg <- tune_xgboost(df, ml_rec, trees_val = trees_val,
                            mtry_val = mtry_val, min_n_val = min_n_val,
                            tree_depth_val = tree_depth_val,
                            learn_rate_val = learn_rate_val,
                            loss_reduction_val = loss_reduction_val, rep = 5,
                            seed = seed, parallel_comp = TRUE, verbose = TRUE,
                            k = 5, save = TRUE, entropy_grid = F,
                            save_name = paste("XGBoost", tune_name, sep = ""))

# Set evaluation metric to choose best hyper parameter values
metric = "roc_auc"

plot_tune_res(tune_res_xg$tune_res, save= T, metric = metric,
              save_name = paste("XGBoost", tune_name, sep = ""))



# Train XGBoost
model_xg <- train_xgboost(df, ml_rec, save = TRUE,
                          save_name = paste("XGBoost", train_name, sep = ""),
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

proba_plot_xg <- plot_proba_truth(model_xg, df, test_rec, type = "violin",
                                  "XGBoost - predicted probabilities",
                                  save = TRUE,
                                  save_name =
                                    paste("XGBoost_predicted_probabilities_",
                                          plot_name, sep = ""))

# Spider/radar chart
radar_chart <- plot_radar(list("XGBoost" = model_xg),
                          test = list("XGBoost" = df),
                          dep_var = df$y, vlcex = 1.5,
                          rec = prep(test_rec),
                          title = "", save = T,legend_offset = c(0.04, -0.05),
                          save_name = paste("radar_chart", plot_name, sep = ""),
                          alpha = 0.2, chance= FALSE)

# Variable importance
shapley_info <- get_shaply_info(juice(prep(test_rec)), simplify = TRUE,
                                select_best(tune_res_xg$tune_res, "mn_log_loss"),
                                n_features = 4)

shapley_summary <- plot_shaply_summary(juice(prep(test_rec)), save = F,
                                       n_features = 10,
                                       select_best(tune_res_xg$tune_res,
                                                   "mn_log_loss"),
                                       save_name = paste("shapley_summary",
                                                         plot_name, sep = ""))

shapley_feature_imp <- ggplot(shapley_info, aes(x = variable, y = mean_value,
                                                fill = mean_value)) +
  geom_bar(stat = "identity") +
  coord_flip() + ylab("Mean absolute variable importance") + xlab("")+
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
ggsave(shapley_feature_imp, filename = paste("XGBoost_shapley_importances_",
                                             plot_name, ".svg", sep = ""))

roc_model_xg <- pROC::roc(juice(prep(test_rec))$y,
                          predict_classprob.model_fit(model_xg, juice(prep(test_rec)))[[1]],
                          ci = TRUE)

ggroc(roc_model_xg, legacy.axes = TRUE) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey",
               linetype="dashed")
