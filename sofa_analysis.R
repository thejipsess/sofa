################################################################################
### MaastrICCht cohort
### Auteur: Sander van Kuijk & Frank van Rosmalen
###
### Doel project: Predictie mortaliteit na IC-opname o.b.v. SOFA dag 1-5/7
### Doel syntax: model voorspelling uiteindelijk overlijden
###
### Start: 26/11/2021
### Laatste aanpassing: 30/11/2021
###
### sessionInfo()
###
### R version 4.0.4 (2021-02-15)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 19042)
###
### Matrix products: default
###
### locale:
### [1] LC_COLLATE=English_Netherlands.1252  LC_CTYPE=English_Netherlands.1252
### [3] LC_MONETARY=English_Netherlands.1252 LC_NUMERIC=C
### [5] LC_TIME=English_Netherlands.1252
###
### attached base packages:
### [1] stats     graphics  grDevices utils     datasets  methods   base
###
### loaded via a namespace (and not attached):
### [1] compiler_4.0.4
################################################################################

library(lme4)
library(rms)
library(pROC)
library(data.table)
library(tidyverse)
library(tidymodels)
library(xgboost)
library(glmnet)
library(themis)
library(tictoc)
library(ggplot2)
library(SHAPforxgboost)
library(fmsb)
library(svglite)
library(doRNG)

# Load local scripts
source("R/XGBoost.r")
source("R/LogisticRegression.r")
source("R/Parallel.r")
source("R/Visualisation.r")
source("R/FeatureImportance.r")
source("R/FeatureImportance.r")
source("R/CustomModels.R")


#==============================================================================#
#
#   This Section contains the code from Sander
#
#==============================================================================#

## Databestand gecreeerd met sofa_data.R inlezen
load("../../Data/sofa_data.Rda")

# Set on what day you want to evaluate the samples
day <- 7

## Persoonslevel regressie om dag 1 en t-dagen verandering sofa score met
## minimale meetfout te bepalen
d <- subset(d, d$dag < day)

# Order samples on ID and day
d <- arrange(d, Record.Id, dag)

models <- lmList(SOFA_score ~ dag | Record.Id, data = d, na.action = na.omit)
coef(models)
res <- data.frame(Record.Id = rownames(coef(models)),
                  coef(models), check.names = FALSE)
names(res)[2] <- "Intercept"
names(res)[3] <- "Slope"

res$day1  <- round(res$Intercept + 1*res$Slope, 1)
res$day5  <- round(res$Intercept + 5*res$Slope, 1)
res$delta <- round(4*res$Slope, 1)

## Databestand maken waarbij iedere patient één observatie bijdraagt
dp <- d[!duplicated(d$Record.Id, fromLast = TRUE), ]
dp <- merge(dp, res, by = "Record.Id")

empty <- c("X", "X.y", "X.x")
dp <- dp[, !names(dp) %in% empty]

constant <- c("CHC.Dementia", "CHC.Connective_tissue_disease", "CHC.Hemiplegia",
              "CHC.AIDS", "ckd_status.Creatinine_265_mmolL", "adrenaline",
              "dobutamine", "dopamine", "anti_viral.Isavuconazol", "isOnMediumCare")
dp <- dp[, !names(dp) %in% constant]

## Model om mortaliteit te voorspellen
table(dp$ICU_mortality)
dp$event <- ifelse(dp$ICU_mortality == "Death", 1, 0)
table(dp$event)

## Afgeleiden bepalen
dp$sofa_stijging <- ifelse(dp$delta > 0, 1, 0)
dp$adv_age <- ifelse(dp$age > 65, 1, 0); table(dp$adv_age)

## Selectie?
## dp <- subset(dp, dp$ECMO != 1)

dd <- datadist(dp)
options(datadist = "dd")

## Model enkel op basis van SOFA score dag 1 en delta 1 tot 5
modela <- lrm(event ~ delta, data  = dp, x = TRUE, y = TRUE)
modela
modelb <- lrm(event ~ day1 + delta, data = dp, x = TRUE, y = TRUE)
modelb

## Model aangevuld met geslacht en leeftijd, geslacht discutabel
model2 <-  lrm(event ~ day1 + delta + age, data = dp, x = TRUE, y = TRUE)
model2

model3 <- lrm(event ~ day1 + delta + age + gender, data = dp, x = TRUE, y = TRUE)
model3

## ROC curve
r <- pROC::roc(dp$event, predict(model3, type = "fitted"),  ci = TRUE)

png("auc.png", width = 500, height = 500, pointsize = 16)
ggroc(r, legacy.axes = TRUE) + 
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey",
                     linetype="dashed")
dev.off()

set.seed(070181)
roc_auc <- round(0.5 + 0.5*validate(model2, B = 1000)[1, 5], 2)
roc_auc

## Afkappunten bepalen
construct_decision_matrix(dp, model3, model_type = "Sander") # This function comes from Visualisation.R

## Alternatief: >70 komt er niet in
risk.FN  <- ifelse(dp$age < 75, 0, 1)
sensl <- round(table(risk.FN, dp$event)[2,2]/(table(risk.FN, dp$event)[2,2] +
                  table(risk.FN, dp$event)[1,2])*100, 1)
specl <- round(table(risk.FN, dp$event)[1,1]/(table(risk.FN, dp$event)[1,1] +
                  table(risk.FN, dp$event)[2,1])*100, 1)
ppvl  <- round(table(risk.FN, dp$event)[2,2]/(table(risk.FN, dp$event)[2,2] +
                  table(risk.FN, dp$event)[2,1])*100, 1)
npvl  <- round(table(risk.FN, dp$event)[1,1]/(table(risk.FN, dp$event)[1,1] +
                  table(risk.FN, dp$event)[1,2])*100, 1)
abstl <- sum(risk.FN == 1, na.rm = TRUE)
miscl <- table(risk.FN, dp$event)[2, 1]
ligdl <- sum(ifelse(dp$ICU_LoS[risk.FN == 1] - 7 < 1, 0,
                       dp$ICU_LoS[risk.FN == 1] - 7), na.rm = TRUE)
ligpl <- round(ligdl/sum(dp$ICU_LoS, na.rm = TRUE)*100, 1)

afkapl <- "lft>75"
data.frame(afkapl, sensl, specl, ppvl, npvl, abstl, miscl, ligdl, ligpl)

df_NICE <- NICE_triage(dp)
dm_NICE <- construct_decision_matrix(df_NICE, model_type = "custom", day = 0)

#==============================================================================#
#
#   This Section contains the code from Jip
#
#==============================================================================#

#=======================#
#                       #
#       XGBOOST         #
#                       #
#=======================#
#!!! still need to implement train/test split beforehand

# Set a random seed
seed <- 7649

# Define repeted cross-validation setup for hyperparameter optimisation
k = 5
repeats <- 5

# Set evaluation metric to choose best hyper parameter values
metric = "mn_log_loss"

# Define Feature set to use. Options are "SOFA", "extra" and "SOFA_separated".
feature_set <- "SOFA"

# Create necessary directories if they do not yet exist
dir.create("Figures", showWarnings = FALSE)
dir.create("Objects", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)

# Set rownames to Record.ID
df <- dp
row.names(df) <- df$Record.Id

if(feature_set == "SOFA"){
        df <- df %>%
                select(event, day1, delta, age, gender) %>% # select variables of interest
                drop_na()
} else if(feature_set == "extra"){
        df <- df %>%
                select(event, day1, delta, age, gender, BMI, sepsis,
                       systolicBP_high, temp_high_admission, temp_high, systolicBP_high,
                       systolicBP_high_admission, systolicBP_low_admission, systolicBP_low,
                       ventilation_admission, APACHE_II, temp_low_admission, 
                       temp_low, bilirubine) %>% # select variables of interest
                drop_na()
        
        # Fix some data types
        df$temp_low <- as.numeric(df$temp_low)
        df$temp_high <- as.numeric(df$temp_high)
} else if(feature_set == "SOFA_separated"){
        SOFA_variables <- c("event", "PF_low", "vent_mode", "trombocytes", "bilirubine",
                            "vasopressors", "nor", "MAP_low", "GCS", "sedation",
                            "GCS_admission", "dialysis", "creatinine", "urine_output")

        df <- dp %>%
                select(all_of(SOFA_variables)) %>%
                drop_na()
}

# Convert character columns to factors
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.factor)

# Set dependent variable name to y
names(df)[1] <- "y"

# Factorise dependent variable
df$y <- as.factor(df$y)

# Make sure there are no missings
if(any(is.na(df))) print("WARNING!! There are still incomplete samples present")

# Create preprocessing recipe for tuning and training
ml_rec <- recipe(y ~ ., data = df) %>%
        #step_range(all_numeric()) %>% # Min-max normalisation
        step_dummy(all_predictors() & where(is.factor)) %>% # Convert to dummy variables
        themis::step_upsample(y)

# Create preprocessing recipe for testing/predicting
test_rec <- recipe(y ~ ., data = df) %>%
        step_dummy(all_predictors() & where(is.factor))# Convert to dummy variables

# Set some naming to base the exported file names on
tune_name <- "_tuning_extra_variables"
train_name <- "_training_extra_variables"
plot_name <- "_plot_extra_variables"

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
                            loss_reduction_val = loss_reduction_val,
                            rep = repeats, seed = seed, parallel_comp = TRUE,
                            verbose = TRUE, k = k, save = TRUE, entropy_grid = F,
                            save_name = paste("XGBoost", tune_name, sep = ""))


# Train XGBoost
model_xg <- train_xgboost(df, ml_rec, save = TRUE,
                          save_name = paste("XGBoost", train_name, sep = ""),
                          mtry = select_best(tune_res_xg$tune_res,
                                           metric = metric)$mtry,
                          trees = select_best(tune_res_xg$tune_res,
                                            metric = metric)$trees,
                          min_n = select_best(tune_res_xg$tune_res,
                                            metric = metric)$min_n,
                          tree_depth = select_best(tune_res_xg$tune_res,
                                                 metric = metric)$tree_depth,
                          learn_rate = select_best(tune_res_xg$tune_res,
                                                 metric = metric)$learn_rate,
                          loss_reduction = select_best(tune_res_xg$tune_res,
                                                     metric = metric)$loss_reduction)


#=======================#
#                       #
#  Logistic Regression  #
#                       #
#=======================#

# Create preprocessing recipe for tuning and training
ml_rec_lr <- recipe(y ~ ., data = df) %>%
        step_range(all_numeric()) %>% # Min-max normalisation
        step_dummy(all_predictors() & where(is.factor)) %>% # Convert to dummy variables
        themis::step_upsample(y)

# Create preprocessing recipe for testing/predicting
test_rec_lr <- recipe(y ~ ., data = df) %>%
        step_range(all_numeric()) %>% # Min-max normalisation
        step_dummy(all_predictors() & where(is.factor))# Convert to dummy variables

# Tune hyperparamters Logistic Regression
penalty_val <- c(0, 0.001, 0.01, 0.1, 1, 5, 10, 2, 20, 50, 100, 80) # Set penalty values to evaluate
mixture_val <- c(0, 0.25, 0.5, 0.75, 1) # Set mixture values to evaluate
tune_res_lr <- tune_logistic_regression(df, ml_rec_lr, penalty_val,
                                        mixture_val, seed = seed, k = k,
                                        save = TRUE, rep = repeats,
                                        save_name = paste("LR", tune_name,
                                                          sep = ""))

# Train Logistic Regression model
model_lr <- logistic_regression(df, ml_rec_lr, save = TRUE,
                                save_name = paste("LR", train_name, sep = ""),
                                penalty = select_best(tune_res_lr$tune_res,
                                                    metric = metric)$penalty,
                                mixture = select_best(tune_res_lr$tune_res,
                                                    metric = metric)$mixture)

#=======================#
#                       #
#     Visualisation     #
#                       #
#=======================#
# Specify whether plots should be saved locally
save_plots <- "FALSE"

#XGBoost
tune_res_plot_xg <- plot_tune_res(tune_res_xg$tune_res, save = save_plots, 
                               save_name = paste("XGBoost", tune_name,
                                                 sep = ""))
#LR
tune_res_plot_lr <- plot_tune_res(tune_res_lr$tune_res, save = save_plots, 
                               save_name = paste("LR", tune_name,
                                                 sep = ""))

#XGBoost
proba_plot_xg <- plot_proba_truth(model_xg, df, test_rec, type = "violin",
                                  "XGBoost - predicted probabilities",
                                  save = save_plots,
                                  save_name =
                                          paste("XGBoost_predicted_probabilities_",
                                                plot_name, sep = ""))
# LR
proba_plot_lr <- plot_proba_truth(model_lr, df, test_rec_lr, type = "violin",
                                  "LR - predicted probabilities",
                                  save = save_plots,
                                  save_name = paste("LR_predicted_probabilities_",
                                                    plot_name, sep = ""))

# Spider/radar chart
radar_chart <- plot_radar(list("XGBoost" = model_xg,
                               "LR" = model_lr),
                          test = list("XGBoost" = df,
                                      "LR" = df),
                          dep_var = df$y, vlcex = 1.5,
                          rec = list("XGBoost" = prep(test_rec),
                                     "LR" = prep(test_rec_lr)),
                          title = "", save = save_plots,
                          legend_offset = c(0.04, -0.05),
                          save_name = paste("radar_chart", plot_name, sep = ""),
                          alpha = 0.2, chance= FALSE)

# Variable importance
shapley_info <- get_shaply_info(juice(prep(test_rec)), simplify = TRUE,
                                select_best(tune_res_xg$tune_res, "mn_log_loss"),
                                n_features = 4)
shapley_summary <- plot_shaply_summary(juice(prep(test_rec)), save = save_plots,
                                       n_features = 4,
                                       params = select_best(tune_res_xg$tune_res,
                                                   "mn_log_loss"),
                                       save_name = paste("shapley_summary",
                                                         plot_name, sep = ""))

shapley_feature_imp <- ggplot(shapley_info, aes(x = variable, y = mean_value,
                                                fill = mean_value)) +
        geom_bar(stat = "identity") +
        coord_flip() + ylab("Mean absolute variable importance") + xlab("")+
        scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
if(save_plots)
        ggsave(shapley_feature_imp,
               filename = paste("XGBoost_shapley_importances_", plot_name,
                                ".svg", sep = ""))

roc_model_lr <- pROC::roc(juice(prep(test_rec_lr))$y,
                          predict_classprob.model_fit(model_lr,
                                                      juice(prep(test_rec_lr)))[[1]],
                          ci = TRUE)
roc_model_xg <- pROC::roc(juice(prep(test_rec))$y,
               predict_classprob.model_fit(model_xg, juice(prep(test_rec)))[[1]],
               ci = TRUE)

roc_plot_xg <- ggroc(roc_model_xg, legacy.axes = TRUE) + 
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey",
                     linetype="dashed")
        
roc_plot_lr <- ggroc(roc_model_lr, legacy.axes = TRUE) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey",
                     linetype="dashed")

if(save_plots){
        ggsave(roc_plot_xg,
               filename = paste("XGBoost_ROC_curve_", plot_name,
                                ".svg", sep = ""))
        ggsave(roc_plot_lr,
               filename = paste("LR_ROC_curve_", plot_name,
                                ".svg", sep = ""))
}

print(paste("XGBoost achieved an auc of", auc(roc_model_xg), "and LR an auc of",
            auc(roc_model_lr)))


dm_xg <- construct_decision_matrix(dp, model_xg, df = df, test_rec = test_rec,
                                   model_type = "XG")
dm_lr <- construct_decision_matrix(dp, model_lr, df = df, test_rec = test_rec_lr,
                                   model_type = "LR")
dm_both <- reshape2::melt(list(LR = dm_lr[,c("misc", "ligd")],
                               XG = dm_xg[,c("misc", "ligd")],
                               NICE = dm_NICE[,c("misc", "ligd")]),
                          id.vars = "misc") %>%
        arrange(value)

# Plot saved bed days over misclassifications for both models.
ggplot(dm_both, aes(x = misc, y = value, group = L1)) +
        geom_line(aes(color = L1), size = 2, alpha = 0.5) +
        geom_point(aes(color = L1), size = 3) +
        xlim(c(0,100)) +
        scale_colour_discrete(name = "model types") +
        ylab("saved bed days") +
        xlab("false positives") +
        coord_flip()
