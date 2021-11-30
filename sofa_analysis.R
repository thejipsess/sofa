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
library(themis)
library(tictoc)

# Load local scripts
source("R/XGBoost.r")
source("R/Parallel.r")



#==============================================================================#
#
#   This Section contains the code from Sander
#
#==============================================================================#

## Databestand gecreeerd met sofa_data.R inlezen
load("Data/sofa_data.Rda")

## Voor deze analyse: subsets o.b.v. dagen
d   <- subset(d, d$dag > 0)
d10 <- subset(d, d$dag < 11)
d8  <- subset(d, d$dag < 9)
d7  <- subset(d, d$dag < 8)
d5  <- subset(d, d$dag < 6)
d3  <- subset(d, d$dag < 4)

table(d7$dag)

## Persoonslevel regressie om dag 1 en t-dagen verandering sofa score met
## minimale meetfout te bepalen

d <- d7

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
r <- roc(dp$event, predict(model3, type = "fitted"), ci = TRUE)
r

setwd("C:/Users/sande/Documents/Werk/sofa/figs")
png("auc.png", width = 500, height = 500, pointsize = 16)
plot(r)
dev.off()

set.seed(070181)
roc_auc <- round(0.5 + 0.5*validate(model2, B = 1000)[1, 5], 2)
roc_auc

## Afkappunten bepalen
dp$kans <- predict(model3, type = "fitted")
round(mean(dp$kans,na.rm = TRUE) - mean(dp$event[!is.na(dp$kans)]), 3) ## Check

range(dp$kans, na.rm = TRUE)
hist(dp$kans)

dp <- subset(dp, !is.na(dp$kans))

afkap <- seq(25, 95, 2.5)
sens  <- rep(NA, length(afkap))
spec  <- rep(NA, length(afkap))
ppv   <- rep(NA, length(afkap))
npv   <- rep(NA, length(afkap))
abst  <- rep(NA, length(afkap))
misc  <- rep(NA, length(afkap))
ligd  <- rep(NA, length(afkap))
ligp  <- rep(NA, length(afkap))

for (i in 1:length(afkap)){

 risk.FN  <- ifelse(dp$kans*100 < afkap[i], 0, 1)
 sens[i] <- round(table(risk.FN, dp$event)[2,2]/(table(risk.FN, dp$event)[2,2] +
                  table(risk.FN, dp$event)[1,2])*100, 1)
 spec[i] <- round(table(risk.FN, dp$event)[1,1]/(table(risk.FN, dp$event)[1,1] +
                  table(risk.FN, dp$event)[2,1])*100, 1)
 ppv[i]  <- round(table(risk.FN, dp$event)[2,2]/(table(risk.FN, dp$event)[2,2] +
                  table(risk.FN, dp$event)[2,1])*100, 1)
 npv[i]  <- round(table(risk.FN, dp$event)[1,1]/(table(risk.FN, dp$event)[1,1] +
                  table(risk.FN, dp$event)[1,2])*100, 1)
 abst[i] <- sum(risk.FN == 1, na.rm = TRUE)
 misc[i] <- table(risk.FN, dp$event)[2, 1]
 ligd[i] <- sum(ifelse(dp$ICU_LoS[risk.FN == 1] - 7 < 1, 0,
                       dp$ICU_LoS[risk.FN == 1] - 7), na.rm = TRUE)
 ligp[i] <- round(ligd[i]/sum(dp$ICU_LoS, na.rm = TRUE)*100, 1)

}

## Afkap in (%), testkenmerken, aantal abstineren (let op: aantal al voor 7 dagen
## van IC af, maar dat weten we bij aanvang niet), aantal onjuist abstineren,
## gespaarde ligdagen, gespaarde ligdagen (%).
data.frame(afkap, sens, spec, ppv, npv, abst, misc, ligd, ligp)

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

# Prepare date for XGBoost
df <- dp %>%
        select(event, day1, delta, age, gender) # select variables of interest

# Set dependent variable name to y
names(df)[1] <- "y"

# Factorise dependent variable
df$y <- as.factor(df$y)

# Make sure there are no missings
if(any(is.na(df))) print("WARNING!! There are still incomplete samples present")

# Create preprocessing recipe
ml_rec <- recipe(y ~ ., data = df) %>%
        #step_range(all_numeric()) %>% # Min-max normalisation
        step_dummy(all_predictors() & where(is.factor)) %>% # Convert to dummy variables
        themis::step_upsample(y)

# Set some naming to base the exported file names on
tune_name <- "_tuning"

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
                            k = 5, save = TRUE, entropy_grid = FALSE,
                            save_name = paste("XGBoost", tune_name, sep = ""))


modela <- lrm(event ~ delta, data  = dp, x = TRUE, y = TRUE)
model3 <- lrm(event ~ day1 + delta + age + gender,
              data = dp, x = TRUE, y = TRUE)
