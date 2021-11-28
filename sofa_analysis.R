################################################################################
### MaastrICCht cohort
### Auteur: Sander van Kuijk
###
### Doel project: Predictie mortaliteit na IC-opname o.b.v. SOFA dag 1-5/7
### Doel syntax: model voorspelling uiteindelijk overlijden
###
### Start: 26/11/2021
### Laatste aanpassing: 28/11/2021
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

## Databestand gecreeerd met sofa_data.R inlezen
setwd("C:/Users/sande/Documents/Werk/sofa/data")
load("sofa_data.Rda")

## Persoonslevel regressie om dag 1 en 5-dagen verandering sofa score met
## minimale meetfout te bepalen
models <- lmList(SOFA_score ~ dag | Record.Id, data = d, na.action = na.omit)
coef(models)
res <- data.frame(Record.Id = rownames(coef(models)),
                  coef(models), check.names = FALSE)
names(res)[2] <- "Intercept"
names(res)[3] <- "Slope"

res$day1  <- round(res$Intercept + 1*res$Slope, 1)
res$delta <- round(4*res$Slope, 1)

## Databestand maken waarbij iedere patient één observatie bijdraagt
dp <- d[!duplicated(d$Record.Id, fromLast = TRUE), ]
dp <- merge(dp, res, by = "Record.Id")

empty <- c("X", "X.y", "X.x")
dp <- dp[, !names(dp) %in% empty]

constant <- c("CHC.Dementia", "CHC.Connective_tissue_disease", "CHC.Hemiplegia",
              "CHC.AIDS", "ckd_status.Creatinine_265_mmolL", "adrenaline",
              "dobutamine", "dopamine", "anti_viral.Isavuconazol")
dp <- dp[, !names(dp) %in% constant]

## Model om mortaliteit te voorspellen
table(dp$ICU_mortality)
dp$event <- ifelse(dp$ICU_mortality == "Death", 1, 0)
table(dp$event)

dd <- datadist(dp)
options(datadist = "dd")

## Model enkel op basis van SOFA score dag 1 en delta 1 tot 5
model <- lrm(event ~ day1 + delta, data = dp, x = TRUE, y = TRUE)
model

## Model aangevuld met geslacht en leeftijd
model2 <-  lrm(event ~ day1 + delta + gender + age + BMI, data = dp, x = TRUE, y = TRUE)
model2

## ROC curve
r <- roc(dp$event, predict(model2, type = "fitted"), ci = TRUE)
r

setwd("C:/Users/sande/Documents/Werk/sofa/figs")
png("auc.png", width = 500, height = 500, pointsize = 16)
plot(r)
dev.off()

set.seed(070181)
roc_auc <- round(0.5 + 0.5*validate(model2, B = 1000)[1, 5], 2)
roc_auc

### Einde file.
