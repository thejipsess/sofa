################################################################################
### MaastrICCht cohort
### Auteur: Sander van Kuijk & Frank van Rosmalen
###
### Doel project: Predictie mortaliteit na IC-opname o.b.v. SOFA dag 1-5/7
### Doel syntax: csv data inlezen en bewerken voor analyse
###
### Start: 26/11/2021
### Laatste aanpassing: 29/11/2021
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

library(reshape2)
library(data.table)

## Verwijzing naar map waarin data staat
setwd("C:/Users/sande/Documents/Werk/sofa/data")

## Admission data uit csv-bestand inlezen
a <- read.csv("COVID19_MUMC__export_20211126.csv", header = TRUE, sep = ";")
names(a)[1] <- "Record.Id"

## Inspectie data: unieke patientnummers
sum(is.na(a$Record.Id))
length(a$Record.Id)
length(unique(a$Record.Id))

## Karakteriestieken
a$gender[a$gender == ""] <- NA
a$gender <- factor(a$gender)
sum(is.na(a$gender))

## Verwijderen die niet compleet zijn ingevoerd? i.e. zelfs geen geslacht
a <- subset(a, !is.na(a$gender))

## Daily data uit csv-bestand inlezen
d <- read.csv("DailyCRFTotal.csv", header = TRUE, sep = ",")
names(d)[1] <- "Record.Id"
length(d$Record.Id)
length(unique(d$Record.Id))

## Tijd sinds starttijd (intubatie voor mechanisch geventileerden)
## day_admission en day_intubation heel vaak 0!
d$date <- as.Date(d$date, format = "%d-%m-%Y")

setDT(d)

setorder(d, Record.Id, date)
d[, order := 1:.N, by = Record.Id]
d$Record.Id
d$order

d <- data.frame(d)
d$dag <- as.numeric(substr(d$ReportNameCustom, 10, 11))
data.frame(d$Record.Id, d$order, d$dag) ## Dubbel-check
paste(round(sum(is.na(d$dag))/length(d$Record.Id)*100), "% missing dag", sep = "")

d <- merge(a, d, by = "Record.Id")

d$PF_low <- as.numeric(d$PF_low)
d$trombocytes <- as.numeric(d$trombocytes)
d$bilirubine <- as.numeric(d$bilirubine)
d$GCS <- as.numeric(d$GCS)
d$creatinine <- as.numeric(d$creatinine)

d[d == -99] <- NA
d[d == -98] <- NA
d[d == -97] <- NA
d[d == -96] <- NA
d[d == -95] <- NA

## Berekening SOFA scores
d$SOFA_resp <- ifelse(d$PF_low < 13.3 & d$vent_mode != 8, 4,
               ifelse(d$PF_low < 26.7 & d$vent_mode != 8, 3,
               ifelse(d$PF_low < 40, 2,
               ifelse(d$PF_low < 53.3, 1,
               ifelse(d$PF_low >= 53.3, 0, NA)))))

d$SOFA_coag <- ifelse(d$trombocytes < 20, 4,
               ifelse(d$trombocytes < 50, 3,
               ifelse(d$trombocytes <100, 2,
               ifelse(d$trombocytes <150, 1,
               ifelse(d$trombocytes >= 150, 0, NA)))))

d$SOFA_live <- ifelse(d$bilirubine > 204, 4,
               ifelse(d$bilirubine >= 102, 3,
               ifelse(d$bilirubine >= 33, 2,
               ifelse(d$bilirubine >= 20, 1,
               ifelse(d$bilirubine < 20, 0, NA)))))

d$SOFA_card <- ifelse(d$vasopressors == 1 & d$nor > 0.1, 4,
               ifelse(d$vasopressors == 1 & d$nor > 0, 3,
               ifelse(d$MAP_low <70, 1,
               ifelse(d$MAP_low >= 70, 0, NA))))

d$SOFA_neur <- ifelse(is.na(d$GCS),
               (ifelse(d$sedation == 1 & d$GCS_admission <= 5, 4,
               ifelse(d$sedation == 1 & d$GCS_admission <= 9, 3,
               ifelse(d$sedation == 1 & d$GCS_admission <= 12, 2,
               ifelse(d$sedation == 1 & d$GCS_admission <= 14, 1,
               ifelse(d$sedation == 1 & d$GCS_admission == 15, 0,
               ifelse(d$sedation == 0 & d$GCS_admission <= 5, 4,
               ifelse(d$sedation == 0 & d$GCS_admission <= 9, 3,
               ifelse(d$sedation == 0 & d$GCS_admission <= 12, 2,
               ifelse(d$sedation == 0 & d$GCS_admission <= 14, 1,
               ifelse(d$sedation == 0 & d$GCS_admission == 15, 0, NA))))))))))),
              (ifelse(d$sedation == 1 & d$GCS <= 5, 4,
               ifelse(d$sedation == 1 & d$GCS <= 9, 3,
               ifelse(d$sedation == 1 & d$GCS <= 12, 2,
               ifelse(d$sedation == 1 & d$GCS <= 14, 1,
               ifelse(d$sedation == 1 & d$GCS == 15, 0,
               ifelse(d$sedation == 0 & d$GCS <= 5, 4,
               ifelse(d$sedation == 0 & d$GCS <= 9, 3,
               ifelse(d$sedation == 0 & d$GCS <= 12, 2,
               ifelse(d$sedation == 0 & d$GCS <= 14, 1,
               ifelse(d$sedation == 0 & d$GCS == 15, 0, NA))))))))))))

d$SOFA_rena <- ifelse(d$dialysis == 1, 4,
               ifelse(d$creatinine > 440 | d$urine_output < 200, 4,
               ifelse(d$creatinine >=300 | d$urine_output < 500, 3,
               ifelse(d$creatinine >= 171, 2,
               ifelse(d$creatinine >= 110, 1,
               ifelse(d$creatinine < 110, 0, NA))))))

d$SOFA_score <- rowSums(data.frame(d$SOFA_resp, d$SOFA_coag, d$SOFA_live,
                                   d$SOFA_card, d$SOFA_neur, d$SOFA_rena))

d$SOFA_score_min_neuro <-  rowSums(data.frame(d$SOFA_resp, d$SOFA_coag,
                                              d$SOFA_live, d$SOFA_card,
                                              d$SOFA_rena))

## Data opslaan om te modelleren
setwd("c:/Users/sande/Documents/Werk/sofa/data")
save(d, file = "sofa_data.Rda")

### Einde file.
