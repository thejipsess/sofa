################################################################################
### MaastrICCht cohort
### Auteur: Sander van Kuijk
###
### Doel project: Predictie mortaliteit na IC-opname o.b.v. SOFA dag 1-5/7
### Doel syntax: csv data inlezen en bewerken voor analyse
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
d <- read.csv("COVID19_MUMC__Daily_CRF_export_20211126.csv", header = TRUE, sep = ";")
names(d)[1] <- "Record.Id"
length(d$Record.Id)
length(unique(d$Record.Id))

## Tweede bestand met additionele patienten
d2 <- read.csv("ExportCRF_1e10dagen1.csv", header = TRUE, sep = ",")

## Beide bestanden samenvoegen
names(d) == names(d2)
d  <- d[, 1:94]
d2 <- d2[, 1:94]
d2$sedation <- ifelse(d2$sedation == 1, "Yes", "No")

d <- rbind(d, d2)

## Tijd sinds starttijd (intubatie voor mechanisch geventileerden)
## day_admission en day_intubation heel vaak 0!
d$date <- as.Date(d$date, format = "%d-%m-%Y")

setDT(d)

setorder(d, Record.Id, date)
d[, order := 1:.N, by = Record.Id]
d$Record.Id
d$order

d <- data.frame(d)
d$dag <- as.numeric(substr(d$Report.Name.Custom, 10, 11))
data.frame(d$Record.Id, d$order, d$dag) ## Dubbel-check

## Voor deze analyse: enkel dag 1 t/m 5 meenemen.
d <- subset(d, d$dag < 6)
table(d$dag)

d <- merge(a, d, by = "Record.Id")

## Berekening SOFA scores
d$SOFA_resp <- ifelse(d$PF_low < 13.3 & d$vent_mode != 'None (Not intubated)', 4,
               ifelse(d$PF_low < 26.7 & d$vent_mode != 'None (Not intubated)', 3,
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

d$SOFA_card <- ifelse(d$vasopressors == "Yes" & d$nor > 0.1, 4,
               ifelse(d$vasopressors == "Yes" & d$nor > 0, 3,
               ifelse(d$MAP_low <70, 1,
               ifelse(d$MAP_low >= 70, 0, NA))))

d$GCS <- as.numeric(d$GCS)
d$SOFA_neur <- ifelse(is.na(d$GCS),
               (ifelse(d$sedation == "Yes" & d$GCS_admission <= 5, 4,
               ifelse(d$sedation == "Yes" & d$GCS_admission <= 9, 3,
               ifelse(d$sedation == "Yes" & d$GCS_admission <= 12, 2,
               ifelse(d$sedation == "Yes" & d$GCS_admission <= 14, 1,
               ifelse(d$sedation == "Yes" & d$GCS_admission == 15, 0,
               ifelse(d$sedation == "No" & d$GCS_admission <= 5, 4,
               ifelse(d$sedation == "No" & d$GCS_admission <= 9, 3,
               ifelse(d$sedation == "No" & d$GCS_admission <= 12, 2,
               ifelse(d$sedation == "No" & d$GCS_admission <= 14, 1,
               ifelse(d$sedation == "No" & d$GCS_admission == 15, 0, NA))))))))))),
              (ifelse(d$sedation == "Yes" & d$GCS <= 5, 4,
               ifelse(d$sedation == "Yes" & d$GCS <= 9, 3,
               ifelse(d$sedation == "Yes" & d$GCS <= 12, 2,
               ifelse(d$sedation == "Yes" & d$GCS <= 14, 1,
               ifelse(d$sedation == "Yes" & d$GCS == 15, 0,
               ifelse(d$sedation == "No" & d$GCS <= 5, 4,
               ifelse(d$sedation == "No" & d$GCS <= 9, 3,
               ifelse(d$sedation == "No" & d$GCS <= 12, 2,
               ifelse(d$sedation == "No" & d$GCS <= 14, 1,
               ifelse(d$sedation == "No" & d$GCS == 15, 0, NA))))))))))))

d$SOFA_rena <- ifelse(d$dialysis == "Yes", 4,
               ifelse(d$creatinine >440 | d$urine_output < 200, 4,
               ifelse(d$creatinine >=300 | d$urine_output <500, 3,
               ifelse(d$creatinine >= 171, 2,
               ifelse(d$creatinine >= 110, 1,
               ifelse(d$creatinine < 110, 0, NA))))))

d$SOFA_score <- rowSums(data.frame(d$SOFA_resp, d$SOFA_coag, d$SOFA_live,
                                   d$SOFA_card, d$SOFA_neur, d$SOFA_rena))

## Data opslaan om te modelleren
setwd("c:/Users/sande/Documents/Werk/sofa/data")
save(d, file = "sofa_data.Rda")

### Einde file.
