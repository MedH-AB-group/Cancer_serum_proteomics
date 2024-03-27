## ----libraries, message=F, cache = FALSE, include=FALSE-----------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(knitr)
library(SomaDataIO)
library(doParallel)



## ----read_data----------------------------------------------------------------------------------------------------------------------------
###Read raw data
raw_data <- read_adat("data/KAR-2342070_2023-07-10/SS-2342070_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat")
is.soma_adat(raw_data)


## ----summary------------------------------------------------------------------------------------------------------------------------------
###Visualize, summarize and control data, calculate min, median, 1Q, 3Q, SD etc. for each protein
summary(raw_data[,34:50])
raw_data[, 34:50] %>%
  split(raw_data$SlideId) %>%
  lapply(summary)


## ----save_prot_metadata-------------------------------------------------------------------------------------------------------------------
### get the metadata for the measured proteins and save them to RData file
attr(raw_data, "Col.Meta")
getAnalyteInfo(raw_data) %>%
  save(file = "results/somalogic_prot_info.RData")


## ----save_data----------------------------------------------------------------------------------------------------------------------------
#### Save the data in tabular form (RData) for downstream usage
data <- raw_data %>% 
  filter(SampleType == "Sample")   ##Remove calibrators
data <- data[, 34:length(data)]
save(data, file = "results/pre_processed_proteomics_data.RData")
###Transform all data to log(10)
data <- data %>% log(10)
save(data, file = "results/pre_processed_log_proteomics_data.RData")

