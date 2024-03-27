## ----libraries, message=F, cache = FALSE, include=FALSE-----------------------------------------------------------------------------------
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(knitr)
library(doParallel) 
library(pheatmap)
library(DT)
library(factoextra)
library(pROC)
library(ComplexHeatmap)
library(gtools)
library(parallel)
library(plotly)
library(caret)




## ----read_data----------------------------------------------------------------------------------------------------------------------------
###getting the data
data <- read.csv2("data/protein_screening.csv")
      data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        mutate_at(c("sample_id", "group"), as.character) %>%
        arrange(group) 

##Scaling the data
scaled_data <- cbind(data[1:2], as.data.frame(scale(data[3:length(data)])))

###Getting results from previous analyses
get(load("results/lasso_elastic_net_res.RData"))
get(load("results/u-test_res.RData")) 
wilcox_res <- as.data.frame(wilcox_res)
get(load("results/t-test_res.RData"))



## ----all_sig_ptn--------------------------------------------------------------------------------------------------------------------------
####selected proteins with (filtering statistical significance p < 0.01 and adding machine learning resutls)
sig_ptn <- 
  as.data.frame(wilcox_res) %>%
  filter(p.value <= 0.01) %>%
  rownames_to_column("proteins") %>%
  select(proteins) %>%
  mutate(Wilcoxon = TRUE) %>%
  full_join(Ttest_res %>%
              filter(p_value <= 0.01) %>%
              rownames_to_column("proteins") %>%
              select(proteins) %>%
              mutate(Ttest = TRUE),
            by = "proteins") %>%
  full_join(sig.feat %>%
              select(probes) %>%
              mutate(ML = TRUE),
            by = c("proteins" = "probes")) %>%
  dplyr::filter(proteins != "sample_id")

###Dataset with only selected proteins
sig_data <- cbind(data$group, scaled_data[,sig_ptn$proteins]) %>% 
  rename("group" = "data$group") %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_at("group", as.character)




## ----roc3, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'----------------------------------------------------------------------
####ROC curves for specific proteins
roc1 <- as.data.frame(scaled_data) %>% 
  roc(group, O14832)
plot(as.data.frame(scaled_data) %>% 
  roc(group, O14832))
power.roc.test(roc1)





## ----roc_comb, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'------------------------------------------------------------------
###ROC for combination of proteins, p < 0.001
predictor <- rowMeans(data.frame(marker1 = scaled_data$P16562, marker2 = scaled_data$Q9NTK1, marker3 = scaled_data$O14832))
roc(scaled_data$group, predictor)
plot(roc(scaled_data$group, predictor))


## ----roc1, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'----------------------------------------------------------------------
###Checking specific proteins
roc1 <- as.data.frame(scaled_data) %>% 
  roc(group, P16562)
roc1
plot(roc1)


## ----combinations-------------------------------------------------------------------------------------------------------------------------
###ROC for combinations of proteins
sig_ptn_int <- sig_ptn %>%
    filter( Wilcoxon == T & ML == T & Ttest == T)
sig_data <- cbind(data$group, scaled_data[,sig_ptn_int$proteins]) %>% 
  rename("group" = "data$group") %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_at("group", as.factor) 

# Function to calculate AUC for a subset of proteins
calculate_auc_for_subset <- function(subset_indices) {
  # Fit a logistic regression model using the subset of proteins
  model <- glm(group ~ ., data = sig_data[, c("group", subset_indices)], family = binomial(link = "logit"))
  # Make predictions
  predictions <- predict(model, newdata = sig_data[, c("group", subset_indices)], type = "response")
  # Calculate AUC
  roc_obj <- roc(sig_data$group, predictions)
  auc_score <- auc(roc_obj)
  return(auc_score)
}

# Get column names of the protein variables
protein_names <- names(sig_data)[-1]

# Initialize variables to store all combinations and AUC values
all_combinations <- data.frame(Subset = character(0), AUC = numeric(0))

# Generate all possible combinations of protein names correctly
for (subset_size in 1:length(protein_names)) {
  combinations <- combn(protein_names, subset_size)
  subset_results <- mclapply(
    1:ncol(combinations),
    function(i) {
      subset_indices <- combinations[, i]
      auc_score <- calculate_auc_for_subset(subset_indices)
      subset_info <- data.frame(
        Subset = paste(subset_indices, collapse = ", "),
        AUC = auc_score
      )
      return(subset_info)
    },
    mc.cores = 6  # Adjust the number of CPU cores to utilize
  )
  all_combinations <- rbind(all_combinations, do.call(rbind, subset_results))
}
   for (i in 1:ncol(combinations)) {
     subset_indices <- combinations[, i]
    auc_score <- calculate_auc_for_subset(subset_indices)
    subset_info <- data.frame(
          Subset = paste(subset_indices, collapse = ", "),
      AUC = auc_score
   )
    all_combinations <- rbind(all_combinations, subset_info)
  }
 }
all_combinations <- all_combinations[order(-as.numeric(all_combinations$AUC)), ]
print(all_combinations)


## ----rocs---------------------------------------------------------------------------------------------------------------------------------
###ROC for proteins in intersection of ML and stat significance
sig_data_f <- sig_data[, c("group", "Q93084", "Q9H4D0", "O14832", "Q13561", "Q9NTK1", "P16562")] %>%
             mutate_at("group", as.factor)
model <- glm(group ~ ., data = sig_data_f, family = binomial(link = "logit"))
predictor <- predict(model, newdata = sig_data_f, type = "response")
roc(sig_data_f$group, predictor)
plot(roc(sig_data_f$group, predictor))


## ----best_roc-----------------------------------------------------------------------------------------------------------------------------
###ALL 7 proteins in combinaiton in ROC curve
sig_data_f <- sig_data[, c("group", "Q93084", "Q9H4D0", "O14832", "Q13561", "P16562", "P06276.1", "P01732.1", "P04196", "A2RU54", "P18440")] %>%
             mutate_at("group", as.factor)
model <- glm(group ~ ., data = sig_data_f, family = binomial(link = "logit"))
predictor <- predict(model, newdata = sig_data_f, type = "response")
roc(sig_data_f$group, predictor)
plot(roc(sig_data_f$group, predictor))

