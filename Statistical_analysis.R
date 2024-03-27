## ----libraries, message=F, cache = FALSE, include=FALSE-----------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(knitr)
library(doParallel) 
library(limma)
library(plotly)
library(QFeatures)




## ----read_data----------------------------------------------------------------------------------------------------------------------------
###Get data in tabular form
data <- read.csv2("data/protein_screening.csv")
data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        mutate_at(c("sample_id", "group"), as.character) %>%
        arrange(group) 

###results from mann whitney test (see Mann_Whitney_test.R)
get(load("results/u-test-res.RData"))


## ----intensities--------------------------------------------------------------------------------------------------------------------------
###Intensities of the signifciant proteins with 0.005
sig_prot_05 <- as.data.frame(wilcox_res) %>%
               rownames_to_column("proteins") %>%
                filter(p.value < 0.005)
##log2 transformation necessary for limma density plot
data_log2 <- cbind(data[1:2], as.data.frame(log2(data[3:length(data)])))

limma::plotDensities(data_log2[,sig_prot_05$proteins])



## ----protein_distrib----------------------------------------------------------------------------------------------------------------------
####Boxplots of the significant proteins with 0.005

boxplot(data_log2[,sig_prot_05$proteins],
    col = palette()[-1],
    main = "Protein distribtutions", ylab = "Protein concentration (log2)"
)



## ----MDS----------------------------------------------------------------------------------------------------------------------------------
###MDS of proteins with 0.005
sig_data <- cbind(data_log2[,1:2], data_log2[,sig_prot_05$proteins])
limma::plotMDS(data_log2[,sig_prot_05$proteins])


## ----cal_pval-----------------------------------------------------------------------------------------------------------------------------
###Calculating t tests for all proteins
t_data <- t(data %>% 
              select(-group, -sample_id))
t.test(t_data[1,1:38], t_data[1,39:82])
ttest <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value}

##add p values
p_value = apply(t_data, 1, ttest, grp1 = c(1:38), grp2 = c(39:82))

##VIsualize histograms
hist(p_value)


## -----------------------------------------------------------------------------------------------------------------------------------------
###Add adjusted p values and visualize histograms
adj_p_value = p.adjust(p_value, method = "bonferroni")
hist(adj_p_value)


## ----cal_FC-------------------------------------------------------------------------------------------------------------------------------
###add fold change and visualize histograms
t_data = log2(t_data)
ctrl = apply(t_data[,1:38], 1, mean)
cancer = apply(t_data[,39:82], 1, mean)
FC <- ctrl - cancer
hist(FC, xlab = "log2 Fold Change of benign vs GBC")


## ----comb_results-------------------------------------------------------------------------------------------------------------------------
###Combine all calculated parameters in one table and save it
Ttest_res = as.data.frame(cbind(FC, p_value, adj_p_value))
save(Ttest_res, file = "results/t-test_res.RData")


## ----volcano------------------------------------------------------------------------------------------------------------------------------
###Volcano plots for T-test and U-test

ggplot(Ttest_res,
    aes(x = FC, y = -log10(p_value), color = p_value < 0.05)) +
    xlab("Fold change") + ylab ("-log10 (p value)") +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("royalblue4", "hotpink"), 0.5)) +
    theme_minimal() +
    ggtitle("Protein frequency difference between benign and GBC (U-test)")

ggplot(as.data.frame(wilcox_res),
    aes(x = Fold.change, y = -log10(p.value), color = p.value < 0.05)) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("royalblue4", "hotpink"), 0.5)) +
    theme_minimal() +
    ggtitle("Mann Whitney test based volcano plot")

