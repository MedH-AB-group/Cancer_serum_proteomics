## ----libraries, message=F, cache = FALSE, include=FALSE-----------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(gplots)
library(dplyr)
library(knitr)
library(RColorBrewer)
library(stringr)
library("StepReg")
library(ComplexHeatmap)
library(ggpubr)
library(plotly)
library(doParallel) 
library(factoextra)
library(FactoMineR)
library(irlba)


## ----settings-----------------------------------------------------------------------------------------------------------------------------
###Palette of colors pre-defined and unified for all plots
Palette = c("#0cc0aa", "#a20655", "#b1d34f", "#553096", "#e4b5ff", "#0b522e", "#7ebef8", "#265582", "#49d261", "#e72525", "#c2dcb8", "#834001", "#f7b8a2")


## ----read_data----------------------------------------------------------------------------------------------------------------------------
###Getting tabular data and results from previous statistic tests and ML
data <- read.csv2("data/protein_screening.csv")
###Data cleaning
data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        # mutate_at("sample_id", as.character) %>%
        arrange(group) 
###Scaling the data for clustering and dimension reduction
scaled_data <- cbind(data[1:2], as.data.frame(scale(data[3:length(data)])))
###Importing results from machine learning, statistical tests to get the names of proteins to use in the anlysis
get(load("results/lasso_elastic_net_res.RData"))  
get(load("results/u-test_res.RData")) 
wilcox_res <- as.data.frame(wilcox_res)
get(load("results/t-test_res.RData"))

## ----first_iteration----------------------------------------------------------------------------------------------------------------------
###PCA analysis and plot to see how data looks like in all dimensions (7K proteins)
pca_res <- prcomp(scaled_data %>% select(-group, -sample_id), scores = TRUE)
fviz_eig(pca_res, scale = FALSE, addlabels=TRUE)  ###Eigenvalues
fviz_pca_ind(pca_res,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = Palette,
              addEllipses = TRUE,
              ellipse.type = "confidence", #confidence, t, norm, euclid, convex depending on the data type
             repel = TRUE # Avoid text overlapping)

## ----marginal_densities-------------------------------------------------------------------------------------------------------------------
###Marginal densities plots to add to PCA
components <- components %>%
  select(-sample_id) %>%
  pivot_longer(1:2)
ggplot(components, aes(x = value, color = group)) +
  geom_density(alpha = 2) + 
  scale_color_manual(values = c("#0b5313", "#ec4dd8")) +
  labs(color = "Patients group") +    # Customize legend label
  theme_minimal() 


## ----sig_ptn------------------------------------------------------------------------------------------------------------------------------
###Defining the significant proteins form the results of previous analysis
   sig_ptn <- 
  sig.feat %>% 
          filter(Lasso == "TRUE") %>%
          filter(probes != "index") %>%    ###lasso_coef from EN analysis have this row
          select(probes) %>%
          rename("proteins" = "probes") %>%
           mutate(Lasso = TRUE) %>%
  full_join(as.data.frame(wilcox_res) %>%
           filter(p.value <= 0.01) %>%
            rownames_to_column("proteins") %>%
           select(proteins) %>%
           mutate(Wilcoxon = TRUE),
           by = "proteins") %>%
  full_join(Ttest_res %>%
              filter(p_value <= 0.01) %>%
              rownames_to_column("proteins") %>%
              select(proteins) %>%
              mutate(Ttest = TRUE),
            by = "proteins") %>%
  dplyr::filter(proteins != "sample_id")
sig_data <- cbind(data$group, scaled_data[,sig_ptn$proteins]) %>% 
  dplyr::rename("group" = "data$group") %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_at("group", as.factor)


## ----pca_sel_ptn--------------------------------------------------------------------------------------------------------------------------
###PCA analysis and plot to see how data looks like in the selected proteins dimensions
sub <- PCA(sig_data %>% select(-group))
fviz_pca_ind(sub, pointshape = 21, 
                    fill = as.factor(sig_data$group),
                    repel = TRUE, # Avoid text overlapping (slow if many points)
                    geom = c("text","point"), 
                    addEllipses = TRUE,
                    xlab = "PC1", ylab = "PC2",label = "sample_id")


## ----interactive_PCA, fig.cap="PSC patients groups"---------------------------------------------------------------------------------------
###Another PCA with nicer plot and eigenvalues
prin_comp <- prcomp(sig_data %>% select(-group), scores = TRUE)
components <- data.frame(prin_comp[["x"]]) %>% 
  mutate (data %>%                             ###add group information for the interaction
          mutate(group = replace(group, group == "0", "Benign tumor")) %>%
          mutate(group = replace(group, group == "4", "Gallbladder cancer")) %>%
          select(group, sample_id))  
fviz_eig(prin_comp, scale = TRUE, addlabels=TRUE)  ##Eigen values

plot_ly(components, x = ~PC1, y = ~PC2, color = ~components$group, 
        colors = c("#0b5313", "#ec4dd8"), type = 'scatter', mode = 'markers', showlegend = T) %>%  
  layout(
    legend = list(title = "Patients"),
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = "PC1 (18%)",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "PC2 (7%)",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))
             

## ----marginal_desnisties_sel_ptn----------------------------------------------------------------------------------------------------------
###marginal densitites to add to PCA plto
components <- components %>%
  select(-sample_id) %>%
  pivot_longer(1:2)

ggplot(components, aes(x = value, color = group)) +
  geom_density(alpha = 2) + 
  scale_color_manual(values = c("#0b5313", "#ec4dd8")) +
  labs(color = "Patients group") +    # Customize legend label
  theme_minimal() 


## ----other_PCA_pltos----------------------------------------------------------------------------------------------------------------------
###Other PCA plots, different kinds
pca_res <- prcomp(sig_data %>% select(-group), scores = TRUE)
fviz_eig(pca_res, scale = FALSE)

fviz_pca_var(pca_res,
             col.var = rev("contrib"), # Color by contributions to the PC
             gradient.cols = rev(brewer.pal(3, "Set1")),   #Dark2, Set1
             repel = TRUE) +
  labs(color = "Contribution to the data variance")

fviz_mfa_ind(pca_res,
             habilage = "group",
             palette = Palette,
              addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
             )

fviz_pca_biplot(pca_res, label = "var", habillage=as.factor(sig_data$group),
                palette = c("#0b5313", "#ec4dd8"),
                col.var = "mediumblue",
               addEllipses=TRUE, ellipse.level=0.95,
               title = "PCA of the predictive proteins",
               ggtheme = theme_minimal())



## ----kmean clustering--------------------------------------------------------------------------------------------------------------------------------
###Kmean anlaysis using significant/chosen proteins

res_km_r <- eclust(sig_data %>% select(-group), "kmeans", k = 2)   ###use FUNcluster = "hclust" to get hierarchical clustering
fviz_cluster(res_km_r, sig_data %>% select(-group), ellipse.type = "convex", palette = c("#0b5313", "#ec4dd8"), ggtheme = theme_minimal())
             
# res_km <- eclust(scale(sig_data %>% select(-group)), "kmeans", k = 3)  ###The best K is actually 3 not 2, where the third cluster is in between the bening and GBC groups and has data points that are very difficult to differentiate
# res_km_r <- eclust(sig_data %>% select(-group), FUNcluster = "hclust", k = 2, graph = TRUE)



## ----3d_PCA, message=F--------------------------------------------------------------------------------------------------------------------
###3D PCA plot
pca_res <- prcomp(sig_data %>% select(-group), scores = TRUE)
prin_comp <- pca_res
components <- data.frame(prin_comp[["x"]]) 
fviz_eig(prin_comp, scale = TRUE)

tot_explained_variance_ratio <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)

plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, 
        color = ~sig_data$group, colors = brewer.pal(8,"Dark2") ) %>%
  add_markers(size = 12) %>%
  layout(
    title = 'Total Explained Variance = 99.48',
    scene = list(bgcolor = "#e5ecf6"))
head(sort(pca_res$rotation[,1], decreasing=TRUE))


## ----correlation_heatmap------------------------------------------------------------------------------------------------------------------
####determine clusters in raws and columns using spearman, calculated as distance from each-other and in proportion
sig_data <- data[,sig_ptn$proteins] %>% mutate_if(is.character, as.numeric)

pheatmap(as.matrix(sig_data %>%
                        cor(method = "spearman") %>%
                        round(2)), cutree_rows = 3)



## ----heatmap, message=FALSE, fig.cap="Protein expression is differential between groups"------------------------------------------
###Using chosen proteins, we reduce the dataset and make heatmap
sig_data <- cbind(data$group, scaled_data[,sig_ptn$proteins]) %>% 
  dplyr::rename("group" = "data$group") %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_at("group", as.numeric)


Heatmap(as.matrix(sig_data ), 
        row_labels = rownames(sig_data ),
        column_labels = colnames(sig_data ),
        cluster_rows = FALSE,
        row_names_max_width = unit(6, "cm"),
        cluster_columns = TRUE,
        row_names_rot = 15,
        # col = colorRamp2(c(1, 10, 20), brewer.pal(n=3, name="Set2")),  ##Paired, PiYG, PuRd is also nice ###this is when you use counts
        row_names_gp = gpar(fontsize = 3)) # Text size for row names


