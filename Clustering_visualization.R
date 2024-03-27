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
data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        # mutate_at("sample_id", as.character) %>%
        arrange(group) 
scaled_data <- cbind(data[1:2], as.data.frame(scale(data[3:length(data)])))
get(load("results/GBC_benign.rfe_lasso_sig.RData"))
get(load("results/wilcox.RData")) 
wilcox_res <- as.data.frame(wilcox_res)
get(load("results/t-test_results.RData"))


## ----first_iteration----------------------------------------------------------------------------------------------------------------------
###PCA plto
pca_res <- prcomp(scaled_data %>% select(-group, -sample_id), scores = TRUE)
fviz_eig(pca_res, scale = FALSE, addlabels=TRUE)

# pdf(file = "results/pca1.pdf", width = 4, height = 4) 
fviz_pca_ind(pca_res,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = Palette,
              addEllipses = TRUE,
              ellipse.type = "confidence", #confidence, t, norm, euclid, convex depending on the data type
             repel = TRUE # Avoid text overlapping
             )
# while (!is.null(dev.list()))  dev.off()
# save(pca_res, file = "results/pca_results.RData")
# save(pca_res, file = "results/pca_results_iterative.RData")
# save(pca_res, file = "results/pca_results_stattests_ML.RData")


## ----marginal_densities-------------------------------------------------------------------------------------------------------------------
###Marginal densities plots to add to PCA
components <- components %>%
  select(-sample_id) %>%
  pivot_longer(1:2)
# pdf(file = "results/Paper_illustratios/PCA_allptn_marginal_densities.pdf")
ggplot(components, aes(x = value, color = group)) +
  geom_density(alpha = 2) + 
  scale_color_manual(values = c("#0b5313", "#ec4dd8")) +
  labs(color = "Patients group") +    # Customize legend label
  theme_minimal() 
 # while (!is.null(dev.list()))  dev.off()




## ----sig_ptn------------------------------------------------------------------------------------------------------------------------------
###Defining the significant proteins form the results of previous analysis
   sig_ptn <- 
  sig.feat %>% 
          filter(Lasso == "TRUE") %>%
          filter(probes != "index") %>%    ###lasso_coef from EN analysis has this row
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
###PCA for significant proteins
sub <- PCA(sig_data %>% select(-group))
fviz_pca_ind(sub, 
                    # pointsize = "cos2", 
                    pointshape = 21, 
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
  # add_text(text=~components$group, textposition="top center", showlegend = F) %>% 
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
 # save_image(p, "results/Paper_illustratios/PCA_sel_ptn01.pdf")


## ----marginal_desnisties_sel_ptn----------------------------------------------------------------------------------------------------------
###marginal densitites to add to PCA plto
components <- components %>%
  select(-sample_id) %>%
  pivot_longer(1:2)
# pdf(file = "results/Paper_illustratios/PCA_sel_ptn01_marginal_densities.pdf")
ggplot(components, aes(x = value, color = group)) +
  geom_density(alpha = 2) + 
  scale_color_manual(values = c("#0b5313", "#ec4dd8")) +
  labs(color = "Patients group") +    # Customize legend label
  theme_minimal() 
 # while (!is.null(dev.list()))  dev.off()


## ----other_PCA_pltos----------------------------------------------------------------------------------------------------------------------
###Other PCA plots, different kinds
pca_res <- prcomp(sig_data %>% select(-group), scores = TRUE)
fviz_eig(pca_res, scale = FALSE)
# pdf(file = "results/Paper_illustratios/PCA_car_lasso_univ01_ptn.pdf")
fviz_pca_var(pca_res,
             col.var = rev("contrib"), # Color by contributions to the PC
             gradient.cols = rev(brewer.pal(3, "Set1")),   #Dark2, Set1
             repel = TRUE) +
  labs(color = "Contribution to the data variance")
 # while (!is.null(dev.list()))  dev.off()
fviz_mfa_ind(pca_res,
             habilage = "group",
             palette = Palette,
              addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
             )
# pdf(file = "results/P# pdf(file = "results/P# pdf(file = "results/Paper_illustratios/PCA_biplot_lasso_ptn.pdf")
fviz_pca_biplot(pca_res, label = "var", habillage=as.factor(sig_data$group),
                palette = c("#0b5313", "#ec4dd8"),
                col.var = "mediumblue",
               addEllipses=TRUE, ellipse.level=0.95,
               title = "PCA of the predictive proteins",
               ggtheme = theme_minimal())
# while (!is.null(dev.list()))  dev.off()


## ----kmean--------------------------------------------------------------------------------------------------------------------------------
###Kmean anlaysis using significant/chosen proteins
# pdf(file = "results/Paper_illustratios/kmean_lasso_grpcolor.pdf")
res_km_r <- eclust(sig_data %>% select(-group), "kmeans", k = 2)   ###FUNcluster = "hclust"
fviz_cluster(res_km_r, sig_data %>% select(-group), ellipse.type = "convex", palette = c("#0b5313", "#ec4dd8"), ggtheme = theme_minimal())
# while (!is.null(dev.list()))  dev.off()
# res_km <- eclust(scale(sig_data %>% select(-group)), "kmeans", k = 3)
# res_km_r <- eclust(sig_data %>% select(-group), FUNcluster = "hclust", k = 2, graph = TRUE)
# res_km_r$silinfo
# save(res_km, file = "results/kmean_results.RData")
# save(res_km, file = "results/kmean_results_stats_ML.RData")
# save(res_km, file = "results/kmean_results_allsig.RData")
# save(res_km, file = "results/kmean_iterative_PCA.RData")


## ----3d_PCA, message=F--------------------------------------------------------------------------------------------------------------------
###3D PCA plot
# sig_data <- data  %>% 
#   select(-sample_id) %>%
#   mutate_if(is.character, as.numeric)
pca_res <- prcomp(sig_data %>% select(-group), scores = TRUE)
####We removed the control, if you want it back, just comment line #2 and add : "Control", "Control", "Control" to the group array at line #4 of the code
prin_comp <- pca_res
components <- data.frame(prin_comp[["x"]]) 
# fviz_eig(prin_comp, scale = TRUE)

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
# pdf(file = "results/Paper_illustratios/heatmap_ML_ptn.pdf")
pheatmap(as.matrix(sig_data %>%
                        cor(method = "spearman") %>%
                        round(2)), cutree_rows = 3)
# while (!is.null(dev.list()))  dev.off()


## ----zoom-in_heatmap, message=FALSE, fig.cap="Protein expression is differential between groups"------------------------------------------
###Using chosen proteins, we reduce the dataset and make heatmap
sig_data <- cbind(data$group, scaled_data[,sig_ptn$proteins]) %>% 
  dplyr::rename("group" = "data$group") %>%
  mutate_if(is.character, as.numeric) %>%
  mutate_at("group", as.numeric)

# pdf(file = "results/Paper_illustratios/heatmap1_EN_univ_sel_ptn.pdf")
Heatmap(as.matrix(sig_data ), 
        row_labels = rownames(sig_data ),
        column_labels = colnames(sig_data ),
        cluster_rows = FALSE,
        row_names_max_width = unit(6, "cm"),
        cluster_columns = TRUE,
        row_names_rot = 15,
        # col = colorRamp2(c(1, 10, 20), brewer.pal(n=3, name="Set2")),  ##Paired, PiYG, PuRd is also nice ###this is when you use counts
        row_names_gp = gpar(fontsize = 3)) # Text size for row names
# while (!is.null(dev.list()))  dev.off()

