## ----libraries, message=F, cache = FALSE, include=FALSE-----------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(knitr)
library(RColorBrewer)
library(ComplexHeatmap)
library(plotly)
library(factoextra)
library(combiroc)
library(doParallel) 
library(pROC)


## ----settings-----------------------------------------------------------------------------------------------------------------------------
###Defining a palette for visualizations
Palette = c("#49edc9", "#781486", "#0b5313", "#5cf070", "#ec4dd8", "#2e30e7")


## ----read_data, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
###Importing tabular data and scaling it
data <- read.csv2("data/protein_screening.csv")
data <- data %>%
        mutate_if(is.character, as.numeric) %>%
        mutate_at("sample_id", as.character) %>%
        arrange(group) 
scaled_data <- cbind(data[1:2], as.data.frame(scale(data[3:length(data)])))






## ----interactive_PCA, message=F, fig.cap="GBC and benign groups are very similar, in general"---------------------------------------------
### Principal component analysis with all the proteins
prin_comp <- prcomp(scaled_data %>% select(-group, -sample_id), scores = TRUE)
components <- data.frame(prin_comp[["x"]]) %>% 
  mutate (data %>%                             ###add group information for the interaction
          mutate(group = replace(group, group == "0", "Benign tumor")) %>%
          mutate(group = replace(group, group == "4", "Gallbladder cancer")) %>%
          select(group, sample_id))  
fviz_eig(prin_comp, scale = TRUE, addlabels=TRUE)  ##Eigen values
###Visualize and plot the PCA to see how data points look like according to all proteomic data
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
      title = "PC2 (6%)",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff')) 
 # save_image("results/Paper_illustratios/PCA_allptn.pdf")


## ----marginal_desnisties------------------------------------------------------------------------------------------------------------------
### Visualize marginal densities to add to PCA plots
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


## ----Mann_whitney, message=F, warning=FALSE-----------------------------------------------------------------------------------------------
###Wilcoxon test, Not paired (independant data, no one to one match like longitudinal data)
lev1 <- which(data$group %in% "4")
lev2 <- which(data$group %in% "0")
t_data <- t(scaled_data %>% select(-sample_id, -group))
eset <- t_data[,c(lev1, lev2)]

wilcox_res <- apply(eset, 1, function (x) {
                wilcox.test(x[seq_along(lev1)], x[(length(lev1) + 1):(length(lev1) + length(lev2))], paired = FALSE)$p.value })
wilcox_res <- cbind(wilcox_res,
                p.adjust(wilcox_res, method = "BH"),
                apply(eset, 1, function (x) {
                      mean(x[seq_along(lev1)]) - mean(x[(length(lev1) + 1):ncol(eset)])
                      }),
                rowMeans(eset[, seq_along(lev1)]),
                rowMeans(eset[, (length(lev1) + 1):ncol(eset)]))
rownames(wilcox_res) <- rownames(eset)
colnames(wilcox_res) <- c("p.value", "adj.p.val", "Fold.change",
                    paste("Mean", paste("GBC", collapse = "_"), sep = "_"),
                    paste("Mean", paste("Benign", collapse = "_"), sep = "_"))
# save(wilcox_res, file = "results/wilcox.RData")
###Save significant proteins with differen tthresholds
proteins <- rownames(wilcox_res)
sig_proteins <- as_data_frame(wilcox_res) 
sig_proteins$proteins <- proteins
sig_proteins5 <- sig_proteins %>%
                relocate(proteins) %>%
                filter(p.value < 0.05)
sig_proteins1 <- sig_proteins %>%
                relocate(proteins) %>%
                filter(p.value < 0.01)
sig_proteins01 <- sig_proteins %>%
                relocate(proteins) %>%
                filter(p.value < 0.001)
sig_proteins05 <- sig_proteins %>%
                relocate(proteins) %>%
                filter(p.value < 0.005)


## ----protein_subset_PCA01, fig.cap="PSC patients groups"----------------------------------------------------------------------------------
####Adding the information about the patients groups and making an interactive PCA, use only proteins with p < 0.001
sig_data <- scaled_data[,c("sample_id", "group", sig_proteins01$proteins)] %>% mutate_if(is.character, as.numeric)
prin_comp <- prcomp(sig_data %>% dplyr::select(-group))
components <- data.frame(prin_comp[["x"]]) %>% 
               mutate (data %>% 
          mutate(group = replace(group, group == "0", "B")) %>%
          mutate(group = replace(group, group == "4", "GBC")) %>%
          select(group, sample_id))  
fviz_eig(prin_comp, scale = TRUE)  ##Eigen values
plot_ly(components, x = ~PC1, y = ~PC2, color = ~components$group, 
        colors = c("#0b5313", "#ec4dd8"), type = 'scatter', mode = 'markers', showlegend = T) %>%  
  # add_text(text=~components$group, textposition="top center",
  #           showlegend = F) %>% 
  layout(
    legend = list(title = "Patients"),
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = "PCA1",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "PCA2",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))
components <- components %>%
  select(-sample_id) %>%
  pivot_longer(1:2)
# Create a density plot with faceting
ggplot(components, aes(x = value, color = group)) +
  geom_density(alpha = 0.5) + 
  scale_color_manual(values = c("#0b5313", "#ec4dd8")) +
  labs(color = "Variable") +    # Customize legend label
  theme_minimal() 


## ----roc1, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'----------------------------------------------------------------------
###ROC for P16562 (most significant protein)
roc1 <- as.data.frame(scaled_data) %>% 
  roc(group, P16562)
roc1
plot(roc1)


## ----roc2, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'----------------------------------------------------------------------
###Significant proteins, ROC
as.data.frame(scaled_data) %>% 
  roc(group, Q9NTK1)
plot(as.data.frame(scaled_data) %>% 
  roc(group, Q9NTK1))


## ----roc3, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'----------------------------------------------------------------------
as.data.frame(scaled_data) %>% 
  roc(group, O14832)
plot(as.data.frame(scaled_data) %>% 
  roc(group, O14832))
power.roc.test(roc1)


## ----roc_comb, message=FALSE, warning=FALSE, error=FALSE, fig.keep='all'------------------------------------------------------------------
###ROC with all 3 most significant proteins, combined
predictor <- rowMeans(data.frame(marker1 = scaled_data$P16562, marker2 = scaled_data$Q9NTK1, marker3 = scaled_data$O14832)) ##, marker4 = scaled_data$Q9H4D0, marker5 = scaled_data$Q93084, marker6 = scaled_data$P18440, marker7 = scaled_data$Q13561
roc(scaled_data$group, predictor)
plot(roc(scaled_data$group, predictor))


## ----protein_subset_PCA05, fig.cap="PSC patients groups"----------------------------------------------------------------------------------
####PCA with significant proteins with p < 0.01
sig_data <- scaled_data[,c("sample_id", "group", sig_proteins1$proteins)] %>% mutate_if(is.character, as.numeric)
prin_comp <- prcomp(sig_data %>% dplyr::select(-group))
components <- data.frame(prin_comp[["x"]]) %>% 
               mutate (data %>% select(group))
fviz_eig(prin_comp, scale = TRUE)
plot_ly(components, x = ~PC1, y = ~PC2, color = ~components$group, 
        colors = c("#0b5313", "#ec4dd8"), type = 'scatter', mode = 'markers', showlegend = T) %>%  
  # add_text(text=~components$group, textposition="top center",
  #           showlegend = F) %>% 
  layout(
    legend = list(title = "Patients"),
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = "PCA1",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "PCA2",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))


## ----protein_subset_PCA1, fig.cap="PSC patients groups"-----------------------------------------------------------------------------------
####Another PCA with p < 0.001
sig_data <- scaled_data[,c("sample_id", "group", sig_proteins01$proteins)] %>% mutate_if(is.character, as.numeric)
prin_comp <- prcomp(sig_data %>% dplyr::select(-group))
components <- data.frame(prin_comp[["x"]]) %>% 
               mutate (data %>% select(group))
fviz_eig(prin_comp, scale = TRUE)
plot_ly(components, x = ~PC1, y = ~PC2, color = ~components$group, 
        colors = c("#0b5313", "#ec4dd8"), type = 'scatter', mode = 'markers', showlegend = T) %>%  
  # add_text(text=~components$group, textposition="top center",
  #           showlegend = F) %>% 
  layout(
    legend = list(title = "Patients"),
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = "PCA1",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "PCA2",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))

