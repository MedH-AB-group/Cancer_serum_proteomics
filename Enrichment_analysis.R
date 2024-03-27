## ----libraries, message=F, cache = FALSE, include=FALSE------------------------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(knitr)
library(doParallel) 
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(forcats)
library(DT)
library(rWikiPathways)


## ----palett--------------------------------------------------------------------------------------------------------------------------------------------
##Predefined, unified palette
Palette = c("#0cc0aa", "#a20655", "#b1d34f", "#553096", "#e4b5ff", "#0b522e", "#7ebef8", "#265582", "#49d261", "#e72525", "#c2dcb8", "#834001", "#f7b8a2")


## ----read_data-----------------------------------------------------------------------------------------------------------------------------------------
###Getting data from previous analysis
metadata_prot <- get(load("results/somalogic_prot_info.RData"))
get(load("results/wilcox.RData"))
get(load("results/t-test_results.RData"))
get(load("results/GBC_benign.rfe_lasso_sig.RData"))
get(load("results/pca_results.RData"))
ptn608 <- get(load("results/probes_lasso_608ptn.RData"))
gda <- read.delim("/Users/ghada.nouairia/projects/Oscar_Jungholm_Fyfa/data/curated_gene_disease_associations.tsv.gz")
wikipathways <- readPathwayGMT("data/wikipathways/wikipathways-20240210-gmt-Homo_sapiens.gmt")


## ----venndiag------------------------------------------------------------------------------------------------------------------------------------------
###combine and VennDiagram of all analyses
ptn_lasso_comp <- ptn608 %>%
                  select(probes) %>%
                  mutate(ptn608 = "TRUE") %>%
                  full_join(sig.feat %>%
                            select(probes, RFE, Lasso),
                    by = "probes") %>%
              full_join(as.data.frame(wilcox_res) %>%
           filter(p.value <= 0.05) %>%
            rownames_to_column("probes") %>%
           select(probes) %>%
           mutate(Wilcoxon = TRUE),
           by = "probes") %>%
  full_join(Ttest_res %>%
              filter(p_value <= 0.05) %>%
              rownames_to_column("probes") %>%
              select(probes) %>%
              mutate(Ttest = TRUE),
            by = "probes") %>%
            mutate_all(~ ifelse(is.na(.), FALSE, .)) %>%
            mutate_at(c("ptn608"), as.logical)
ptn_lasso_comp <- ptn_lasso_comp %>%
                  rename("Elastic Net" = "ptn608") %>%
                  rename("U-test" = "Wilcoxon")
# pdf(file = "results/Paper_illustratios/Venndiag_1027ptn.pdf")
ggvenn(ptn_lasso_comp, c("Elastic Net", "U-test"), 
  fill_color = c("#e72525", "#b1d34f"),
  stroke_size = 0.5, set_name_size = 4)
# while (!is.null(dev.list()))  dev.off()


## ----sig_ptn-------------------------------------------------------------------------------------------------------------------------------------------
###Choose proteins by signficance (p < 0.05) and add ml results
sig_ptn <- 
  as.data.frame(wilcox_res) %>%
  filter(p.value <= 0.05) %>%
  rownames_to_column("proteins") %>%
  select(proteins) %>%
  mutate(Wilcoxon = TRUE) %>%
  full_join(Ttest_res %>%
              filter(p_value <= 0.05) %>%
              rownames_to_column("proteins") %>%
              select(proteins) %>%
              mutate(Ttest = TRUE),
            by = "proteins") %>%
  full_join(sig.feat %>%
              select(probes) %>%
              mutate(Lasso = TRUE),
            by = c("proteins" = "probes")) %>%
  full_join(ptn980 %>%
              select(probes) %>%
              mutate(ptn980 = "TRUE"),
            by = c("proteins" = "probes")) %>%
  dplyr::filter(proteins != "sample_id") %>%
  inner_join(metadata_prot %>%
               select(UniProt, EntrezGeneID, EntrezGeneSymbol),
             by = c("proteins" = "UniProt"))


## ----sig_Ab, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------------
###Setting universe and protein list
sig.ptn.list <- sig_ptn %>%
              dplyr::select(proteins)
sig.ptn <- sig.ptn.list[[1]]
ptn.Universe <- metadata_prot %>%
               select(UniProt)
ptn.Universe <- ptn.Universe[[1]]
ptn.Universe <- unlist(mget(ptn.Universe, envir=org.Hs.egSYMBOL2EG, ifnotfound = NA))


## ----go_bp_enrichment, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------
###Computing enrichment with Gene ontology, Biological pathway
ans.go <- enrichGO(gene = sig.ptn, 
                   keyType = "UNIPROT",
                   universe = ptn.Universe,
                   ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05)
datatable(as.data.frame(ans.go))

## ----go_bp_plots, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------
###Plotting in different ways and saving
barplot(ans.go, showCategory=15)
dotplot(ans.go, showCategory=15) + ggtitle("Gene Ontology: Biological Process")
upsetplot(ans.go)
ans.go1 <- pairwise_termsim(ans.go)
# pdf(file = "results/Paper_illustratios/GO_BP_network_1027ptn.pdf")
emapplot(ans.go1, showCategory = 20)
# while (!is.null(dev.list()))  dev.off()
# pdf(file = "results/Paper_illustratios/GO_BP_1027ptn.pdf")
ego <- ans.go %>% 
       mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ggplot(ego, showCategory = 15, 
       aes(richFactor, fct_reorder(Description, Count))) +
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10",
                              guide=guide_colorbar(reverse=TRUE, order=1)) +
       scale_size_continuous(range=c(2, 10)) +
       theme_dose(12) +
       xlab("Rich Factor") +
       ylab(NULL) +
       ggtitle("Gene Ontology based GSEA: Biological Process")
# while (!is.null(dev.list()))  dev.off()


## ----go_enrichment_CC, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------
###Computing enrichment with Gene ontology, Cellular components (no results expected)
ans.go <- enrichGO(gene = sig.ptn, 
                   keyType = "UNIPROT",
                   ont = "CC",
                   OrgDb ="org.Hs.eg.db",
                   universe = ptn.Universe,
                   readable=TRUE,
                   pvalueCutoff = 0.05)
datatable(as.data.frame(ans.go))


## ----ans_plots_CC, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------
##plots
barplot(ans.go, showCategory=15)
dotplot(ans.go, showCategory=15) + ggtitle("GO:CC")
upsetplot(ans.go)
ans.go1 <- pairwise_termsim(ans.go)
emapplot(ans.go1)
ego <- ans.go %>% 
       mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
# pdf(file = "results/Paper_illustratios/GO_CC_1027ptn.pdf")
ggplot(ego, showCategory = 10, 
       aes(richFactor, fct_reorder(Description, Count))) +
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10",
                              guide=guide_colorbar(reverse=TRUE, order=1)) +
       scale_size_continuous(range=c(2, 10)) +
       theme_dose(12) +
       xlab("Rich Factor") +
       ylab(NULL) +
       ggtitle("Gene Ontology GSEA in Cell Composits")
# while (!is.null(dev.list()))  dev.off()
ggplot(ego, showCategory = 15, 
       aes(Count, fct_reorder(Description, Count), fill=qvalue)) +
       geom_col() +
       scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                            guide=guide_colorbar(reverse=TRUE)) +
       theme_dose(12) +
       xlab("Normalized Enrichment Score") +
       ylab(NULL) +
       ggtitle("Gene Ontology enrichment in Cell Composits")
# cowplot::plot_grid(p3, p4, ncol=2, labels=LETTERS[1:3])


## ----go_enrichment_MF, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------
###Computing enrichment with Gene ontology, Molecular functions
ans.go <- enrichGO(gene = sig.ptn, 
                   keyType = "UNIPROT",
                   ont = "MF",
                   OrgDb ="org.Hs.eg.db",
                   universe = ptn.Universe,
                   readable = TRUE,
                   pvalueCutoff = 0.05)
datatable(as.data.frame(ans.go))


## ----ans_plots_MF, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------
###Plots and image saving as pdf
barplot(ans.go, showCategory=15)
dotplot(ans.go, showCategory=15) + ggtitle("GO based enrichment analysis (Molecular Function)")
upsetplot(ans.go)
ans.go1 <- pairwise_termsim(ans.go)
# pdf(file = "results/Paper_illustratios/GO_MF_network_1027ptn.pdf")
emapplot(ans.go1)
# while (!is.null(dev.list()))  dev.off()
ego <- ans.go %>% 
       mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
# pdf(file = "results/Paper_illustratios/GO_MF_1027ptn.pdf")
ggplot(ego, showCategory = 10, 
       aes(richFactor, fct_reorder(Description, Count))) +
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10",
                              guide=guide_colorbar(reverse=TRUE, order=1)) +
       scale_size_continuous(range=c(2, 10)) +
       theme_dose(12) +
       xlab("Rich Factor") +
       ylab(NULL) +
       ggtitle("Gene Ontology based GSEA (Molecular Function)")
# while (!is.null(dev.list()))  dev.off()
ggplot(ego, showCategory = 10, 
       aes(Count, fct_reorder(Description, Count), fill=qvalue)) +
       geom_col() +
       scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                            guide=guide_colorbar(reverse=TRUE)) +
       theme_dose(12) +
       xlab("Normalized Enrichment Score") +
       ylab(NULL) +
       ggtitle("GO based enrichment analysis (Molecular Function)")
# cowplot::plot_grid(p3, p4, ncol=2, labels=LETTERS[1:3])


## ----kegg_enrichment, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------
###Pathway enrichment using KEGG database
ans.kegg <- enrichKEGG(gene = sig.ptn,
                       organism = 'hsa',
                       keyType = "uniprot",
                       universe = ptn.Universe,
                       pvalueCutoff = 1)
datatable(as.data.frame(ans.kegg))


## ----kegg_plots, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------
###Plotting
barplot(ans.kegg, showCategory=15)
dotplot(ans.kegg, showCategory=15) + ggtitle("KEGG")
upsetplot(ans.kegg)
ans.kegg1 <- pairwise_termsim(ans.kegg)
emapplot(ans.kegg1)
ego <- ans.kegg %>% 
       mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
# pdf(file = "results/Paper_illustratios/KEGG_1027ptn.pdf")
ggplot(ego, showCategory = 15, 
       aes(richFactor, fct_reorder(Description, Count))) +
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10",
                              guide=guide_colorbar(reverse=TRUE, order=1)) +
       scale_size_continuous(range=c(2, 10)) +
       theme_dose(12) +
       xlab("Rich Factor") +
       ylab(NULL) +
       ggtitle("KEGG based pathway analysis")
# while (!is.null(dev.list()))  dev.off()
# pdf(file = "results/Paper_illustratios/KEGG_barplot_1027ptn.pdf")
ggplot(ego, showCategory = 15,   
       aes(Count, fct_reorder(Description, Count), fill=p.adjust)) +
       geom_col() +
       scale_fill_gradientn(colours=c("#f7ca64", "#b3eebe", "#7e62a3", "#371ea3"),
                            guide=guide_colorbar(reverse=TRUE)) +
       theme_dose(12) +
       xlab("Protein count") +
       ylab(NULL) +
       ggtitle("Pathway analysis based on KEGG")
# while (!is.null(dev.list()))  dev.off()


## ----sig_genes, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------
###preparing protein list for wikipathway (get gene IDs)
sig.ptn.list <- sig_ptn %>%
              dplyr::select(EntrezGeneID)
sig.ptn <- sig.ptn.list[[1]]
ptn.Universe <- metadata_prot %>%
               select(EntrezGeneID)
ptn.Universe <- ptn.Universe[[1]]
ptn.Universe <- unlist(mget(ptn.Universe, envir=org.Hs.egSYMBOL2EG, ifnotfound = NA))


## ----wikipathways, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------------------------
###Importing Wikipathways data
disease2gene <- wikipathways[, c("wpid", "gene")]
disease2name <- wikipathways[, c("wpid", "name")]
###Computing enrichment 
ans.wkp <- enricher(sig.ptn, TERM2GENE=disease2gene, TERM2NAME=disease2name, pvalueCutoff = 1, universe = ptn.Universe)
datatable(as.data.frame(ans.wkp))


## ----wkp_plots-----------------------------------------------------------------------------------------------------------------------------------------
###plotting results
barplot(ans.wkp, showCategory=15)
dotplot(ans.wkp, showCategory=15) + ggtitle("Disease db")
upsetplot(ans.wkp)
ans.wkp1 <- pairwise_termsim(ans.wkp)
emapplot(ans.wkp1)


## ----signaling, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------------
###More plotting
ego <- ans.wkp %>% 
       # filter(str_detect(Description, "Calcium")) %>%
       # filter(p.adjust < 0.001, Count > 10) %>%
       mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
# pdf(file = "results/Paper_illustratios/Wikipathwyas_1027ptn.pdf")
ggplot(ego, showCategory = 15, 
       aes(richFactor, fct_reorder(Description, Count))) +
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10",
                              guide=guide_colorbar(reverse=TRUE, order=1)) +
       scale_size_continuous(range=c(2, 10)) +
       theme_dose(12) +
       xlab("Rich Factor") +
       ylab(NULL) +
       ggtitle("Wikipathways based enrichment analysis")
# while (!is.null(dev.list()))  dev.off()

ggplot(ego, showCategory = 15, 
       aes(Count, fct_reorder(Description, Count), fill=qvalue)) +
       geom_col() +
       scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                            guide=guide_colorbar(reverse=TRUE)) +
       theme_dose(12) +
       xlab("Normalized Enrichment Score") +
       ylab(NULL) +
       ggtitle("Wikipathway based analysis")


## ----disesase_association_db, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------
###Preparing list for disease database (DisGeNet: https://www.disgenet.org/) containing > 1M gene-disease association (GDA) and > 30K diseases/disorders
disease2gene <- gda[, c("diseaseId", "geneId")]
disease2name <- gda[, c("diseaseId", "diseaseName")]
sig.ptn.list <- sig_ptn %>%
                 select(EntrezGeneID)
sig.ptn <- sig.ptn.list[[1]]
##Computing enrichment
ans.dis <- enricher(sig.ptn, TERM2GENE=disease2gene, TERM2NAME=disease2name)
datatable(as.data.frame(ans.dis))




## ----plot2, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------
###Plotting
ego <- ans.dis %>% mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>% filter(p.adjust < 0.032)
# pdf(file = "results/Paper_illustratios/DisGeNet_1027ptn.pdf")
ggplot(ego, showCategory = 15, 
       aes(richFactor, fct_reorder(Description, Count))) +
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = Count)) +
       scale_color_gradientn (colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10",
                              guide=guide_colorbar(reverse=TRUE, order=1)) +
       scale_size_continuous(range=c(2, 10)) +
       theme_dose(12) +
       xlab("Rich Factor") +
       ylab(NULL) +
       ggtitle("Disease association enrichment analysis")
# while (!is.null(dev.list()))  dev.off()
# pdf(file = "results/Paper_illustratios/DisGeNet_barplot_1027ptn.pdf")
ggplot(ego, showCategory = 15,   
       aes(Count, fct_reorder(Description, Count), fill=p.adjust)) +
       geom_col() +
       scale_fill_gradientn(colours=c("#f7ca64", "#b3eebe", "#7e62a3", "#371ea3"),
                            guide=guide_colorbar(reverse=TRUE)) +
       theme_dose(12) +
       xlab("Protein count") +
       ylab(NULL) +
       ggtitle("Disease association analysis")
# while (!is.null(dev.list()))  dev.off()

