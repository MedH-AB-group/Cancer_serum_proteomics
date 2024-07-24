# Cancer_serum_proteomics
Proteomics data analysis in the serum of gallbladder cancer patients and benign cholecystitis pre-operation


## Description of the data
Study design: Serum samples from patients who underwent tumor resection surgery for suspicion of Gallbladder Cancer (GBC) were analyzed using the 7 K panel of proteomic screening from SomaLogic. We collected the clinical data about these patients and in this document, we provide a description of the cohort.
Data: row data was delivered by SomaLogic in form of .adat file. It was preprocessed using the R package SomaDataIO (internal control of the assay removed, 3 repetitions of samples removed). 

Frequentist statistical analysis (using U-test) have been performed, lasso regression and elastic net. We also perform unsupervised and semi-supervised clustering (k-mean and hierarchical) as well as dimension reduction (Principal Component Analysis) to visualize and test the data. Finally, we perform different enrichment analysis.
