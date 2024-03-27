all: set_up, Cohort_description, data_preprocess, Statistics, ML, diag_test, Enrichment

####get clinical data and raw data and place them in a folder called data/

set_up:
	mkdir results/

Cohort_description:
	scripts/cohort_description.R

data_preprocess:
	scripts/Data_preprocessing.R


Statistics:
	scripts/Mann_Whitney_test.R
	scripts/Statistical_analysis.R

ML:
	scripts/Machine_learning.R

diag_test:
	scripts/ROC_combinations.R
	scripts/Clustering_visualization.R

Enrichment:
	scripts/Enrichment_analysis.R

