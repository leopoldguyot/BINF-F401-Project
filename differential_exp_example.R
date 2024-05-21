source("desq2_workflow.R")
sample_data <- read.table("data_output/clinical_data_normalized.tsv", sep = "\t", header = TRUE)
sample_data$SEX <- as.factor(sample_data$SEX)
sample_data$DTHVNT <- as.factor(sample_data$DTHVNT)
sample_data$DTHHRDY <- as.factor(sample_data$DTHHRDY)
sample_data$COHORT <- as.factor(sample_data$COHORT)
sample_data$SMPTHNTS <- as.factor(sample_data$SMPTHNTS)
row.names(sample_data) <- sample_data$SMPLID

features_count <- read.table("data/morphological_counts_lunit_dino.tsv", sep = "\t", header = TRUE, row.names = 1)
features_count <- t(as.matrix(features_count))
features_count <- features_count + 1 # Adding a pseudo count
design <- as.formula(" ~ SEX")

global_res <- diff_exp(features_count, sample_data, design)

count_depth(global_res)
resultsNames(global_res)
volcano_plot(global_res, contrast = c("SEX", "2", "1"), title = "sex volcano plot")
significant_features_count(global_res, contrast = c("SEX", "2", "1"))
