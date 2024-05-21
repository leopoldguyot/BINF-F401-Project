source("desq2_workflow.R")
library(tidyverse)
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
#levels(sample_data$DTHHRDY)
#sample_data$DTHHRDY <- relevel(sample_data$DTHHRDY, ref = "3")

global_res <- diff_exp(features_count, sample_data, design)

tibble(sample = colnames(global_res), SeqDepth = colSums(assay(global_res)))

count_depth(global_res)
resultsNames(global_res)
#results(global_res, name = "AGE")
#results(global_res, contrast = c("DTHHRDY", "2", "3"))
#volcano_plot(global_res, name = "AGE", title = "age volcano plot")
volcano_plot(global_res, name = "SEX_2_vs_1", title = "sex volcano plot")
significant_features_count(global_res, name = "SEX_2_vs_1")

results(global_res, contrast = c("SEX", "2", "1")) %>% 
    as.tibble(rownames = "gene") %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
    View()

#saveRDS(global_res, "data_output/global_res.rds")
