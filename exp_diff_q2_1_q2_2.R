setwd("C:/Users/Draguet/Documents/GitHub/BINF-F401-Project")

source("desq2_workflow.R")
library(tidyverse)

sample_data <- read.table("data_output/clinical_data_normalized.tsv", sep = "\t", header = TRUE)
sample_data$SEX <- as.factor(sample_data$SEX)
sample_data$DTHVNT <- as.factor(sample_data$DTHVNT)
sample_data$DTHHRDY <- as.factor(sample_data$DTHHRDY)
sample_data$COHORT <- as.factor(sample_data$COHORT)
sample_data$SMPTHNTS <- as.factor(sample_data$SMPTHNTS)
row.names(sample_data) <- sample_data$SMPLID

#reshaping the class in DTHHRDY
hardy_scale <- c()
for (i in sample_data$DTHHRDY) {
  if (i == 0 | i == 4) {
    hardy_scale <- c(hardy_scale, 0) #slow deaths
  }
  if (i == 1 | i == 2) {
    hardy_scale <- c(hardy_scale, 2) #fast death
  }
  if (i == 3) {
    hardy_scale <- c(hardy_scale, 1) #intermediate
  }
}
sample_data$DTHHRDY <- hardy_scale

#discarding samples with DTHVNT unknown
discard_id <- sample_data$SMPLID[sample_data$DTHVNT %in% c(99)]
sample_data <- sample_data[!(sample_data$DTHVNT %in% c(99)),]
sample_data$DTHHRDY <- as.factor(sample_data$DTHHRDY)

features_count <- read.table("data/morphological_counts_lunit_dino.tsv", sep = "\t", header = TRUE, row.names = 1)
features_count <- features_count[!(row.names(features_count) %in% discard_id),]
features_count <- t(as.matrix(features_count))
features_count <- features_count + 1 # Adding a pseudo count

var_names <- c(names(sample_data)[3], names(sample_data)[7:14])

asssociationsVar <- function(name){
  design <- as.formula(paste(" ~", name))
  global_res <- diff_exp(features_count, sample_data, design)
  saveRDS(global_res, paste("data_output/desq2_outputs/", name, ".rds", sep = ''))
  return(NULL)
}

lapply(var_names, asssociationsVar)



count_depth(DTHHRDY)
resultsNames(DTHHRDY)
volcano_plot(DTHHRDY, name = "DTHHRDY_1_vs_0", title = "DTHHRDY volcano plot")
volcano_plot(DTHHRDY, name = "DTHHRDY_2_vs_0", title = "DTHHRDY volcano plot")
significant_features_count(DTHHRDY, name = "DTHHRDY_1_vs_0")
significant_features_count(DTHHRDY, name = "DTHHRDY_2_vs_0")

results(DTHHRDY, contrast = c("DTHHRDY", '1', '0')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

results(DTHHRDY, contrast = c("DTHHRDY", '2', '0')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

count_depth(HGHT)
resultsNames(HGHT)
volcano_plot(HGHT, name = "HGHT", title = "Height volcano plot")
significant_features_count(HGHT, name = "HGHT")
results(HGHT, name = "HGHT") %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()