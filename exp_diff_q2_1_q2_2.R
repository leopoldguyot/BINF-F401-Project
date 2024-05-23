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


resultsNames(AGE)
results(HGHT, name = "AGE") %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(SEX)
results(SEX, contrast = c("SEX", '2', '1')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(HGHT)
results(HGHT, name = "HGHT") %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(WGHT)
results(WGHT, name = "WGHT") %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(BMI)
results(BMI, name = "BMI") %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(COHORT)
results(COHORT, contrast = c("COHORT", 'Postmortem', 'Organ Donor (OPO)')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(TRISCHD)
results(TRISCHD, name = "TRISCHD") %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

resultsNames(DTHVNT)
results(DTHVNT, contrast = c("DTHVNT", '1', '0')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()


resultsNames(DTHHRDY)
results(DTHHRDY, contrast = c("DTHHRDY", '1', '0')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

results(DTHHRDY, contrast = c("DTHHRDY", '2', '0')) %>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()

