source("desq2_workflow.R")
library(tidyverse)
sample_data <- read.table("data_output/clinical_data_normalized.tsv", sep = "\t", header = TRUE)
sample_data$DTHHRDY_3 <- recode(sample_data$DTHHRDY,
                                  "0" = "slow", "1" = "fast", "2" = "fast",
                                  "3" = "intermediate", "4" = "slow" )
sample_data$SEX <- as.factor(sample_data$SEX)
sample_data$DTHVNT <- as.factor(sample_data$DTHVNT)
sample_data$DTHHRDY <- as.factor(sample_data$DTHHRDY)
sample_data$COHORT <- as.factor(sample_data$COHORT)
sample_data$SMPTHNTS <- as.factor(sample_data$SMPTHNTS)
row.names(sample_data) <- sample_data$SMPLID

features_count <- read.table("data/morphological_counts_lunit_dino.tsv", sep = "\t", header = TRUE, row.names = 1)
features_count <- t(as.matrix(features_count))
features_count <- features_count + 1 # Adding a pseudo count

non_tech_var <- c('AGE','SEX','HGHT','WGHT','BMI') 
confounding_technical_variables <- c('COHORT','DTHHRDY_3','DTHHRDY_3',
                                     'DTHHRDY_3',
                                     'DTHHRDY_3 + COHORT')

design <- as.formula(" ~ DTHHRDY_3")
vec_global_res <- c()
'''
lapply(non_tech_var,function(i){
  X <- paste("~",i)
  design <- as.formula(X)
  global_res <- diff_exp(features_count, sample_data, design)
  vec_global_res <- c(vec_global_res,global_res)
  
  return(null)
})
for(i in vec_global_res ){
  print(resultsNames(i))
}
'''

for(i in 1:5) {
  print(i)
  X <- paste("~",non_tech_var[[i]],"+",confounding_technical_variables[[i]])
  design <- as.formula(X)
  print(design)
  global_res <- diff_exp(features_count, sample_data, design)
  vec_global_res <- c(vec_global_res,global_res)
}

for(i in 1:5){
  print(resultsNames(vec_global_res[[i]]))
  print(volcano_plot(vec_global_res[[i]], 
                      name = resultsNames(vec_global_res[[i]])[2],
                      title = non_tech_var[i]))
  print(paste("significant_features_count:",
              significant_features_count(vec_global_res[[i]],
                             name = resultsNames(vec_global_res[[i]])[2])))
  
}



results(vec_global_res[[4]],name = resultsNames(vec_global_res[[4]])[2])%>% 
  as.tibble(rownames = "gene") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  View()  

