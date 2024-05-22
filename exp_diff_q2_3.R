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



# avec mapply et apply 

asssociationsVar <- function(name,name2){
  
  design <- as.formula(paste("~",name,"+",name2))
  global_res <- diff_exp(features_count, sample_data, design)
  saveRDS(global_res, paste("data_output/desq2_outputs/",name,"_Q2_3", ".rds", sep = ''))
  return(NULL)
  
}

getfile <- function(name){
  
  file_diff_expret <- readRDS(paste0("data_output/desq2_outputs/",
                                     name,"_Q2_3",".rds"), refhook = NULL)
  
  print(resultsNames(file_diff_expret[0]))
  print(volcano_plot(file_diff_expret, 
                     name = resultsNames(file_diff_expret[0])[2],
                     title =paste(resultsNames(file_diff_expret[0])[2],
                                  "with adjustment for the confounding technical variables") ))
  print(paste("significant_features_count:",
              significant_features_count(file_diff_expret,
                                         name = resultsNames(file_diff_expret[0])[2])))
  results(file_diff_expret,name = resultsNames(file_diff_expret[0])[2])%>% 
    as.tibble(rownames = "gene") %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
    View(title = paste("Significant Results for", name)) 
  
  return(NULL)
}


mapply(asssociationsVar,non_tech_var,confounding_technical_variables)


lapply(non_tech_var, getfile)







