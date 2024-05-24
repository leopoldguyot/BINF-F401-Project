library(dplyr)
clinical_data_clean <- read.table("data/clinical_data.tsv", sep = "\t", header = TRUE)
clinical_data_normalized <- read.table("data_output/clinical_data_normalized.tsv", sep = "\t", header = TRUE)

AGE_normalized <- clinical_data_normalized$AGE
HGHT_normalized <- clinical_data_normalized$HGHT
WGHT_normalized <- clinical_data_normalized$WGHT
BMI_normalized <- clinical_data_normalized$BMI
TRISCHD_normalized <- clinical_data_normalized$TRISCHD

normalized <- data.frame(AGE_normalized = AGE_normalized, HGHT_normalized = HGHT_normalized, WGHT_normalized = WGHT_normalized, BMI_normalized = BMI_normalized, TRISCHD_normalized = TRISCHD_normalized)

complete_df <- cbind(clinical_data_clean, normalized)
complete_df$DTHHRDY_3 <- recode(complete_df$DTHHRDY,
    "0" = "slow", "1" = "fast", "2" = "fast", "3" = "intermediate", "4" = "slow"
)

breaks_age <- c(-Inf, 55, Inf)
breaks_hght <- c(-Inf, 69, Inf)
breaks_wght <- c(-Inf, 180, Inf)
breaks_bmi <- c(-Inf, 27, Inf)
breaks_trischd <- c(-Inf, 500, Inf)

# Create new categorical variables
complete_df$AGE_cat <- cut(complete_df$AGE, breaks = breaks_age, labels = c("0-55", "56+"))
complete_df$HGHT_cat <- cut(complete_df$HGHT, breaks = breaks_hght, labels = c("0-70", "71+"))
complete_df$WGHT_cat <- cut(complete_df$WGHT, breaks = breaks_wght, labels = c("0-180", "181+"))
complete_df$BMI_cat <- cut(complete_df$BMI, breaks = breaks_bmi, labels = c("0-27", "28+"))
complete_df$TRISCHD_cat <- cut(complete_df$TRISCHD, breaks = breaks_trischd, labels = c("0-500", "501+"))

# Save the new data
write.table(complete_df, "data_output/clinical_data_preprocessed.tsv", sep = "\t", row.names = FALSE)
