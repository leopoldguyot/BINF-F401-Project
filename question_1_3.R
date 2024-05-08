library(ggplot2)
library(dplyr)

# Replace with the path of the output data
clinical_data <- read.table("data_output/clinical_data_normalized.tsv",
    sep = "\t",
    header = TRUE
)

clinical_data$SEX <- as.factor(clinical_data$SEX)
clinical_data$DTHVNT <- as.factor(clinical_data$DTHVNT)
clinical_data$DTHHRDY <- as.factor(clinical_data$DTHHRDY)
clinical_data$COHORT <- as.factor(clinical_data$COHORT)
clinical_data$SMPTHNTS <- as.factor(clinical_data$SMPTHNTS)

clinical_data$DTHHRDY_3 <- recode(clinical_data$DTHHRDY,
    "0" = "slow", "1" = "fast", "2" = "fast", "3" = "intermediate", "4" = "slow"
)

technical_var_formula <- paste(" ~ ", paste(c("COHORT", "TRISCHD", "DTHHRDY"), collapse = " + "))

technical_var_formula_3 <- paste(" ~ ", paste(c("COHORT", "TRISCHD", "DTHHRDY_3"), collapse = " + "))

result <- lapply(c("AGE", "SEX", "HGHT", "WGHT", "BMI"), function(x) {
    print(paste(x, technical_var_formula))
    if (is.numeric(clinical_data[[x]])) {
        model <- lm(data = clinical_data, formula = paste(x, technical_var_formula))
    } else {
        model <- glm(data = clinical_data, formula = paste(x, technical_var_formula), family = "binomial")
    }
    list(variable = x, model = model)
})

result_3 <- lapply(c("AGE", "SEX", "HGHT", "WGHT", "BMI"), function(x) {
    print(paste(x, technical_var_formula))
    if (is.numeric(clinical_data[[x]])) {
        model <- lm(data = clinical_data, formula = paste(x, technical_var_formula_3))
    } else {
        model <- glm(data = clinical_data, formula = paste(x, technical_var_formula_3), family = "binomial")
    }
    list(variable = x, model = model)
})
summary(result[[1]]$model) # AGE is dependent of the cohort
summary(result[[2]]$model) # no significant association for the SEX
summary(result[[3]]$model) # no significant association for the HGHT
summary(result[[4]]$model) # no significant association for the WGHT
summary(result[[5]]$model) # no significant association for the BMI

summary(result_3[[1]]$model) # AGE is dependent of the cohort
summary(result_3[[2]]$model) # SEX is dependent of the DTHHRDY *
summary(result_3[[3]]$model) # HGHT is dependent of the DTHHRDY **
summary(result_3[[4]]$model) # WGHT is dependent of the DTHHRDY ***
summary(result_3[[5]]$model) # BMI is dependent of the DTHHRDY and cohort *

# Maybe we could add interaction between technical variables to the models

# individual modelling
result_age <- lm(data = clinical_data, formula = AGE ~ COHORT)
result_sex <- glm(data = clinical_data, formula = SEX ~ DTHHRDY_3, family = "binomial")
result_hght <- lm(data = clinical_data, formula = HGHT ~ DTHHRDY_3)
result_wght <- lm(data = clinical_data, formula = WGHT ~ DTHHRDY_3)
result_bmi <- lm(data = clinical_data, formula = BMI ~ DTHHRDY_3 + COHORT) # COHORT is almost significant
result_bmi_2 <- lm(data = clinical_data, formula = BMI ~ DTHHRDY_3) # The R square is stronlgy reduced, we will keep the first model


ggplot(data = clinical_data, aes(x = COHORT, y = AGE, color = COHORT)) +
    geom_point()

clinical_data$RES_AGE <- result_age$residuals
clinical_data$RES_SEX <- result_sex$residuals
clinical_data$RES_HGHT <- result_hght$residuals
clinical_data$RES_WGHT <- result_wght$residuals
clinical_data$RES_BMI <- result_bmi$residuals

ggplot(data = clinical_data, aes(x = seq_along(RES_AGE), y = RES_AGE)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")

write.table(clinical_data, "data_output/clinical_data_corrected.tsv", sep = "\t")
