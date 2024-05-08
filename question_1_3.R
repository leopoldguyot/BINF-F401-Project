# Replace with the path of the output data
clinical_data <- read.table("data/data_output.tsv",
 sep = "\t",
 header = TRUE)

clinical_data$SEX <- as.factor(clinical_data$SEX)
clinical_data$DTHVNT <- as.factor(clinical_data$DTHVNT)
clinical_data$DTHHRDY <- as.factor(clinical_data$DTHHRDY)
clinical_data$COHORT <- as.factor(clinical_data$COHORT)
clinical_data$SMPTHNTS <- as.factor(clinical_data$SMPTHNTS)

technical_var_formula <- paste(" ~ ", paste(c("COHORT","TRISCHD", "DTHHRDY"), collapse = " * "))

result <- lapply(c("AGE", "SEX", "HGHT", "WGHT", "BMI"), function(x){
    print(paste(x, technical_var_formula))
    lm(data = clinical_data, formula = paste(x, technical_var_formula))
})
summary(result[[1]]) # AGE is dependent of the cohort
summary(result[[2]]) # no significant association for the SEX
summary(result[[3]]) # no significant association for the HGHT
summary(result[[4]]) # no significant association for the WGHT
summary(result[[5]]) # no significant association for the BMI

# Overall it seems that there is an association between the COHORT and the TRISCHD

# Since only age seems to be affected by a technical variable, we model it as a function of the cohort
result_age <- lm(data = clinical_data, formula = AGE ~ COHORT)
plot(result_age)

ggplot(data = clinical_data, aes(x = COHORT, y = AGE, color = COHORT)) +
    geom_point()
### !!! AGE is not normally distributed, we should use the square transfo

clinical_data$res_age <- result_age$residuals
ggplot(data = clinical_data, aes(x = seq_along(res_age),y = res_age)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")

clinical_data$test <- clinical_data$AGE - clinical_data$res_age
