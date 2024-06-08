library(ggplot2)
library(corrplot)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)

clinical <- read.table(file = "clinical_data.tsv", sep = "\t", header = TRUE)


# Factorial Analysis of Mixed Data

val_matrix <- cbind(clinical[3:8], clinical[10:12])

sex_genre <- c()
for (i in clinical$SEX) {
    if (i == 1) {
        sex_genre <- c(sex_genre, "Male")
    }
    if (i == 2) {
        sex_genre <- c(sex_genre, "Female")
    }
}
val_matrix$SEX <- sex_genre

vent <- c()
for (i in clinical$DTHVNT) {
    if (i == 0) {
        vent <- c(vent, "No")
    }
    if (i == 1) {
        vent <- c(vent, "Yes")
    }
    if (i == 98) {
        vent <- c(vent, "Not Reported")
    }
    if (i == 99) {
        vent <- c(vent, "Unknown")
    }
}
val_matrix$DTHVNT <- vent

hardy_scale <- c()
for (i in clinical$DTHHRDY) {
    if (i == 0) {
        hardy_scale <- c(hardy_scale, "Ventilator Case")
    }
    if (i == 1) {
        hardy_scale <- c(hardy_scale, "Violent and Fast Death")
    }
    if (i == 2) {
        hardy_scale <- c(hardy_scale, "Fast Death of Natural Causes")
    }
    if (i == 3) {
        hardy_scale <- c(hardy_scale, "Intermediate Death")
    }
    if (i == 4) {
        hardy_scale <- c(hardy_scale, "Slow Death")
    }
}
val_matrix$DTHHRDY <- hardy_scale

clinical_famd <- FAMD(val_matrix, ncp = 10, graph = TRUE)


# Correlation between quantitative variables

val_matrix <- cbind(clinical[3:8], clinical[10:12])
val_matrix$COHORT <- as.numeric(factor(val_matrix$COHORT))

cor.test(val_matrix$AGE, val_matrix$HGHT, method = "spearman") # not significant
cor.test(val_matrix$AGE, val_matrix$WGHT, method = "spearman") # not significant
cor.test(val_matrix$AGE, val_matrix$BMI, method = "spearman") # not significant
cor.test(val_matrix$AGE, val_matrix$TRISCHD, method = "spearman") # significant

cor.test(val_matrix$HGHT, val_matrix$WGHT, method = "spearman") # significant
cor.test(val_matrix$HGHT, val_matrix$BMI, method = "spearman") # significant
cor.test(val_matrix$HGHT, val_matrix$TRISCHD, method = "spearman") # significant

cor.test(val_matrix$WGHT, val_matrix$BMI, method = "spearman") # significant
cor.test(val_matrix$WGHT, val_matrix$TRISCHD, method = "spearman") # significant

cor.test(val_matrix$BMI, val_matrix$TRISCHD, method = "spearman") # not significant

p_val_M_corr <- c(
  0.6228, 0.7948, 0.9113, 0.001322, 2.2e-16, 0.001282, 2.693e-05,
  2.2e-16, 0.00193, 0.5186
)
p.adjust(p_val_M_corr, method = "bonferroni") # stay the same

# correlation matrix with mixed-data (or maybe need to make ANOVA, chi-square and Pearson's correlation test)

test <- model.matrix(~ 0 + ., data = val_matrix)
test_corr <- cor(test, use = "pairwise.complete.obs")
ggcorrplot(test_corr, show.diag = FALSE, type = "lower", lab = TRUE, lab_size = 2)


# Chi-square between categories
val_matrix_categ <- cbind(clinical[3:4], clinical[11:12])

test_data <- table(val_matrix$COHORT, val_matrix$SEX)
print(test_data)
test_chisq <- chisq.test(table(val_matrix$COHORT, val_matrix$SEX))
test_chisq$p.value

col_names <- colnames(val_matrix_categ)
chisq_pval <- c()
chisq_names <- c()
for (i in 1:3) {
    for (j in i + 1:4) {
        if (j < 5) {
            chisq_calc <- chisq.test(table(val_matrix_categ[[i]], val_matrix_categ[[j]]), simulate.p.value = TRUE)
            chisq_pval <- c(chisq_pval, chisq_calc$p.value)
            chisq_names <- c(chisq_names, paste(col_names[i], "x", col_names[j]))
        }
    }
}
corrected_pval <- p.adjust(chisq_pval, method = "bonferroni")
results_categ <- data.frame(
    Comparaison = chisq_names,
    Pval = corrected_pval
)
results_categ # devrait pas avoir correlation entre sex et les autres ?


# Correlation between categorical variables => Cramer's V
library(rcompanion)

chisq.test(table(val_matrix$SEX, val_matrix$COHORT)) # p-val ptet à corriger, ici sont bien dépendantes
cramerV(table(val_matrix$SEX, val_matrix$COHORT), bias.correct = TRUE) # correlation value, allowed because base on X2, non-param
# 0.1999
# 0 - 0.2 -> weak association
# 0.2 - 0.6 -> middle association
# 0.6 - 1 -> strong association

chisq.test(table(val_matrix$SEX, val_matrix$DTHVNT)) # may be incorrect because some cells are very low => only consider No and Yes
chisq.test(
    val_matrix$SEX[val_matrix$DTHVNT != "Unknown"],
    val_matrix$DTHVNT[val_matrix$DTHVNT != "Unknown"]
)
DTHVNT_remov <- as.data.frame.matrix(table(val_matrix$SEX, val_matrix$DTHVNT))
DTHVNT_remov[2] <- NULL
cramerV(as.matrix(DTHVNT_remov), bias.correct = TRUE) # 0.1868

chisq.test(table(val_matrix$SEX, val_matrix$DTHHRDY)) # may delete last column
chisq.test(
    val_matrix$SEX[val_matrix$DTHHRDY != "Violent and Fast Death"],
    val_matrix$DTHHRDY[val_matrix$DTHHRDY != "Violent and Fast Death"]
)
DTHHRDY_table <- as.data.frame.matrix(table(val_matrix$SEX, val_matrix$DTHHRDY))
DTHHRDY_table[5] <- NULL
cramerV(as.matrix(DTHHRDY_table), bias.correct = TRUE) # 0.2255

chisq.test(table(val_matrix$COHORT, val_matrix$DTHVNT)) # may be incorrect
chisq.test(
    val_matrix$COHORT[val_matrix$DTHVNT != "Unknown"],
    val_matrix$DTHVNT[val_matrix$DTHVNT != "Unknown"]
)
cramerV(table(val_matrix$COHORT, val_matrix$DTHVNT), bias.correct = TRUE) # 0.8613

chisq.test(table(val_matrix$COHORT, val_matrix$DTHHRDY)) # may be incorrect
chisq.test(
    val_matrix$COHORT[val_matrix$DTHHRDY != "Violent and Fast Death"],
    val_matrix$DTHHRDY[val_matrix$DTHHRDY != "Violent and Fast Death"]
)
DTHHRDY_table2 <- as.data.frame.matrix(table(val_matrix$COHORT, val_matrix$DTHHRDY))
DTHHRDY_table2[5] <- NULL
cramerV(as.matrix(DTHHRDY_table2), bias.correct = TRUE) # 0.8813

chisq.test(table(val_matrix$DTHVNT, val_matrix$DTHHRDY)) # may be incorrect
n_DTHVNT <- val_matrix$DTHVNT[val_matrix$DTHHRDY != "Violent and Fast Death"]
n_DTHHRDY <- val_matrix$DTHHRDY[val_matrix$DTHHRDY != "Violent and Fast Death"]
chisq.test(
    n_DTHVNT[n_DTHVNT != "Unknown"],
    n_DTHHRDY[n_DTHVNT != "Unknown"]
)
DTHHRDY_table <- as.data.frame.matrix(table(val_matrix$COHORT, val_matrix$DTHHRDY))
DTHHRDY_table[5] <- NULL
cramerV(table(n_DTHVNT, n_DTHHRDY), bias.correct = TRUE) # 0.6855

chisq_pval <- c(0.0006642102, 0.001469577, 0.0006295674, 1.727819e-47, 1.589829e-47, 6.091008e-57)
corr_chisq <- p.adjust(chisq_pval, method = "bonferroni") # still below 0.05
cramerCorr <- data.frame(
    SEX = c(1, 0.1999, 0.1868, 0.2255),
    COHORT = c(0.1999, 1, 0.8613, 0.8813),
    DTHVNT = c(0.1868, 0.8613, 1, 0.6855),
    DTHHRDY = c(0.2255, 0.8813, 0.6855, 1)
)
rownames(cramerCorr) <- c("SEX", "COHORT", "DTHVNT", "DTHHRDY")
corrplot(as.matrix(cramerCorr))


# logistic regression with binomial variables and quantitatives variables
summary(glm(as.factor(SEX) ~ AGE, family = binomial, data = val_matrix)) # 0.02397 slope not significant
summary(glm(as.factor(SEX) ~ HGHT, family = binomial, data = val_matrix)) # 0.8233
summary(glm(as.factor(SEX) ~ WGHT, family = binomial, data = val_matrix)) # 0.037508
summary(glm(as.factor(SEX) ~ BMI, family = binomial, data = val_matrix)) # 0.06833 #not significant
summary(glm(as.factor(SEX) ~ TRISCHD, family = binomial, data = val_matrix)) # 0.0012615
summary(glm(as.factor(COHORT) ~ AGE, family = binomial, data = val_matrix)) # 0.06545 slope
summary(glm(as.factor(COHORT) ~ HGHT, family = binomial, data = val_matrix)) # 0.10607
summary(glm(as.factor(COHORT) ~ WGHT, family = binomial, data = val_matrix)) # 0.006206 #not significant
summary(glm(as.factor(COHORT) ~ BMI, family = binomial, data = val_matrix)) # 0.002597 not significant
summary(glm(as.factor(COHORT) ~ TRISCHD, family = binomial, data = val_matrix)) # 0.009701
test <- c(0.0448, 6.15e-16, 1.7e-13, 0.0286, 0.00024, 6.82e-07, 0.000683, 0.0522, 0.928, 2e-16)
test <- p.adjust(test, method = "bonferroni")
test

log_reg <- data.frame(
    SEX = c(0, 0.8233, 0.037508, 0, 0.0012615),
    COHORT = c(0.06545, 0.10607, 0, 0, 0.009701)
)
row.names(log_reg) <- c("AGE", "HGHT", "WGHT", "BMI", "TRSCHD")
corrplot(as.matrix(log_reg))

# Intraclass correlation coefficient (between categ and quanti)
library(irr)

kruskal.test(AGE ~ DTHVNT, data = val_matrix) # 0.0003567, not the same distrib per group
# using corr coef calc (sqrt(X2 / n_obs)) coef = 0.4815845
sqrt(15.877 / 288) # 0.2347945

kruskal.test(HGHT ~ DTHVNT, data = val_matrix) # 0.001951
sqrt(12.479 / 288) # 0.2081583

kruskal.test(WGHT ~ DTHVNT, data = val_matrix) # 0.07141

kruskal.test(BMI ~ DTHVNT, data = val_matrix) # 0.9119

kruskal.test(TRISCHD ~ DTHVNT, data = val_matrix) # 2.2e-16
sqrt(173.33 / 288) # 0.7757837

kruskal.test(AGE ~ DTHHRDY, data = val_matrix) # 2.372e-05
sqrt(26.62 / 288) # 0.3040239

kruskal.test(HGHT ~ DTHHRDY, data = val_matrix) # 9.619e-05
sqrt(23.597 / 288) # 0.2862412

kruskal.test(WGHT ~ DTHHRDY, data = val_matrix) # 0.0004999
sqrt(19.998 / 288) # 0.26351

kruskal.test(BMI ~ DTHHRDY, data = val_matrix) # 0.1106

kruskal.test(TRISCHD ~ DTHHRDY, data = val_matrix) # 2.2e-16
sqrt(182.7 / 288) # 0.7964766

icc_pval <- c(0.252, 0.374, 0.518, 0.614, 0.515, 0.108, 0.0802, 0.475, 0.649, 0.468)
p.adjust(icc_pval, method = "bonferroni")

# Final matrix with all the different correlations
AllCorr <- data.frame(
    COHORT = c(1, 0.1999, 0.06545, 0.1060, 0, 0, 0.009701, 0.8613, 0.8813),
    SEX = c(0, 1, 0, 0.8233, 0.03750, 0, 0.001201, 0.1868, 0.2255),
    AGE = c(0, 0, 1, 0, 0, 0, 0.188337165, 0.2347945, 0.3040239),
    HGHT = c(0, 0, 0, 1, 0.68294460, 0.18884719, 0.24466244, 0.2081583, 0.2862412),
    WGHT = c(0, 0, 0, 0, 1, 0.82440147, 0.18198056, 0, 0.26351),
    BMI = c(0, 0, 0, 0, 0, 1, 0, 0, 0),
    TRISCHD = c(0, 0, 0, 0, 0, 0, 1, 0.7757837, 0.7964766),
    DTHVNT = c(0, 0, 0, 0, 0, 0, 0, 1, 0.6855),
    DTHHRDY = c(0, 0, 0, 0, 0, 0, 0, 0, 1)
)
rownames(AllCorr) <- c("COHORT", "SEX", "AGE", "HGHT", "WGHT", "BMI", "TRISCHD", "DTHVNT", "DTHHRDY")
corrplot(as.matrix(AllCorr))
