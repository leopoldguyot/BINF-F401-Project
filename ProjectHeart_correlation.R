setwd('C:/Users/Draguet/Desktop/Cours/MA-BINF 2023-2024/BINF-F401 Computational Methods for Functional Genomics/project-BINF-F401-2024/Heart')

library(ggplot2)
library(corrplot)

clinical <- read.table(file = 'clinical_data.tsv', sep = '\t', header = TRUE)
morpho <- read.table(file = 'morphological_counts_lunit_dino.tsv', sep = '\t', header = TRUE)
RNAcount <- read.table(file = 'RNA_read_counts.tsv', sep = '\t', header = TRUE)

#Morphological atlas with the different clusters
#https://www-hpda.ulb.ac.be/iribhm/ai/atlases/lunit_dino/Heart/graphics/sankey.html 


#Q1 Explore clinical variables

#AGE = Age at death
#SEX = 1 (Male) & 2 (Female)
#HGHT = height of donor (inch)
#WGHT = weight of donor (pounds)
#BMI = general body fat indicator (ratio of weight to height)
#COHORT = Organ donor & Postmortem
#TRISCHD = Ischemic Time in minutes (time btw death and tissue collection)
#DTHHRDY = Hardy scale (type of death)
  #0 = Ventilator Case
  #1 = Violent and fast death
  #2 = Fast death of natural causes
  #3 = Intermediate death
  #4 = Slow death
#DTHVNT Was donor on a ventilator immediately prior to Death.
  #0 = No
  #1 = Yes
  #98 = Not Reported
  #99 = Unknown  

#Correlation
val_matrix <- cbind(clinical[3:8], clinical[10:12])
val_matrix$COHORT <- as.numeric(factor(val_matrix$COHORT))
M_corr <- cor(val_matrix)
corrplot(M_corr, method="circle")


#Q1.2. Are the clinical variables correlated ?

library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)

#Checking for null values
colSums(is.na(clinical)) #no null values

#Normalizing the data
quant_variables <- cbind(clinical[5:8], clinical[10])
data_normalized <- scale(quant_variables)

#Computing correlation matrix
corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)

#Applying PCA
data_pca <- princomp(corr_matrix)
summary(data_pca)

#Visualization of the principal components
fviz_eig(data_pca, addlabels = TRUE)
fviz_pca_var(data_pca, col.var = "black")
fviz_cos2(data_pca, choice = "var", axes = 1:2)
fviz_pca_var(data_pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)



#Considering every variables Factorial Analysis of Mixed Data 
val_matrix <- cbind(clinical[3:8], clinical[10:12])

sex_genre <- c()
for (i in clinical$SEX){
  if (i == 1){sex_genre <- c(sex_genre, 'Male')}
  if (i == 2){sex_genre <- c(sex_genre, 'Female')}
}
val_matrix$SEX <- sex_genre

vent <- c()
for (i in clinical$DTHVNT){
  if (i == 0){vent <- c(vent, 'No')}
  if (i == 1){vent <- c(vent, 'Yes')}
  if (i == 98){vent <- c(vent, 'Not Reported')}
  if (i == 99){vent <- c(vent, 'Unknown')}
}
val_matrix$DTHVNT <- vent

hardy_scale <- c()
for (i in clinical$DTHHRDY){
  if (i == 0){hardy_scale <- c(hardy_scale, 'Ventilator Case')}
  if (i == 1){hardy_scale <- c(hardy_scale, 'Violent and Fast Death')}
  if (i == 2){hardy_scale <- c(hardy_scale, 'Fast Death of Natural Causes')}
  if (i == 3){hardy_scale <- c(hardy_scale, 'Intermediate Death')}
  if (i == 4){hardy_scale <- c(hardy_scale, 'Slow Death')}
}
val_matrix$DTHHRDY <- hardy_scale

clinical_famd <- FAMD(val_matrix, ncp = 10, graph = TRUE)

get_eigenvalue(clinical_famd)
fviz_screeplot(clinical_famd) #scree plot (the percentages of inertia explained by each FAMD dimensions)
fviz_famd_var(clinical_famd, repel=TRUE) #to plot both quantitative and qualitative variables
fviz_contrib(clinical_famd, "var", axes = 1) #to visualize the contribution of variables to the principal dimensions
fviz_contrib(clinical_famd, "var", axes = 2)
fviz_famd_var(clinical_famd, "quanti.var", repel = TRUE, col.var = "black") #correlation circle for quantitative data
fviz_famd_var(clinical_famd, "quali.var", col.var = "black", repel=TRUE)


#correlation matrix with mixed-data (or maybe need to make ANOVA, chi-square and Pearson's correlation test)
library(ggcorrplot)
test <- model.matrix(~0+., data=val_matrix)
test_corr <- cor(test, use="pairwise.complete.obs")
ggcorrplot(test_corr, show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)


#Chi-square between categories
val_matrix_categ <- cbind(clinical[3:4], clinical[11:12])

test_data <- table(val_matrix$COHORT, val_matrix$SEX)
print(test_data)
test_chisq <- chisq.test(table(val_matrix$COHORT, val_matrix$SEX))
test_chisq$p.value

col_names <- colnames(val_matrix_categ)
chisq_pval <- c()
chisq_names <- c()
for (i in 1:3){
  for (j in i+1:4){
    if (j < 5){
      chisq_calc <- chisq.test(table(val_matrix_categ[[i]], val_matrix_categ[[j]]), simulate.p.value = TRUE)
      chisq_pval <- c(chisq_pval, chisq_calc$p.value)
      chisq_names <- c(chisq_names, paste(col_names[i], 'x', col_names[j]))
    }
  }
}
corrected_pval <- p.adjust(chisq_pval, method = 'bonferroni')
results_categ <- data.frame(Comparaison = chisq_names,
                            Pval = corrected_pval)
results_categ #devrait pas avoir correlation entre sex et les autres ?


#ANOVA between categ and continuous, dont know how to actually do something
test_aov <- table()
val_matrix_cont <- cbind(clinical[5:8], clinical[10])

#Levene Test of homogeneity, maybe not useful
library(car)
leveneTest(AGE ~ COHORT, data = val_matrix)

test_aov <- aov(AGE ~ COHORT, data = val_matrix)
summary(test_aov) #signif difference
leveneTest(log(val_matrix$AGE), factor(val_matrix$COHORT), data=val_matrix)

model  <- lm(AGE ~ COHORT, data = val_matrix)
shapiro.test(residuals(model))
