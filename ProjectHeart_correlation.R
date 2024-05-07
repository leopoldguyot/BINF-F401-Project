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
M_corr <- cor(val_matrix, method = 'spearman')
corrplot(M_corr, method="circle")

p_val_M_corr <- c(0.6228, 0.7948, 0.9113, 0.001322, 2.2e-16, 0.001282, 2.693e-05, 
                  2.2e-16, 0.00193, 0.5186)
p.adjust(p_val_M_corr, method='bonferroni') #stay the same
cor.test(val_matrix$AGE,val_matrix$HGHT, method = 'spearman') #not significant
cor.test(val_matrix$AGE,val_matrix$WGHT, method = 'spearman') #not significant
cor.test(val_matrix$AGE,val_matrix$BMI, method = 'spearman') #not significant
cor.test(val_matrix$AGE,val_matrix$TRISCHD, method = 'spearman') #significant

cor.test(val_matrix$HGHT,val_matrix$WGHT, method = 'spearman') #significant
cor.test(val_matrix$HGHT,val_matrix$BMI, method = 'spearman') #significant
cor.test(val_matrix$HGHT,val_matrix$TRISCHD, method = 'spearman') #significant

cor.test(val_matrix$WGHT,val_matrix$BMI, method = 'spearman') #significant
cor.test(val_matrix$WGHT,val_matrix$TRISCHD, method = 'spearman') #significant

cor.test(val_matrix$BMI,val_matrix$TRISCHD, method = 'spearman') # not significant

#Q1.2. Are the clinical variables correlated ?

#library(corrr)
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
fviz_famd_var(clinical_famd, "quanti.var", repel = TRUE, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) #correlation circle for quantitative data
fviz_famd_var(clinical_famd, "quali.var", col.var = "black", repel=TRUE)
fviz_cos2(clinical_famd, choice = "var", axes = 1:2)

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





#Correlation between binary variables and quantitative variables

shap_pval <- c()
shap_name <- c()
wilcox_pval <- c()
wilcox_name <- c()

#Sex and Age

res <- shapiro.test(val_matrix$AGE[val_matrix$SEX=="Male"]) #not a normal distribution
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "AgeMale")
res <- shapiro.test(val_matrix$AGE[val_matrix$SEX=="Female"]) #not normal
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "AgeFemale")

res <- wilcox.test(AGE~SEX, data=val_matrix) #pval = 0.1524 => H0: pas de différence entre les deux groupes
wilcox_pval <- c(wilcox_pval, res$p.value) #pas de corrélation entre l'âge et le sexe
wilcox_name <- c(wilcox_name, 'AgeSex')
#If we consider that the distrib are normal because central limit theorem
t.test(AGE~SEX, data=val_matrix) #0.06 pas de réelle différence entre les deux moyennes

#Sex and Height

res <- shapiro.test(val_matrix$HGHT[val_matrix$SEX=="Male"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "HeightMale")
res <- shapiro.test(val_matrix$HGHT[val_matrix$SEX=="Female"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "HeightFemale")

res <- wilcox.test(HGHT~SEX, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'HeightSex')

#Sex and Weight

res <- shapiro.test(val_matrix$WGHT[val_matrix$SEX=="Male"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "WeightMale")
res <- shapiro.test(val_matrix$WGHT[val_matrix$SEX=="Female"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "WeightFemale")

res <- wilcox.test(WGHT~SEX, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'WeightSex')

#Sex and BMI

res <- shapiro.test(val_matrix$BMI[val_matrix$SEX=="Male"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "BMIMale")
res <- shapiro.test(val_matrix$BMI[val_matrix$SEX=="Female"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "BMIFemale")

res <- wilcox.test(BMI~SEX, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'BMISex')

#Sex and Trisch

res <- shapiro.test(val_matrix$TRISCHD[val_matrix$SEX=="Male"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "TRISCHDMale")
res <- shapiro.test(val_matrix$TRISCHD[val_matrix$SEX=="Female"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "TRISCHDFemale")

res <- wilcox.test(TRISCHD~SEX, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'TRISCHDSex')

#Cohort and Age

res <- shapiro.test(val_matrix$AGE[val_matrix$COHORT=="Postmortem"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "AgePM")
res <- shapiro.test(val_matrix$AGE[val_matrix$COHORT=="Organ Donor (OPO)"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "AgeOD")

res <- wilcox.test(AGE~COHORT, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'AgeCohort')

#Sex and Height

res <- shapiro.test(val_matrix$HGHT[val_matrix$COHORT=="Postmortem"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "HeightPM")
res <- shapiro.test(val_matrix$HGHT[val_matrix$COHORT=="Organ Donor (OPO)"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "HeightOD")

res <- wilcox.test(HGHT~COHORT, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'HeightCohort')

#Sex and Weight

res <- shapiro.test(val_matrix$WGHT[val_matrix$COHORT=="Postmortem"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "WeightPM")
res <- shapiro.test(val_matrix$WGHT[val_matrix$COHORT=="Organ Donor (OPO)"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "WeightOD")

res <- wilcox.test(WGHT~COHORT, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'WeightCohort')

#Sex and BMI

res <- shapiro.test(val_matrix$BMI[val_matrix$COHORT=="Postmortem"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "BMIPM")
res <- shapiro.test(val_matrix$BMI[val_matrix$COHORT=="Organ Donor (OPO)"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "BMIOD")

res <- wilcox.test(BMI~COHORT, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'BMICohort')

#Sex and Trisch

res <- shapiro.test(val_matrix$TRISCHD[val_matrix$COHORT=="Postmortem"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "TRISCHDPM")
res <- shapiro.test(val_matrix$TRISCHD[val_matrix$COHORT=="Organ Donor (OPO)"])
shap_pval <- c(shap_pval, res$p.value)
shap_name <- c(shap_name, "TRISCHDOD")

res <- wilcox.test(TRISCHD~COHORT, data=val_matrix)
wilcox_pval <- c(wilcox_pval, res$p.value)
wilcox_name <- c(wilcox_name, 'TRISCHDCohort')

corrected_shap <- p.adjust(shap_pval, method = 'bonferroni')
corrected_wilc <- p.adjust(wilcox_pval, method = 'bonferroni')
shap_data <- data.frame(pval = shap_pval,
                        pval_corrected = corrected_shap,
                        names = shap_name)
#only weight variable has a normal distribution in the different groups (pval > 0.05)
t.test(WGHT~SEX, data=val_matrix) #means are different, then some correlation
t.test(WGHT~COHORT, data=val_matrix) #means not different, no correlation

wilcox_data <- data.frame(quant_var = c('Age', 'Height', 'Weight', 'BMI', 'Trisch'),
                          sex = corrected_wilc[1:5],
                          cohort = corrected_wilc[6:10])
#H0 => no significant differences between the two groups
#Thus, there is a correlation between the Sex variable and HGHT, WGHT and TRISCH
#Also correlation between COHORT and AGE, HGHT and TRISCHT


#Correlation between categorical variables => Cramer's V
library(rcompanion)

chisq.test(table(val_matrix$SEX, val_matrix$COHORT)) #p-val ptet à corriger, ici sont bien dépendantes
cramerV(table(val_matrix$SEX, val_matrix$COHORT), bias.correct = TRUE) #correlation value, allowed because base on X2, non-param
#0.1999
#0 - 0.2 -> weak association
#0.2 - 0.6 -> middle association
#0.6 - 1 -> strong association

chisq.test(table(val_matrix$SEX, val_matrix$DTHVNT)) #may be incorrect because some cells are very low => only consider No and Yes
chisq.test(val_matrix$SEX[val_matrix$DTHVNT!='Unknown'],
           val_matrix$DTHVNT[val_matrix$DTHVNT!='Unknown'])
DTHVNT_remov <- as.data.frame.matrix(table(val_matrix$SEX, val_matrix$DTHVNT))
DTHVNT_remov[2] <- NULL
cramerV(as.matrix(DTHVNT_remov), bias.correct = TRUE) #0.1868

chisq.test(table(val_matrix$SEX, val_matrix$DTHHRDY)) #may delete last column
chisq.test(val_matrix$SEX[val_matrix$DTHHRDY!="Violent and Fast Death"],
           val_matrix$DTHHRDY[val_matrix$DTHHRDY!="Violent and Fast Death"])
DTHHRDY_table <- as.data.frame.matrix(table(val_matrix$SEX, val_matrix$DTHHRDY))
DTHHRDY_table[5] <- NULL
cramerV(as.matrix(DTHHRDY_table), bias.correct = TRUE) #0.2255

chisq.test(table(val_matrix$COHORT, val_matrix$DTHVNT)) #may be incorrect
chisq.test(val_matrix$COHORT[val_matrix$DTHVNT!='Unknown'],
           val_matrix$DTHVNT[val_matrix$DTHVNT!='Unknown'])
cramerV(table(val_matrix$COHORT,val_matrix$DTHVNT), bias.correct = TRUE) #0.8613

chisq.test(table(val_matrix$COHORT, val_matrix$DTHHRDY)) #may be incorrect
chisq.test(val_matrix$COHORT[val_matrix$DTHHRDY!="Violent and Fast Death"],
           val_matrix$DTHHRDY[val_matrix$DTHHRDY!="Violent and Fast Death"])
DTHHRDY_table2 <- as.data.frame.matrix(table(val_matrix$COHORT, val_matrix$DTHHRDY))
DTHHRDY_table2[5] <- NULL
cramerV(as.matrix(DTHHRDY_table2), bias.correct = TRUE) #0.8813

chisq.test(table(val_matrix$DTHVNT, val_matrix$DTHHRDY)) #may be incorrect
n_DTHVNT <- val_matrix$DTHVNT[val_matrix$DTHHRDY!="Violent and Fast Death"]
n_DTHHRDY <- val_matrix$DTHHRDY[val_matrix$DTHHRDY!="Violent and Fast Death"]
chisq.test(n_DTHVNT[n_DTHVNT!="Unknown"],
           n_DTHHRDY[n_DTHVNT!="Unknown"])
DTHHRDY_table <- as.data.frame.matrix(table(val_matrix$COHORT, val_matrix$DTHHRDY))
DTHHRDY_table[5] <- NULL
cramerV(table(n_DTHVNT, n_DTHHRDY), bias.correct = TRUE) #0.6855

chisq_pval <- c(0.0006642102,0.001469577,0.0006295674,1.727819e-47,1.589829e-47,6.091008e-57)
corr_chisq <- p.adjust(chisq_pval, method = 'bonferroni') #still below 0.05
cramerCorr <- data.frame(SEX = c(1, 0.1999, 0.1868,0.2255),
                         COHORT = c(0.1999,1,0.8613,0.8813),
                         DTHVNT = c(0.1868,0.8613,1,0.6855),
                         DTHHRDY = c(0.2255, 0.8813, 0.6855, 1))
rownames(cramerCorr) <- c('SEX', 'COHORT', 'DTHVNT', 'DTHHRDY')
corrplot(as.matrix(cramerCorr))


#logistic regression with binomial variables and quantitatives variables
summary(glm(as.factor(SEX) ~ AGE, family=binomial, data=val_matrix)) #0.02397 slope not significant
summary(glm(as.factor(SEX) ~ HGHT, family=binomial, data=val_matrix)) #0.8233
summary(glm(as.factor(SEX) ~ WGHT, family=binomial, data=val_matrix)) #0.037508
summary(glm(as.factor(SEX) ~ BMI, family=binomial, data=val_matrix)) #0.06833 #not significant
summary(glm(as.factor(SEX) ~ TRISCHD, family=binomial, data=val_matrix)) #0.0012615
summary(glm(as.factor(COHORT) ~ AGE, family=binomial, data=val_matrix)) #0.06545 slope
summary(glm(as.factor(COHORT) ~ HGHT, family=binomial, data=val_matrix)) #0.10607
summary(glm(as.factor(COHORT) ~ WGHT, family=binomial, data=val_matrix)) # 0.006206 #not significant
summary(glm(as.factor(COHORT) ~ BMI, family=binomial, data=val_matrix)) #0.002597 not significant
summary(glm(as.factor(COHORT) ~ TRISCHD, family=binomial, data=val_matrix)) #0.009701
test <- c(0.0448, 6.15e-16, 1.7e-13, 0.0286, 0.00024, 6.82e-07, 0.000683, 0.0522, 0.928, 2e-16)
test <- p.adjust(test, method = 'bonferroni')
test

log_reg <- data.frame(SEX = c(0, 0.8233, 0.037508, 0, 0.0012615),
                      COHORT = c(0.06545,0.10607,0,0,0.009701))
row.names(log_reg) <- c('AGE', 'HGHT', 'WGHT', 'BMI', 'TRSCHD')
corrplot(as.matrix(log_reg))

#Intraclass correlation coefficient (between categ and quanti)
library(irr)

kruskal.test(AGE ~ DTHVNT, data = val_matrix) #0.0003567, not the same distrib per group 
#using corr coef calc (sqrt(X2 / n_obs)) coef = 0.4815845
sqrt(15.877/288) #0.2347945

kruskal.test(HGHT ~ DTHVNT, data = val_matrix) #0.001951
sqrt(12.479/288) #0.2081583

kruskal.test(WGHT ~ DTHVNT, data = val_matrix) #0.07141

kruskal.test(BMI ~ DTHVNT, data = val_matrix) #0.9119

kruskal.test(TRISCHD ~ DTHVNT, data = val_matrix) # 2.2e-16
sqrt(173.33/288) #0.7757837

kruskal.test(AGE ~ DTHHRDY, data = val_matrix) # 2.372e-05
sqrt(26.62/288) #0.3040239

kruskal.test(HGHT ~ DTHHRDY, data = val_matrix) # 9.619e-05
sqrt(23.597/288) #0.2862412

kruskal.test(WGHT ~ DTHHRDY, data = val_matrix) # 0.0004999
sqrt(19.998/288) #0.26351

kruskal.test(BMI ~ DTHHRDY, data = val_matrix) # 0.1106

kruskal.test(TRISCHD ~ DTHHRDY, data = val_matrix) # 2.2e-16
sqrt(182.7/288) #0.7964766

icc_pval <- c(0.252, 0.374, 0.518, 0.614, 0.515, 0.108, 0.0802, 0.475, 0.649, 0.468)
p.adjust(icc_pval, method='bonferroni')

#Final matrix with all the different correlations
AllCorr <- data.frame(COHORT = c(1, 0.1999, 0.06545,0.1060,0,0,0.009701,0.8613,0.8813),
                      SEX = c(0,1,0,0.8233,0.03750,0,0.001201,0.1868,0.2255),
                      AGE = c(0,0,1,0,0,0,0.188337165,0.2347945,0.3040239),
                      HGHT = c(0,0,0,1,0.68294460, 0.18884719, 0.24466244, 0.2081583,0.2862412),
                      WGHT = c(0,0,0,0,1, 0.82440147, 0.18198056, 0,0.26351),
                      BMI = c(0,0,0,0,0, 1, 0, 0,0),
                      TRISCHD = c(0,0,0,0,0, 0, 1, 0.7757837,0.7964766),
                      DTHVNT = c(0,0,0,0,0, 0, 0, 1,0.6855),
                      DTHHRDY = c(0,0,0,0,0, 0, 0, 0,1))
rownames(AllCorr) <- c('COHORT', 'SEX','AGE','HGHT','WGHT','BMI','TRISCHD', 'DTHVNT', 'DTHHRDY')
corrplot(as.matrix(AllCorr))
