######### Packages Import #########
library("ggplot2")
library("GGally")
library("corrplot")
library("reshape2")

######### Functions Definition #########
distribution_plot <- function(data, variable, title){
  if(is.numeric(data[[variable]])) {
    plot <- ggplot(data, aes(x = .data[[variable]])) +
      geom_histogram(fill = "skyblue", color = "black")
  } else {
    plot <- ggplot(data, aes(x = .data[[variable]])) +
      geom_bar(fill = "skyblue", color = "black")
  }

  plot + labs(title = title,
              x = variable,
              y = "Frequency")
}

######### Body #########
clinical_data <- read.table("data/clinical_data.tsv",
 sep = "\t",
 header = TRUE)

clinical_data$SEX <- as.factor(clinical_data$SEX)
clinical_data$DTHVNT <- as.factor(clinical_data$DTHVNT)
clinical_data$DTHHRDY <- as.factor(clinical_data$DTHHRDY)
clinical_data$COHORT <- as.factor(clinical_data$COHORT)
clinical_data$SMPTHNTS <- as.factor(clinical_data$SMPTHNTS)


plot <- lapply(colnames(clinical_data)[-c(1:2, 13:14)], function(x){
    plot <- distribution_plot(clinical_data, x, x)
    path <- file.path("figures", "clinical_data_var_dist", paste0(x, ".pdf"))
    ggsave(filename = path,
    plot,
     width = 10,
      height = 10)
    })

ggpairs_plot <- GGally::ggpairs(clinical_data[-c(1,2,13,14)])
ggsave(filename = file.path("figures",
                            "clinical_data_exploration",
                            "ggpairs_plot.pdf"),
       ggpairs_plot,
       width = 15,
       height = 10)

numeric_clinical_data <- apply(clinical_data[-c(1:2, 13:14)],
                               2,
                               as.numeric)

cor_clinical_data <- cor(numeric_clinical_data)
cov_clinical_data <- cov(numeric_clinical_data)

corrplot(cor_clinical_data[-c(1,7),-c(1,7)], method = "circle")
melted_cor <- melt(cor_clinical_data)
melted_cov <- melt(cov_clinical_data)
ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()
ggplot(data = melted_cov, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()

# mise a l'echelle de certaines variable 
clinical_data_scale <- data.frame(clinical_data$SMPLID,clinical_data$SUBJID,
                                  clinical_data$COHORT,clinical_data$SMPTHNTS,
                                  clinical_data$SMPLID.1,clinical_data$IMGURL,
                                  clinical_data$SEX,clinical_data$DTHHRDY,
                                  clinical_data$DTHVNT)
scaleAGE <- scale(clinical_data$AGE)
scaleHGHT <- scale(clinical_data$HGHT)
scaleWGHT <- scale(clinical_data$WGHT)
scaleBMI <- scale(clinical_data$BMI)
scaleTRISCHD <- scale(clinical_data$TRISCHD)

clinical_data_scale$AGE_scale <- scaleAGE 
clinical_data_scale$AGE <- clinical_data$AGE
clinical_data_scale$HGHT_scale <- scaleHGHT
clinical_data_scale$WGHT_scale <- scaleWGHT
clinical_data_scale$BMI_scale <- scaleBMI
clinical_data_scale$TRISCHD_scale <- scaleTRISCHD
hist(scaleBMI)

plot <- lapply(colnames(clinical_data_scale)[-c(1:2, 5:6)], function(x){
  plot <- distribution_plot(clinical_data_scale, x, x)
  path <- file.path("figures", "clinical_data_scale_var_dist", paste0(x, ".pdf"))
  ggsave(filename = path,
         plot,
         width = 10,
         height = 10)
})
# transformation Racine carré 
clinical_data_racine_caree <- data.frame(clinical_data$SMPLID,clinical_data$SUBJID,
                                  clinical_data$COHORT,clinical_data$SMPTHNTS,
                                  clinical_data$SMPLID.1,clinical_data$IMGURL,
                                  clinical_data$SEX,clinical_data$DTHHRDY,
                                  clinical_data$DTHVNT)
racineAGE <- sqrt(clinical_data$AGE)
racineHGHT <- sqrt(clinical_data$HGHT)
racineWGHT <- sqrt(clinical_data$WGHT)
racineBMI <- sqrt(clinical_data$BMI)
racineTRISCHD <- sqrt(clinical_data$TRISCHD)
clinical_data_racine_caree$AGE_racine <- racineAGE 
clinical_data_racine_caree$HGHT_racine <- racineHGHT
clinical_data_racine_caree$WGHT_racine <- racineWGHT
clinical_data_racine_caree$BMI_racine <- racineBMI
clinical_data_racine_caree$TRISCHD_racine <- racineTRISCHD

plot <- lapply(colnames(clinical_data_racine_caree)[-c(1:2, 5:6)], function(x){
  plot <- distribution_plot(clinical_data_racine_caree, x, x)
  path <- file.path("figures", "clinical_data_racine_caree_var_dist", paste0(x, ".pdf"))
  ggsave(filename = path,
         plot,
         width = 10,
         height = 10)
})
#transoformation log 
clinical_data_log <- data.frame(clinical_data$SMPLID,clinical_data$SUBJID,
                                         clinical_data$COHORT,clinical_data$SMPTHNTS,
                                         clinical_data$SMPLID.1,clinical_data$IMGURL,
                                         clinical_data$SEX,clinical_data$DTHHRDY,
                                         clinical_data$DTHVNT)
logAGE <- log(clinical_data$AGE)
logHGHT <- log(clinical_data$HGHT)
logWGHT <- log(clinical_data$WGHT)
logBMI <- log(clinical_data$BMI)
logTRISCHD <- log(clinical_data$TRISCHD)
clinical_data_log$AGE_log <- logAGE 
clinical_data_log$HGHT_log <- logHGHT
clinical_data_log$WGHT_log <- logWGHT
clinical_data_log$BMI_log <- logBMI
clinical_data_log$TRISCHD_log <- logTRISCHD

plot <- lapply(colnames(clinical_data_log)[-c(1:2, 5:6)], function(x){
  plot <- distribution_plot(clinical_data_log, x, x)
  path <- file.path("figures", "clinical_data_log_var_dist", paste0(x, ".pdf"))
  ggsave(filename = path,
         plot,
         width = 10,
         height = 10)
})
# Inverse Normal Transformation (rank based)
clinical_data_Inverse_Normal_Transformation <- clinical_data[c("SMPLID",
 "SUBJID",
  "COHORT",
   "SMPTHNTS",
    "SMPLID.1",
     "IMGURL",
      "SEX",
       "DTHHRDY",
        "DTHVNT")]

INT_AGE <-  qnorm((rank(clinical_data$AGE,na.last="keep")-0.375)/(sum(!is.na(clinical_data$AGE))-2*0.375+1))
INT_HGHT <-  qnorm((rank(clinical_data$HGHT,na.last="keep")-0.375)/(sum(!is.na(clinical_data$HGHT))-2*0.375+1))
INT_WGHT <-  qnorm((rank(clinical_data$WGHT,na.last="keep")-0.5)/(sum(!is.na(clinical_data$WGHT))-2*0.375+1))
INT_BMI <-  qnorm((rank(clinical_data$BMI,na.last="keep")-0.5)/(sum(!is.na(clinical_data$BMI))-2*0.375+1))
INT_TRISCHD <-  qnorm((rank(clinical_data$TRISCHD,na.last="keep")-0.5)/(sum(!is.na(clinical_data$TRISCHD))-2*0.375+1))



clinical_data_Inverse_Normal_Transformation$AGE_INT <- INT_AGE 
clinical_data_Inverse_Normal_Transformation$HGHT_INT <- INT_HGHT
clinical_data_Inverse_Normal_Transformation$WGHT_INT <- INT_WGHT
clinical_data_Inverse_Normal_Transformation$BMI_INT <- INT_BMI
clinical_data_Inverse_Normal_Transformation$TRISCHD_INT <- INT_TRISCHD
hist(clinical_data_Inverse_Normal_Transformation$AGE_INT)

plot <- lapply(colnames(clinical_data_Inverse_Normal_Transformation)[-c(1:2, 5:6)], function(x){
  plot <- distribution_plot(clinical_data_Inverse_Normal_Transformation, x, x)
  path <- file.path("figures", "clinical_data_Inverse_Normal_Transformation_var_dist", paste0(x, ".pdf"))
  ggsave(filename = path,
         plot,
         width = 10,
         height = 10)
})


# test de normalité 
shapiro.test(clinical_data$AGE)
shapiro <- lapply(clinical_data, function(x){
  if(is.numeric(x)){
    shapiro.test(x)
  }
})
names(shapiro) <- names(clinical_data)
shapiro
#test de normalité sur donné transformé avec racine 
shapiro_racine <- lapply(clinical_data_racine_caree, function(x){
  if(is.numeric(x)){
    shapiro.test(x)
  }
})
shapiro_racine
#test de normalité sur donnée transformé avec log 
shapiro_log <- lapply(clinical_data_log, function(x){
  if(is.numeric(x)){
    shapiro.test(x)
  }
})
shapiro_log
#test de normalité sur donnée scale 
shapiro_scale <- lapply(clinical_data_scale, function(x){
  if(is.numeric(x)){
    shapiro.test(x)
  }
})
shapiro_scale
# test de normalité  sur INT 
shapiro_INT <- lapply(clinical_data_Inverse_Normal_Transformation, function(x){
  if(is.numeric(x)){
    shapiro.test(x)
  }
})
shapiro_INT
# save clinical_data_Inverse_Normal_Transformation 
write.table(clinical_data_Inverse_Normal_Transformation, "data_output/clinical_data_normalized.tsv", sep="\t")
test <- clinical_data_test <- read.table("data/clinical_data_Inverse_Normal_Transformation.tsv",
                                    sep = "\t",
                                    header = TRUE)
