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

# test de normalité 
shapiro.test(clinical_data$AGE)
shapiro <- lapply(clinical_data, function(x){
  if(is.numeric(x)){
    shapiro.test(x)
  }
})
names(shapiro) <- names(clinical_data)
shapiro
