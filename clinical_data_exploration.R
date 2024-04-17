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


