library(ggplot2)
library(ggcorrplot)

clusters_df <- read.table("data/morphological_counts_lunit_dino.tsv", header=TRUE, sep="\t")
rownames(clusters_df) <- clusters_df$SMPLID
clusters_df$SMPLID <- NULL 
summary(clusters_df)

corr_matrix <- cor(clusters_df)
corr_plot <- ggcorrplot(corr_matrix)
ggsave("figures/clusters_correlation/corrplot.pdf")
corrplot(corr_matrix)
