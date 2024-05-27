library(ggplot2)
library(ggcorrplot)
library(corrplot)

clusters_df <- read.table("data/morphological_counts_lunit_dino.tsv", header=TRUE, sep="\t")
rownames(clusters_df) <- clusters_df$SMPLID
clusters_df$SMPLID <- NULL 
summary(clusters_df)

clusters_df <- log10(clusters_df + 1)
corr_matrix <- cor(clusters_df)
corr_plot <- ggcorrplot(corr_matrix)
ggsave("figures/clusters_correlation/corrplot.pdf")
corrplot(corr_matrix)


ggplot(clusters_df, aes(x=Mophological.cluster.G4_23, y=Mophological.cluster.G4_25)) + geom_point() +geom_point()
