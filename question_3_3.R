library(ggplot2)

file_names <- list.files("~/data_output/significant_q3")

sigTSV <- read.table(file = 'significant_count.tsv', sep = '\t', header = TRUE)
sigCSV <- rbind(test5, c('Total', sum(test5$count_non_confounders), sum(test5$count_confounders)))
write.csv(sigCSV, file = "~/data_output/significant_q3/significant_total.csv", row.names = FALSE)


# Genes without considering confounders

gene_list <- c()
occ_list <- c()

for (i in file_names){
  if(!(grepl('confounders', i))){
    table <- read.table(file = i, sep = '\t', header = TRUE)
    for (j in table$gene){
      if (j %in% gene_list){
        occ_list[which(gene_list == j)] <- occ_list[which(gene_list == j)] + 1}
      else{
        gene_list <- c(gene_list, j)
        occ_list <- c(occ_list, 1)}
    }
  }
}

nonConfounders <- data.frame(Genes = gene_list[occ_list > 1],
                             Occurences = occ_list[occ_list > 1])

ggplot(data=nonConfounders, aes(x=Genes, y=Occurences)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  ggtitle('Occurences of significant genes in the clusters')


# Genes with confounders

gene_list <- c()
occ_list <- c()

for (i in file_names){
  if(grepl('confounders', i)){
    table <- read.table(file = i, sep = '\t', header = TRUE)
    for (j in table$gene){
      if (j %in% gene_list){
        occ_list[which(gene_list == j)] <- occ_list[which(gene_list == j)] + 1}
      else{
        gene_list <- c(gene_list, j)
        occ_list <- c(occ_list, 1)}
    }
  }
}

confounders <- data.frame(Genes = gene_list[occ_list > 1],
                             Occurences = occ_list[occ_list > 1])

ggplot(data=confounders, aes(x=Genes, y=Occurences)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  ggtitle('Significant genes considering confounders effect')


#Counting similarities between clusters

all_sim <- c()

for (i in 0:31){
  fileA <- paste('Mophological.cluster.G4_', i, '_confounders.tsv', sep = '')
  tableA <- read.table(file = fileA, sep = '\t', header = TRUE)
  simil_list <- c()
  
  for (j in 0:31){
    fileB <- paste('Mophological.cluster.G4_', j, '_confounders.tsv', sep = '')
    tableB <- read.table(file = fileB, sep = '\t', header = TRUE)
    n_simil <- 0
    
    for (k in tableA$gene){
      if (k %in% tableB$gene){
        n_simil <- n_simil + 1
      }
    }
    simil_list <- c(simil_list, n_simil)
  }
  all_sim <- c(all_sim, simil_list)
}

test <- c(0:31)
test_list <- c()
for(i in 0:31){
  test_list <- c(test_list, test)
}
test_listOrd <- sort(test_list)

confSimdf <- data.frame(res = all_sim,
                        Var1 = test_list,
                        Var2 = test_listOrd)

ggplot(confSimdf, aes(x = factor(Var1), y = factor(Var2), fill = res)) +
  geom_tile(color = "white") + # Create heatmap
  scale_fill_gradient(low = "white", high = "steelblue") + # Set color gradient
  theme_minimal() + # Set theme
  labs(x = "", y = "", title = "Number of similar genes per cluster and considering confounders") + # Labels
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + # Rotate x-axis labels
  geom_text(aes(label = res), color = "black", size = 3)


#Similarities without considering confounders

all_sim <- c()

for (i in 0:31){
  fileA <- paste('Mophological.cluster.G4_', i, '.tsv', sep = '')
  tableA <- read.table(file = fileA, sep = '\t', header = TRUE)
  simil_list <- c()
  
  for (j in 0:31){
    fileB <- paste('Mophological.cluster.G4_', j, '.tsv', sep = '')
    tableB <- read.table(file = fileB, sep = '\t', header = TRUE)
    n_simil <- 0
    
    for (k in tableA$gene){
      if (k %in% tableB$gene){
        n_simil <- n_simil + 1
      }
    }
    simil_list <- c(simil_list, n_simil)
  }
  all_sim <- c(all_sim, simil_list)
}

test <- c(0:31)
test_list <- c()
for(i in 0:31){
  test_list <- c(test_list, test)
}
test_listOrd <- sort(test_list)

confSimdf <- data.frame(res = all_sim,
                        Var1 = test_list,
                        Var2 = test_listOrd)

ggplot(confSimdf, aes(x = factor(Var1), y = factor(Var2), fill = res)) +
  geom_tile(color = "white") + # Create heatmap
  scale_fill_gradient(low = "white", high = "steelblue") + # Set color gradient
  theme_minimal() + # Set theme
  labs(x = "", y = "", title = "Number of similar genes per cluster and without confounders") + # Labels
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + # Rotate x-axis labels
  geom_text(aes(label = res), color = "black", size = 3)


#Pathways

file_names <- list.files("~/data_output/table_reactome")
data <- read.table("Mophological.cluster.G4_0_conf_fgseaRes.txt",sep="\t",header=T)
data2 <- read.table("Mophological.cluster.G4_0_conf_fgseaRes.txt",sep="\t",header=T)

path_list <- c()
occ_list <- c()

for (i in file_names){
  if(grepl('conf', i)){
    table <- read.table(file = i, sep = '\t', header = TRUE)
    for (j in table$pathway[1:10]){
      path <- substr(j, 10, nchar(j))
      if (path %in% path_list){
        occ_list[which(path_list == path)] <- occ_list[which(path_list == path)] + 1}
      else{
        path_list <- c(path_list, path)
        occ_list <- c(occ_list, 1)}
    }
  }
}

pathConfounders <- data.frame(Pathways = gsub('_', ' ', path_list[occ_list > 3]),
                             Occurences = occ_list[occ_list > 3])

ggplot(data=pathConfounders, aes(x=Pathways, y=Occurences)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,9,by=1)) +
  ggtitle('Occurences of significant pathways in the clusters considering confonders effects') +
  coord_flip() +
  theme(plot.title = element_text(hjust = 1.3))


all_sim <- c()

for (i in 0:31){
  fileA <- paste('Mophological.cluster.G4_', i, '_conf_fgseaRes.txt', sep = '')
  tableA <- read.table(file = fileA, sep = '\t', header = TRUE)
  simil_list <- c()
  
  for (j in 0:31){
    fileB <- paste('Mophological.cluster.G4_', j, '_conf_fgseaRes.txt', sep = '')
    tableB <- read.table(file = fileB, sep = '\t', header = TRUE)
    n_simil <- 0
    
    for (k in tableA$pathway[1:10]){
      if (k %in% tableB$pathway[1:10]){
        n_simil <- n_simil + 1
      }
    }
    simil_list <- c(simil_list, n_simil)
  }
  all_sim <- c(all_sim, simil_list)
}

test <- c(0:31)
test_list <- c()
for(i in 0:31){
  test_list <- c(test_list, test)
}
test_listOrd <- sort(test_list)

confSimdf <- data.frame(res = all_sim,
                        Var1 = test_list,
                        Var2 = test_listOrd)

ggplot(confSimdf, aes(x = factor(Var1), y = factor(Var2), fill = res)) +
  geom_tile(color = "white") + # Create heatmap
  scale_fill_gradient(low = "white", high = "steelblue") + # Set color gradient
  theme_minimal() + # Set theme
  labs(x = "", y = "", title = "Similar pathways between clusters considering confounders") + # Labels
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + # Rotate x-axis labels
  geom_text(aes(label = res), color = "black", size = 3)
