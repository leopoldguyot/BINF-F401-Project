#GSEA
#https://biostatsquid.com/fgsea-tutorial-gsea/
library('fgsea')
library('ggplot2')
library('data.table')


matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

prepare_gmt <- function(genes_data){
  # read gmt file 
  gmt <-gmtPathways("data/c2.cp.reactome.v7.5.1.symbols.gmt")
  # convert list pathways into a vector and without duplicate (vector of gene  )
  hidden <- unique(unlist(gmt))
  
  # convert into a matrix col = pathways and row = gene 
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  # check if hidden in gmt[[i]] (list of the gene for the i pathways ) put 1 if yes 
  
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  # take only the gene that are in genes_data 
  hidden1 <- intersect(genes_data, hidden)
  print(class(hidden1))
  # filter for gene sets with more than 5 genes annotated
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]]
  
  final_list <- matrix_to_list(mat)


}

GSEA <- function(name,transcript_count){
  # prepare bg_gene 
  genes_diff <- read.table(paste0("data_output/desq2_outputs_q3/",name,"_confounders.tsv"),
                           sep = "\t", header = TRUE, stringsAsFactors = TRUE)
  genes_diff$Name <- row.names(genes_diff)
  merged_df <- merge(genes_diff, transcript_count[, c("Description", "Name")], by = "Name", all.x = FALSE)
  genes1 <- merged_df$Description
  bg_genes <- prepare_gmt(genes1)
  # ranking 
  rankings <- sign(genes_diff$log2FoldChange)*(-log10(genes_diff$padj))  # ou juste utulise log2foldchange 
  names(rankings) <- merged_df$Description
  rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
  #plot(rankings)
  print(max(rankings))# if infiniti does not work with fgsea 
  print(min(rankings))
  if (max(rankings) == "Inf" || min(rankings) == "INF" ){
    max_ranking <- max(rankings[is.finite(rankings)])
    min_ranking <- min(rankings[is.finite(rankings)])
    rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
    rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
    }

  GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                   stats = rankings,
                   scoreType = 'std', # in this case we have both positif  and negative rankings
                   minSize = 10,
                   maxSize = 1000,
                   nproc = 1) # for parallelisation
  filename = 'data_output/fgsea_result/'
  #filename2= 'data_output/fgsea_result2/'
  #saveRDS(GSEAres, file = paste0(filename2,name,'_gsea_results.rds'))
  topPathways <- GSEAres[ES > 0][head(order(padj), n = 10), pathway]
  GSEA_filtered <- GSEAres[pathway %in% topPathways]
  GSEA_filtered <- GSEA_filtered[GSEA_filtered$padj < (0.05/32)]
  nb_patwhays <- sum(GSEA_filtered[, padj < (0.05/32)])
  print(nb_patwhays)
  GSEA_filtered <- GSEA_filtered[order(GSEA_filtered$NES , decreasing = TRUE)]
  write.table(GSEA_filtered[,c(-2,-4,-7,-8)],file= paste0("data_output/table_reactome/",name,"_conf","_fgseaRes.txt"),
              sep = "\t", row.names = FALSE)
  return(nb_patwhays)
}

GSEA2 <- function(name,transcript_count){
  # prepare bg_gene 
  genes_diff <- read.table(paste0("data_output/desq2_outputs_q3/",name,".tsv"),
                           sep = "\t", header = TRUE, stringsAsFactors = TRUE)
  genes_diff$Name <- row.names(genes_diff)
  merged_df <- merge(genes_diff, transcript_count[, c("Description", "Name")], by = "Name", all.x = FALSE)
  genes1 <- merged_df$Description
  bg_genes <- prepare_gmt(genes1)
  # ranking 
  rankings <- sign(genes_diff$log2FoldChange)*(-log10(genes_diff$padj))  # ou juste utulise log2foldchange 
  names(rankings) <- merged_df$Description
  rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
  #plot(rankings)
  print(max(rankings))# si infini error fgsea 
  print(min(rankings))
  if (max(rankings) == "Inf" || min(rankings) == "INF" ){
    max_ranking <- max(rankings[is.finite(rankings)])
    min_ranking <- min(rankings[is.finite(rankings)])
    rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
    rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
    rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
  }
  
  GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                   stats = rankings,
                   scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                   minSize = 10,
                   maxSize = 1000,
                   nproc = 1) # for parallelisation
  filename = 'data_output/fgsea_result_not_conf/'
  topPathways <- GSEAres[ES > 0][head(order(padj), n = 10), pathway]
  GSEA_filtered <- GSEAres[pathway %in% topPathways]
  GSEA_filtered <- GSEA_filtered[GSEA_filtered$padj < (0.05/32)]
  nb_patwhays <- sum(GSEA_filtered[, padj < (0.05/32)])
  print(nb_patwhays)
  GSEA_filtered <- GSEA_filtered[order(GSEA_filtered$NES , decreasing = TRUE)]

  write.table(GSEA_filtered[,c(-2,-4,-7,-8)],file= paste0("data_output/table_reactome/",name,"_fgseaRes.txt"),
              sep = "\t", row.names = FALSE)
  return(nb_patwhays)
}

features_count <- read.table("data/morphological_counts_lunit_dino.tsv", sep = "\t", header = TRUE, row.names = 1)
transcript_count <- read.table("data/RNA_read_counts.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)

cluster_name<- names(features_count)

list_path <- lapply(cluster_name, GSEA,transcript_count = transcript_count)
significant_pathways <- sapply(list_path[1:32], function(x) x[[1]])
results_df <- data.frame(Cluster = cluster_name, SignificantPathways = significant_pathways)

write.table(results_df,file= "data_output/nb_of_in_clusters/nb_pathways_Clusters_conf.txt",
            sep = "\t", row.names = FALSE)

list_path2 <- lapply(cluster_name, GSEA2,transcript_count = transcript_count)
significant_pathways2 <- sapply(list_path2[1:32], function(x) x[[1]])
results_df2 <- data.frame(Cluster = cluster_name, SignificantPathways = significant_pathways2)
write.table(results_df2,file= "data_output/nb_of_in_clusters/nb_pathways_Clusters.txt",
            sep = "\t", row.names = FALSE)





