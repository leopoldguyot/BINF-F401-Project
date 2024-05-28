#GSEA
library('fgsea')


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
  hidden1 <- intersect(genes, hidden)
  # filter for gene sets with more than 5 genes annotated
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] 
  
  final_list <- matrix_to_list(mat)


}
# prepare background gene 
bg_genes <- prepare_gmt() # mettre gene exp diff 



# prepare ranked list of gene 


# rune GSEA 

GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation





#### test 
hidden <- unique(unlist(gmt))
names(gmt)
mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
              nrow = length(hidden), ncol = length(gmt))
for (i in 1:dim(mat)[2]){
  mat[,i] <- as.numeric(hidden %in% gmt[[i]])
}
mat2 <- mat[hidden, colnames(mat)[which(colSums(mat[hidden1,])>5)]]
final_list <- matrix_to_list(mat2)
