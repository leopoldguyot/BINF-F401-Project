
Gene_count <- read.table("data/RNA_read_counts.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)


sum(duplicated(Gene_count))
which(is.na(Gene_count))
nrow(Gene_count)
### remove genes that do not have a total count over 10 ### 
Gene_count_new <- Gene_count[rowSums(Gene_count[,c(-1,-2)]) >= 10 , ]
Gene_count_new2 <- Gene_count[rowSums(Gene_count[,c(-1,-2)]) >= 500 , ]




#### filter using mad()  ###
mad_values <- apply(Gene_count[,c(-1,-2)], 1, mad) # apply mad on each row 
Gene_count$mad_values <- mad_values

# sort dataframe usnig mad value 
Gene_count <- Gene_count[order(Gene_count$mad_values, decreasing = TRUE),]

# pour retirer la moitier 
num_rows <- nrow(Gene_count)/2

Gene_count_mad_filtering <- Gene_count[1:num_rows,] # reste tj des mad value =  0 

Gene_count_mad_filtering2 <- Gene_count[Gene_count$mad_values > 0, ]

                    
                    