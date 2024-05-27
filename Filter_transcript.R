
transcript_count <- read.table("data/RNA_read_counts.tsv", sep = "\t", header = TRUE, stringsAsFactors = TRUE)


sum(duplicated(transcript_count))
which(is.na(transcript_count))
nrow(transcript_count)

#### filter using mad()  ###
mad_values <- apply(transcript_count[,c(-1,-2)], 1, mad) # apply mad on each row 
transcript_count$mad_values <- mad_values

# sort dataframe usnig mad values
transcript_count <- transcript_count[order(transcript_count$mad_values, decreasing = TRUE),]

# pour retirer la moitier 
num_rows <- nrow(transcript_count)/2

transcript_count_filtered <- transcript_count[1:num_rows,] # reste tj des mad value =  0 

transcript_count_filtered <- transcript_count[transcript_count$mad_values > 0, ]
nrow(transcript_count_filtered)
transcript_count_filtered <- transcript_count_filtered[rowSums(transcript_count_filtered[,c(-1,-2)]) >= 500 , ]
nrow(transcript_count_filtered)

write.table(transcript_count_filtered, "data/RNA_read_counts_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
