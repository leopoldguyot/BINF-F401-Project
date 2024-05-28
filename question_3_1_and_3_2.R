library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

diff_exp <- function(features_count, sample_table, design, name) {
    # features_count: a matrix with the features as rows and the samples as columns
    # sample_table: a DataFrame or data.frame with at least a single column. Rows of colData correspond to columns of countData
    # design: a formula object specifying the design

    stopifnot(is.matrix(features_count))
    stopifnot(is.data.frame(sample_table))
    stopifnot(is(design, "formula"))
    deseq_obj <- DESeqDataSetFromMatrix(
        countData = features_count,
        colData = sample_table,
        design = design
    )

    deseq_obj <- deseq_obj[rowSums(assay(deseq_obj)) > 0, ] # removing rows containing only zeros

    deseq_obj <- DESeq(deseq_obj) # Running the differential analysis (normalization + stats test)

    result_table <- results(deseq_obj, name = resultsNames(deseq_obj)[2])
    filename <- paste("data_output/desq2_outputs_q3/", name, ".tsv", sep = "")
    print(filename)

    write.table(result_table, file = filename, sep = "\t", row.names = TRUE)
    return(NULL)
}

extract_count_significant <- function(result_table, padj_th = 0.05, fold_th = 1) {
    res <- as_tibble(result_table, rownames = "gene")
    return(res %>%
        filter(padj < padj_th & abs(log2FoldChange) > fold_th) %>%
        nrow())
}

extract_count_wrapper <- function(result_folder, padj_th = 0.05, fold_th = 1) {
    result_files <- list.files(result_folder, pattern = "*.tsv", full.names = TRUE)
    significant_results <- list()
    for (file in result_files) {
        result_table <- read.table(file, sep = "\t", header = TRUE)
        name <- gsub(".*\\/(.*)\\.tsv", "\\1", file)
        significant_results[[name]] <- extract_count_significant(result_table)
    }
    return(significant_results)
}

report_up_regulated <- function(result_table, padj_th = 0.05, fold_th = 1, name) {
    res <- as_tibble(result_table, rownames = "gene")
    res <- res %>%
        filter(padj < padj_th & log2FoldChange > fold_th) %>%
        arrange(desc(log2FoldChange)) %>%
        slice_head(n = 10) %>%
        select(gene, log2FoldChange, padj)
    
    write.table(res, paste0("data_output/significant_q3/", name, ".tsv"), sep = "\t", row.names = FALSE)
}


clinical_data <- read.table("data_output/clinical_data_preprocessed.tsv", sep = "\t", header = TRUE)
cluster_data <- read.table("data/morphological_counts_lunit_dino.tsv", sep = "\t", header = TRUE)
cluster_data[-1] <- log10(cluster_data[-1] + 1)
sample_table <- merge(clinical_data, cluster_data, by = "SMPLID")
rownames(sample_table) <- sample_table$SMPLID


technical_confounders <- "COHORT + DTHHRDY_3 + TRISCHD + DTHVNT"
designs <- colnames(cluster_data)[-1] %>% map(~as.formula(paste("~", .x)))
names(designs) <- colnames(cluster_data)[-1]
designs_confounders <- colnames(cluster_data)[-1] %>% map(~as.formula(paste("~", .x, "+", technical_confounders)))
names(designs_confounders) <- paste(colnames(cluster_data)[-1], "_confounders", sep = "")
full_designs <- c(designs, designs_confounders)

rna_counts <- read.table("data/RNA_read_counts_filtered.tsv", sep = "\t", header = TRUE)
rownames(rna_counts) <- rna_counts$Name
rna_counts <- as.matrix(rna_counts[-c(1:2, ncol(rna_counts))])
lapply(names(full_designs), function(x) diff_exp(rna_counts, sample_table, full_designs[[x]], x))

# results exploration
bonferonni <- 0.05/ncol(cluster_data[-1])

significant_count <- extract_count_wrapper("data_output/desq2_outputs_q3/", padj_th = bonferonni, fold_th = 1)
non_confounders_elements <- significant_count[grep("confounders", names(significant_count), invert = TRUE)]
confounders_elements <- significant_count[grep("confounders", names(significant_count))]


significant_count_df <- data.frame(cluster = names(non_confounders_elements),
 count_non_confounders = unlist(non_confounders_elements),
  count_confounders = unlist(confounders_elements))
write.table(significant_count_df, "data_output/significant_q3/significant_count.tsv", sep = "\t", row.names = FALSE)

lapply(list.files("data_output/desq2_outputs_q3/", pattern = "*.tsv", full.names = TRUE), function(x) {
    result_table <- read.table(x, sep = "\t", header = TRUE)
    name <- gsub(".*\\/(.*)\\.tsv", "\\1", x)
    report_up_regulated(result_table, padj_th = bonferonni, fold_th = 1, name)
})
