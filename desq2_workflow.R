library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)

diff_exp <- function(features_count, sample_table, design) {
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

    count_depth <- tibble(sample = colnames(deseq_obj), SeqDepth = colSums(assay(deseq_obj))) %>%
        ggplot(aes(x = sample, y = SeqDepth)) +
        geom_col() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    deseq_obj <- DESeq(deseq_obj) # Running the differential analysis (normalization + stats test)

    return(deseq_obj)
}

count_depth <- function(deseq_obj) {
    tibble(sample = colnames(deseq_obj), SeqDepth = colSums(assay(deseq_obj))) %>%
        ggplot(aes(x = sample, y = SeqDepth)) +
        geom_col() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

volcano_plot <- function(deseq_obj, contrast, title = "") {
    # contrast: cf. DESeq2 documentation (results)
    res <- results(deseq_obj, contrast = contrast)
    res <- as_tibble(res, rownames = "gene")
    res %>%
        filter(!is.na(padj)) %>%
        ggplot(aes(
            x = log2FoldChange, y = -log10(padj),
            color = padj < 0.05 & abs(log2FoldChange) > 1
        )) +
        scale_colour_manual(values = c("gray", "red")) +
        geom_point(size = 0.5) +
        geom_hline(yintercept = -log10(0.05)) +
        geom_vline(xintercept = 1) +
        geom_vline(xintercept = -1) +
        ggtitle(title) +
        theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
}

significant_features_count <- function(deseq_obj, contrast, padj_th = 0.05, fold_th = 1) {
    res <- results(deseq_obj, contrast)
    res <- as_tibble(res, rownames = "gene")
    res %>%
        filter(padj < padj_th & abs(log2FoldChange) > fold_th) %>%
        nrow()
}
