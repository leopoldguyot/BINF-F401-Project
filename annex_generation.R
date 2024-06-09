library(xtable)
library(stringr)

doc_head <- c("\\documentclass{article}",
              "\\usepackage{graphicx}",
              "\\usepackage{longtable}",
              "\\usepackage{booktabs}",
              "\\usepackage{multirow}",
              "\\usepackage{array}",
              "\\usepackage{float}",
              "\\usepackage{caption}",
              "\\usepackage{subcaption}",
              "\\usepackage{hyperref}",
              "\\usepackage{geometry}",
              "\\geometry{a4paper, top=2cm, bottom=2cm, left=2cm, right=2cm}",
              "\\begin{document}",
              "\\title{Annex 1: Top 10 significant results by cluster}",
              "\\maketitle",
              "\\tableofcontents",
              "\\newpage")

doc_end <- c("\\end{document}")

write(doc_head, "annex.tex", append = FALSE, sep = "\n")


# add tables

q3_1_conf <- list.files("data_output/significant_q3", pattern="^Mophological.*\\confounders.tsv$", full.names=TRUE)
q3_1 <- list.files("data_output/significant_q3", pattern="^Mophological.*\\.tsv$", full.names=TRUE)
q3_1 <- q3_1[!q3_1 %in% q3_1_conf]
q3_2_conf <- list.files("data_output/table_reactome", pattern="^Mophological.*\\conf_fgseaRes.txt$", full.names=TRUE)
q3_2 <- list.files("data_output/table_reactome", pattern="^Mophological.*\\.txt$", full.names=TRUE)
q3_2 <- q3_2[!q3_2 %in% q3_2_conf]

q3_1_conf_df <- list()

for (i in 1:length(q3_1_conf)) {
    cluster_number <- str_extract(q3_1_conf[i], "(?<=G4_)\\d+")
    q3_1_conf_df[[cluster_number]] <- read.table(q3_1_conf[i], header=TRUE, sep="\t")
}
q3_1_conf_df <- q3_1_conf_df[order(as.numeric(names(q3_1_conf_df)))]

q3_1_df <- list()

for (i in 1:length(q3_1)) {
    cluster_number <- str_extract(q3_1[i], "(?<=G4_)\\d+")
    q3_1_df[[cluster_number]] <- read.table(q3_1[i], header=TRUE, sep="\t")
}
q3_1_df <- q3_1_df[order(as.numeric(names(q3_1_df)))]

q3_2_conf_df <- list()

for (i in 1:length(q3_2_conf)) {
    cluster_number <- str_extract(q3_2_conf[i], "(?<=G4_)\\d+")
    q3_2_conf_df[[cluster_number]] <- read.table(q3_2_conf[i], header=TRUE, sep="\t")
}
q3_2_conf_df <- q3_2_conf_df[order(as.numeric(names(q3_2_conf_df)))]

q3_2_df <- list()

for (i in 1:length(q3_2)) {
    cluster_number <- str_extract(q3_2[i], "(?<=G4_)\\d+")
    q3_2_df[[cluster_number]] <- read.table(q3_2[i], header=TRUE, sep="\t")
}
q3_2_df <- q3_2_df[order(as.numeric(names(q3_2_df)))]
# write tables

write("\\section{Q3.1 Top 10 up-regulated transcripts without adjustement for confondant effect}", "annex.tex", append=TRUE, sep="\n")
for (i in 1:length(q3_1_df)) {
    write(paste("\\subsection{Cluster", names(q3_1_df)[i], "}"), "annex.tex", append=TRUE, sep="\n")
    if (nrow(q3_1_df[[i]]) == 0) {
        write("No significant transcript", "annex.tex", append=TRUE, sep="\n")
        next
    }
    latex_code <- capture.output(print(xtable(q3_1_df[[i]], caption=paste("Top 10 up-regulated transcripts for cluster", names(q3_1_df)[i]), label=paste0("tab:q3_1_", names(q3_1_df)[i]), digits=4), type = "latex", table.placement = "H"))

    write(latex_code, "annex.tex", append=TRUE, sep="\n")
}

write("\\section{Q3.1 Top 10 up-regulated transcript with adjustement for confondant effect}", "annex.tex", append=TRUE, sep="\n")
for (i in 1:length(q3_1_conf_df)) {
    write(paste("\\subsection{Cluster", names(q3_1_conf_df)[i], "}"), "annex.tex", append=TRUE, sep="\n")
    if (nrow(q3_1_conf_df[[i]]) == 0) {
        write("No significant transcript", "annex.tex", append=TRUE, sep="\n")
        next
    }
    latex_code <- capture.output(print(xtable(q3_1_conf_df[[i]], caption=paste("Top 10 up-regulated transcripts for cluster", names(q3_1_conf_df)[i], "with adjustement for confondant effect"), label=paste0("tab:q3_1_conf_", names(q3_1_conf_df)[i]), digits=4), type = "latex", table.placement = "H"))

    write(latex_code, "annex.tex", append=TRUE, sep="\n")
}

write("\\section{Q3.2 Top 10 up-regulated pathways without adjustement for confondant effect}", "annex.tex", append=TRUE, sep="\n")
for (i in 1:length(q3_2_df)) {
    write(paste("\\subsection{Cluster", names(q3_2_df)[i], "}"), "annex.tex", append=TRUE, sep="\n")
    if (nrow(q3_2_df[[i]]) == 0) {
        write("No significant pathway", "annex.tex", append=TRUE, sep="\n")
        next
    }
    table <- q3_2_df[[i]][, -3]
    table$pathway <- sapply(table$pathway, function(x) gsub("_", " ", x))
    latex_code <- capture.output(print(xtable(table, caption=paste("Top 10 up-regulated pathways for cluster", names(q3_2_df)[i]), label=paste0("tab:q3_2_", names(q3_2_df)[i]), digits=4, align = c("p{0.05\\linewidth}","p{0.7\\linewidth}", "p{0.1\\linewidth}", "p{0.1\\linewidth}")), type = "latex", table.placement = "H"))

    write(latex_code, "annex.tex", append=TRUE, sep="\n")
}

write("\\section{Q3.2 Top 10 up-regulated pathways with adjustement for confondant effect}", "annex.tex", append=TRUE, sep="\n")
for (i in 1:length(q3_2_conf_df)) {
    write(paste("\\subsection{Cluster", names(q3_2_conf_df)[i], "}"), "annex.tex", append=TRUE, sep="\n")
    if (nrow(q3_2_conf_df[[i]]) == 0) {
        write("No significant pathway", "annex.tex", append=TRUE, sep="\n")
        next
    }
    table <- q3_2_conf_df[[i]][, -3]
    table$pathway <- sapply(table$pathway, function(x) gsub("_", " ", x))
    latex_code <- capture.output(print(xtable(table, caption=paste("Top 10 up-regulated pathways for cluster", names(q3_2_conf_df)[i], "with adjustement for confondant effect"), label=paste0("tab:q3_2_conf_", names(q3_2_conf_df)[i]), digits=4, align = c("p{0.05\\linewidth}","p{0.7\\linewidth}", "p{0.1\\linewidth}", "p{0.1\\linewidth}")), type = "latex", table.placement = "H"))

    write(latex_code, "annex.tex", append=TRUE, sep="\n")
}

write(doc_end, "annex.tex", append=TRUE, sep="\n")
