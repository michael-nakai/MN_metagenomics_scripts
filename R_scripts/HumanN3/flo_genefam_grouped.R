library(tidyverse)
library(limma)
source("helper_funcs.R")

setwd("C:/Users/Michael/Desktop/Code/R/Scripts/HumanN")

sample_IDs <- c("HF1", "HF2", "HF3",
                "LF1", "LF2", "LF3")
sample_groups <- c("HF", "HF", "HF",
                   "LF", "LF", "LF")
folder_with_tsvs <- "flo_agg/genefam_grouped"
output_folder <- "outputs/flo/outputs_genefam_grouped"


# Get folders ready
dir.create(output_folder)

# Import and merge the data. Also clean it up a bit
dflist <- generate_dataframe_list(folder_with_tsvs)
genefam_df <- merge_humanN_outputs("genefam", dflist)
colnames(genefam_df) <- c("Gene_Family", sample_IDs)
saveRDS(genefam_df, file = paste0(output_folder, "/original_merged_table.RDS"))

# Rename to two groups
genefam_twogroups <- genefam_df
genefam_twogroups$Gene_Family <- NULL
colnames(genefam_twogroups) <- sample_groups

# Limma DE analysis
fit2 <- DE_analysis_HFvLF(genefam_twogroups)

# Add some columns back to genefam
genefam_df$Pvalue <- as.vector(fit2$p.value)
genefam_df$Log_odds_of_DE <- as.vector(fit2$lods)
write_tsv(genefam_df, file = paste0(output_folder, "/intermediate.tsv"))

# Remove rows with NAs in Pvalue column
genefam_df_no_NA <- drop_na(genefam_df, "Pvalue")

# Add gene names to fit2
fit2$genes <- genefam_df$Gene_Family

# Return table of geneID, logFC, AvgExp, pvalue, adj.pvalue, ordered according to adj.pvalue
results <- topTable(fit2, number = nrow(genefam_df))
rownames(results) <- NULL
write_tsv(results, file = paste0(output_folder, "/log_fold_change.tsv"))

# Remove all rows where P.Value is NA in results
results_semiprocessed <- drop_na(results, "P.Value")
number_of_rows_removed_NA <- nrow(results) - nrow(results_no_NA)
message(paste0("Removed ", number_of_rows_removed_NA, " rows that had an NA in the p-value, ",
              "out of a total of ", nrow(results), " rows"))

# Remove all rows where the ID includes "unclassified"
# results_semiprocessed <- results_no_NA[!grepl("unclassified", results_no_NA$ID),]
# number_of_rows_removed_unclassified <- nrow(results_no_NA) - nrow(results_semiprocessed)
# message(paste0("Removed ", number_of_rows_removed_unclassified, " rows that had 'unclassified' in the ID, ",
#               "out of a total of ", nrow(results_no_NA), " rows"))

# Report processing stats
lines_to_write <- c(paste0("Original number of rows: ", nrow(genefam_df)),
                    paste0("Rows remaining after removing NAs: ", nrow(results_semiprocessed)))
fileConn<-file(paste0(output_folder, "/processing_stats.txt"))
writeLines(lines_to_write, fileConn)
close(fileConn)
write_tsv(results_semiprocessed, file = paste0(output_folder, "/results_semiprocessed.tsv"))

# ggplot2 step
# HF if adj -log(adj.p) > -log(0.05) AND logFC < -2
# LF if adj -log(adj.p) > -log(0.05) AND logFC > 2
labelled_results <- ID_logFC_adjP(results_semiprocessed, 'HF', 'LF')

volcano_plot <- create_volcano_plot(labelled_results, '#80cb9f', '#8AB4F8', '#acb4bc')
ggsave(filename = 'genefam_volcano_plot.png',
       dpi = 600,
       height = 6,
       width = 8,
       device = "png",
       path = output_folder)

# Count the number of HFs, LFs, and NSs in labelled_results, then report
ID_numbers <- table(labelled_results$volcano_ID)
num_DE_HF <- ID_numbers[1]
num_DE_LF <- ID_numbers[2]
message(paste0("Of a total of ", sum(ID_numbers), " gene families, ",
               num_DE_LF, " were differentially abundant in LF mice, and ",
               num_DE_HF, " were differentially abundant in HF mice."))

# Save a table with only DE genes
only_DE <- labelled_results[!(labelled_results[["volcano_ID"]] == "NS"), ]
write_tsv(only_DE, file = paste0(output_folder, "/DE_genefam_ECgrouping.tsv"))

# Make a PCoA of all samples, separated by pathways
table_for_PCoA <- dplyr::select(genefam_df, -Gene_Family, -Pvalue, -Log_odds_of_DE)
transposition <- as.data.frame(t(table_for_PCoA))
transposition$labels <- sample_groups
plot <- make_PCoA(transposition, 'labels')
saveRDS(plot, file = paste0(output_folder, "/base_PCoA_plot.RDS"))
PCoA_plot <- plot +
        ggtitle('PCoA based on gene family abundance, grouped by level 4 EC number') +
        geom_point(size = 3) +
        labs(color = element_blank())
ggsave(filename = 'taxa_PCoA_plot.png',
       dpi = 600,
       height = 6,
       width = 8,
       device = "png",
       path = output_folder)
