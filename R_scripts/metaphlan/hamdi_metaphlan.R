library(tidyverse)
library(limma)
library(ALDEx2)
source("helper_funcs.R")

setwd("C:/Users/Michael/Desktop/Code/R/Scripts/HumanN")

sample_IDs <- c("HF1", "HF2", "HF3", "HF4", "HF5", "HF6",
                "LF1", "LF2", "LF3", "LF4", "LF5", "LF6")
sample_groups <- c("HF", "HF", "HF", "HF", "HF", "HF",
                   "LF", "LF", "LF", "LF", "LF", "LF")
output_folder <- "outputs/hamdi/outputs_metaphlan"


# Get folders ready
dir.create(output_folder)

# Import and merge the data. Also clean it up a bit
metaphlan_table <- read_tsv("metaphlan_out_hamdi.tsv")
saveRDS(metaphlan_table, file = paste0(output_folder, "/original_merged_table.RDS"))

# Rename to two groups
metaphlan_twogroups <- metaphlan_table
metaphlan_twogroups$SampleID <- NULL
colnames(metaphlan_twogroups) <- sample_groups

# Limma DE analysis
fit2 <- DE_analysis_HFvLF(metaphlan_twogroups)

# Add some columns back to genefam
metaphlan_table$Pvalue <- as.vector(fit2$p.value)
metaphlan_table$Log_odds_of_DE <- as.vector(fit2$lods)
write_tsv(metaphlan_table, file = paste0(output_folder, "/intermediate.tsv"))

# Remove rows with NAs in Pvalue column
genefam_df_no_NA <- drop_na(metaphlan_table, "Pvalue")

# Add gene names to fit2
fit2$genes <- metaphlan_table$SampleID

# Return table of geneID, logFC, AvgExp, pvalue, adj.pvalue, ordered according to adj.pvalue
results <- topTable(fit2, number = nrow(metaphlan_table))
rownames(results) <- NULL
write_tsv(results, file = paste0(output_folder, "/log_fold_change.tsv"))

# Remove all rows where P.Value is NA in results
results_semiprocessed <- drop_na(results, "P.Value")
number_of_rows_removed_NA <- nrow(results) - nrow(results_semiprocessed)
message(paste0("Removed ", number_of_rows_removed_NA, " rows that had an NA in the p-value, ",
              "out of a total of ", nrow(results), " rows"))

# Remove all rows where the ID includes "unclassified"
# results_semiprocessed <- results_no_NA[!grepl("unclassified", results_no_NA$ID),]
# number_of_rows_removed_unclassified <- nrow(results_no_NA) - nrow(results_semiprocessed)
# message(paste0("Removed ", number_of_rows_removed_unclassified, " rows that had 'unclassified' in the ID, ",
#               "out of a total of ", nrow(results_no_NA), " rows"))

# Report processing stats
lines_to_write <- c(paste0("Original number of rows: ", nrow(metaphlan_table)),
                    paste0("Rows remaining after removing NAs: ", nrow(results_semiprocessed)))
fileConn <- file(paste0(output_folder, "/processing_stats.txt"))
writeLines(lines_to_write, fileConn)
close(fileConn)
write_tsv(results_semiprocessed, file = paste0(output_folder, "/results_semiprocessed.tsv"))

# ggplot2 step
# HF if adj -log(adj.p) > -log(0.05) AND logFC < -2
# LF if adj -log(adj.p) > -log(0.05) AND logFC > 2
labelled_results <- ID_logFC_adjP(results_semiprocessed, 'HF', 'LF')

volcano_plot <- ggplot(data = labelled_results,
                       aes_string(x = 'logFC', y = 'logadjP', color = 'volcano_ID')) +
        geom_point(size = 3) +
        labs(x = 'log Fold Change',
             y = "-log10(adjusted P-value)",
             color = element_blank()) +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        geom_vline(xintercept = c(-2, 2), linetype = 'dashed') +
        scale_color_manual(values = c('#8AB4F8', '#acb4bc'))

ggsave(volcano_plot,
       filename = 'genefam_volcano_plot_limma.png',
       dpi = 600,
       height = 6,
       width = 8,
       device = "png",
       path = output_folder)

# Count the number of HFs, LFs, and NSs in labelled_results, then report
ID_numbers <- table(labelled_results$volcano_ID)
num_DE_LF <- ID_numbers[1]
message(paste0("Of a total of ", sum(ID_numbers), " taxa, ",
               num_DE_LF, " were differentially abundant in LF mice, and ",
               0, " were differentially abundant in HF mice."))

# Save a table with only DE taxa
only_DE <- labelled_results[!(labelled_results[["volcano_ID"]] == "NS"), ]
write_tsv(only_DE, file = paste0(output_folder, "/DE_taxa_limma.tsv"))

# Do another DE analysis, but this time with ALDEx2
aldex_output <- aldex_taxa_analysis(metaphlan_table, sample_groups, "zero",
                                    lowES_cutoff = -1, highES_cutoff =  1,
                                    negES_label = "HF", posES_label = "LF")
write.table(aldex_output, file = paste0(output_folder, "/aldex_raw_output.tsv"),
            quote = FALSE, sep = '\t',
            row.names = FALSE, col.names = TRUE)

aldex_volcano <- ggplot(aldex_output, aes_string(x = 'effect', y = 'logadj_wi.eBH', color = "point_ID")) +
        geom_point() +
        labs(x = 'Effect Size',
             y = "-log10(adjusted P-value)",
             color = element_blank()) +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
        scale_color_manual(values = c('#80cb9f', '#acb4bc'))
ggsave(aldex_volcano,
       filename = 'genefam_volcano_plot_aldex.png',
       dpi = 600,
       height = 6,
       width = 8,
       device = "png",
       path = output_folder)
only_DE_aldex <- aldex_output[!(aldex_output[["point_ID"]] == "NS"), ]
write_tsv(only_DE_aldex, file = paste0(output_folder, "/DE_taxa_aldex.tsv"))

# Make a PCoA of all samples, separated by pathways
table_for_PCoA <- dplyr::select(metaphlan_table, -SampleID, -Pvalue, -Log_odds_of_DE)
transposition <- as.data.frame(t(table_for_PCoA))
transposition$labels <- sample_groups
plot <- make_PCoA(transposition, 'labels')
PCoA_plot <- plot +
        geom_point(size = 3) +
        ggtitle('PCoA based on taxa abundance') +
        labs(color = element_blank())
ggsave(filename = 'taxa_PCoA_plot.png',
       dpi = 600,
       height = 6,
       width = 8,
       device = "png",
       path = output_folder)
