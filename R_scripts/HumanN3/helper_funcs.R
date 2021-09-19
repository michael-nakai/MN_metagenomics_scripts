library(tidyverse)
library(docstring)
library(limma)
library(ggplot2)
library(mgcv)
library(vegan)
library(lattice)
library(labdsv)

merge_humanN_outputs <- function(type, data_list) {

    #' Merges multiple tsv files outputted by HumanN2 into
    #' one dataframe, which is returned. Works with pathabundance,
    #' pathcoverage, and genefamilies.
    #'
    #' @param type The input data type. Can be 'pathabun', 'pathcov' or 'genefam'.
    #' @param data_list The path to the tsv files. Can specify as many files as needed here.

    col1 <- switch(type, "pathabun" = "# Pathway", "pathcov" = "# Pathway", "genefam" = "# Gene Family")

    i <- 1
    for (dataframe in data_list) {
        if (i == 1) {
            intermediate <- dataframe
            i <- i + 1
        } else {
            intermediate <- merge(x = intermediate, y = dataframe, by = col1, all = TRUE)
        }
    }

    intermediate[is.na(intermediate)] <- 0
    return(intermediate)

}

generate_dataframe_list <- function(path_to_folder) {

    #' Returns a list with all tsv files imported
    #' from the folder provided.
    #'
    #' @param path_to_folder The filepath to the folder, as a string

    filelist <- list.files(path = path_to_folder, full.names = TRUE)
    dflist <- vector(mode = "list", length = length(filelist))

    i <- 1
    for (filepath in filelist) {
        dflist[[i]] <- read_tsv(filepath)
        i <- i + 1
    }

    return(dflist)
}

DE_analysis_HFvLF <- function(data_frame) {

    #' Performs lmfit, contrasts.fit, and eBayes on data_frame.
    #' Note that values in the data_frame are log transformed prior
    #' to fitting. Returns the output of eBayes.
    #'
    #' @param data_frame A dataframe containing the columns to be compared.

    TS <- as.factor(colnames(data_frame))
    design <- model.matrix(~0+TS)
    colnames(design) <- levels(TS)
    fit <- lmFit(log(data_frame), design)
    cont.matrix <- makeContrasts(HF_vs_LF = LF-HF, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, trend=F)

    return(fit2)
}

ID_logFC_adjP <- function(data_frame, negFC_label, posFC_label) {

    #' Takes a dataframe with the columns 'adj.P.Val' and 'logFC', and
    #' labels them accordingly.
    #'
    #' @param data_frame The dataframe containing the necessary columns
    #' @param negFC_label The label for samples with a logFC < -2 and adj.P.Val < 0.05
    #' @param posFC_label The label for samples with a logFC > 2 and adj.P.Val < 0.05

    data_frame[['volcano_ID']] <- rep(NA, nrow(data_frame))
    data_frame[['logadjP']] <- rep(NA, nrow(data_frame))
    rownum <- 1
    for (pval in data_frame[['adj.P.Val']]) {
        FC <- data_frame[['logFC']][rownum]
        if (-log10(pval) > -log10(0.05) & FC < -2) {
            data_frame[['volcano_ID']][rownum] <- negFC_label
        } else if (-log10(pval) > -log10(0.05) & FC > 2) {
            data_frame[['volcano_ID']][rownum] <- posFC_label
        } else {
            data_frame[['volcano_ID']][rownum] <- 'NS'
        }
        data_frame[['logadjP']][rownum] <- -log10(pval)
        rownum <- rownum + 1
    }

    return(data_frame)
}

create_volcano_plot <- function(data_frame, HF_color, LF_color, NS_color) {

    #' Used to create a volcano plot from the output of ID_logFC_adjP
    #'
    #' @param data_frame The dataframe outputted by ID_logFC_adjP
    #' @param HF_color The hexcode for the HF group's color
    #' @param LF_color The hexcode for the LF group's color
    #' @param NS_color The hexcode for the non-significant points' color

    ggplot(data = data_frame,
           aes_string(x = 'logFC', y = 'logadjP', color = 'volcano_ID')) +
        geom_point() +
        labs(x = 'log Fold Change',
             y = "-log10(adjusted P-value)",
             color = element_blank()) +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        geom_vline(xintercept = c(-2, 2), linetype = 'dashed') +
        scale_color_manual(values = c(HF_color, LF_color, NS_color)) # HF, LF, NS
}

create_only_DE_gene_tsv <- function(data_frame) {

    #' Used to create a table containing only DE genefams from the
    #' output of ID_logFC_adjP
    #'
    #' @param data_frame The outputted dataframe from ID_logFC_adjP

    only_DE <- data_frame[!(data_frame[["volcano_ID"]] == "NS"), ]
    only_DE[['Name']] <- sub("UniRef90_.*: ", "", only_DE$ID)
    only_DE[['UniRefID']] <- sub(":.*", "", sub("UniRef90_", "", only_DE$ID))
    write_tsv(only_DE, file = paste0(output_folder, "/Differentially_abundant_genefamilies.tsv"))

    return(only_DE)
}

make_PCoA <- function(your_dataframe, label_colname) {

    #' Creates and returns a ggplot2 PCoA out of a provided dataframe. The dataframe
    #' should only include columns with taxa/pathway/gene abundances and a
    #' column containing the group that each sample belongs to. Samples should
    #' be on rows, and taxa/pathways/genes should be on columns, so you might
    #' have to transpose the dataframe first.
    #'
    #' @param data_frame The dataframe containing the abundances and the label column
    #' @param label_colname The column name of the label column

    data_frame <- your_dataframe
    labels <- as.vector(data_frame[[label_colname]])
    data_frame[[label_colname]] <- NULL
    i <- c(1:ncol(data_frame))
    data_frame[ , i] <- apply(data_frame[ , i], 2, function(x) as.numeric(as.character(x)))

    x <- scale(data_frame)
    E <- cov(x)
    evd <- eigen(E)

    ev1 <- evd$vectors[,1]
    ev2 <- evd$vectors[,2]
    pc1 <- x %*% ev1
    pc2 <- x %*% ev2
    pc1.v <- evd$values[1]
    pc2.v <- evd$values[2]
    evd_total <- sum(evd$values)

    df <- data.frame(cbind(pc1, pc2), labels)
    names(df) <- c('PC1', 'PC2', label_colname)
    pc1_var <- pc1.v / evd_total
    pc2_var <- pc2.v / evd_total
    PCoA_plot <- ggplot(df,aes_string(x = 'PC1', y = 'PC2', colour = label_colname)) + geom_point(alpha = 1) +
        xlab(paste0("PC1 (", round(pc1_var*100, 1), "%)")) +
        ylab(paste0("PC2 (", round(pc2_var*100, 1), "%)")) +
        scale_color_manual(values = c('#80cb9f', '#8AB4F8'))
    return(PCoA_plot)
}

aldex_taxa_analysis <- function(data_frame, sample_groups, denom,
                                lowES_cutoff, highES_cutoff,
                                negES_label, posES_label) {

    #' Performs the aldex function on the dataframe, with some formatting.
    #' Also ID's the columns with effect size > highEF_cutoff or
    #' < lowEF_cutoff and BH-adjP < 0.05.
    #' Returns an aldex2-outputted dataframe.
    #'
    #' @param data_frame The dataframe in question
    #' @param sample_groups The sample groups for each column
    #' @param denom The denom type to use for the aldex function
    #' @param lowES_cutoff The negative effect size cutoff
    #' @param highES_cutoff The positive effect size cutoff
    #' @param negES_label The label to use for significant negative ES samples
    #' @param posES_label The label to use for significant positive ES samples

    newtable <- as.data.frame(data_frame)
    newtable <- dplyr::select(newtable, -Pvalue, -Log_odds_of_DE, -SampleID)
    newtable <- newtable * 10^7 #Aldex needs each col to be an int, not a dbl
    newtable <- as.data.frame(sapply(newtable, as.integer))
    aldex_output <- aldex(newtable, conditions = sample_groups, denom = denom)
    aldex_output$point_ID <- rep(NA, nrow(aldex_output))
    aldex_output$logadj_wi.eBH <- rep(NA, nrow(aldex_output))

    rownum <- 1
    for (pval in aldex_output$wi.eBH) {
        effect_size <- aldex_output[['effect']][rownum]
        if (pval < 0.05 & effect_size < lowES_cutoff) {
            aldex_output[['point_ID']][rownum] <- negES_label
        } else if (pval < 0.05 & effect_size > highES_cutoff) {
            aldex_output[['point_ID']][rownum] <- posES_label
        } else {
            aldex_output[['point_ID']][rownum] <- 'NS'
        }
        aldex_output[['logadj_wi.eBH']][rownum] <- -log10(pval)
        rownum <- rownum + 1
    }

    aldex_output$Taxa_name <- data_frame[['SampleID']]

    return(aldex_output)
}
