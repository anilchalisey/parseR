#' Pipeline for DESeq2 analysis
#'
#' Function to perform entire DESeq2 analysis, beginning with
#' a SummarizedExperiment and ending with list of differential genes and
#' generation of diagnostic plots along the way.
#'
#' @param se SummarizedExperiment with an assay slot named counts_fil containing
#'   the filtered counts and colData slot containing sample metadata.
#' @param block.column The column in the samples metadata dataframe specifying
#'   the block/additive effect column.
#' @param adjust.method Method to be used for adjustment of nominal p-values.
#'   May be one of "BH", "bonferroni", "holm", "hochberg", "hommel", "BY".
#' @param p.value Value between 0 and 1.  Adjusted p-value for the differential
#'   expression analysis.
#'
#' @return
#' List containing the DESeq2 DDS container, the raw results for each contrast
#' and the processed differential gene results.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @import SummarizedExperiment
#'
#' @export

deseq2_pipe <- function(se = NULL,
                        block.column = NULL,
                        adjust.method = "BH",
                        p.value = 0.05) {

  ############
  # Create a DESeq2 object
  ############
  design.formula <- formula(paste0("~0+ ", "Contrast",
                                   ifelse(!is.null(block.column),
                                          paste0("+", block.column), "")))
  counts_fil <- SummarizedExperiment::assays(se)$counts_fil
  fil <- counts_fil[rowSums(is.na(counts_fil)) != ncol(counts_fil),]
  ddsMat <- DESeq2::DESeqDataSetFromMatrix(
    countData = fil,
    colData = SummarizedExperiment::colData(se),
    design = design.formula)
  ddsMat$Contrast <- factor(ddsMat$Contrast,
                            levels = levels(SummarizedExperiment::colData(se)$Contrast))

  ############
  # Normalise libraries
  ############
  cat("\nNormalisation Factors:\n")
  deseq2_norm(ddsMat)

  ############
  # Construct design and contrast matrices
  ############
  cat("\nDesign matrix for DESeq2 analysis:\n\n")
  design.matrix <- as.data.frame.array(model.matrix(design.formula,
                                                    data = colData(se)))
  print(design.matrix)
  cat("\nContrasts to be performed:\n\n")
  contrast <- create_contrast(design.matrix)

  ############
  # Create a PCA plot
  ############
  create_pca(ddsMat)

  ############
  # Perform differential analysis
  ############
  deseq2_results <- run_deseq2(dds = ddsMat,
                               p.value = p.value,
                               adjust.method = adjust.method)

  out <- deseq2_results[[2]]
  for (i in 1:length(out)) {
    write.table(out[[i]], row.names = F, col.names = TRUE,
                quote = F, sep = "\t",
                file = paste0("DEgenes_deseq2_", names(out)[i], ".txt"))
  }
  return(deseq2_results)
}
