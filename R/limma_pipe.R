#' Pipeline for limma+voom analysis
#'
#' Function to perform entire limma+voom analysis using voom with quality
#' weights, and beginning with a SummarizedExperiment and ending with list of
#' differential genes and generation of diagnostic plots along the way.
#'
#' @param se SummarizedExperiment with an assay slot named counts_fil containing
#'   the filtered counts and colData slot containing sample metadata.
#' @param design A design matrix specifying the experimental design.  If NULL,
#'   then a design matrix will be created using the values for contrast.levels
#'   and block.levels.
#' @param contrast A matrix of contrasts (as created by \code{makeContrasts} in
#'   the \code{limma} package specifying the comparisons to perform).  If NULL,
#'   then will be created from the design matrix.
#' @param block.column The column in the samples metadata dataframe specifying
#'   the block/additive effect column.  Ignored if design is not NULL.
#' @param norm.method Normalisation method to use - may be "TMM", "RLE",
#'   "upperquartile" or "none".
#' @param adjust.method Method to be used for adjustment of nominal p-values.
#'   May be one of "BH", "bonferroni", "holm", "hochberg", "hommel", "BY".
#' @param p.value Value between 0 and 1.  Adjusted p-value for the differential
#'   expression analysis.
#'
#' @return
#' Results of limma+voom differential analysis.
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom edgeR DGEList
#'
#' @export

limma_pipe <- function(se = NULL,
                       design = NULL,
                       contrast = NULL,
                       block.column = NULL,
                       norm.method = "TMM",
                       adjust.method = "BH",
                       p.value = 0.05) {

  ############
  # Create a dge object
  ############
  counts <- SummarizedExperiment::assays(se)$counts_fil
  counts <- counts[rowSums(is.na(counts)) != ncol(counts),]
  dge <- edgeR::DGEList(counts = counts,
                        genes = rownames(counts),
                        group = SummarizedExperiment::colData(se)$Contrast)

  ############
  # Normalise libraries
  ############
  dge <- edger_norm(dge = dge, norm.method = norm.method)

  ############
  # Create a MDS plot
  ############
  create_mds(dge)

  ############
  # Create experimental design
  ############

  if (is.null(design)) {
    # Create design and contrast matrix.
    cat("\nDesign matrix:\n")
    design.matrix <- create_design(se = se, block.column = block.column)
    cat("\n\nContrast options:\n")
    contrast.matrix <- create_contrast(design = design.matrix)
  } else {
    cat("\nDesign matrix:\n")
    design.matrix <- design
    cat("\n\nContrast options:\n")
    contrast.matrix <- contrast
  }

  ############
  # Perform differential analysis
  ############
  limma_results <- run_limma(dge = dge, design = design.matrix,
                             contrast = contrast.matrix,
                             adjust.method = adjust.method,
                             p.value = p.value)

  out <- limma_results[[2]]
  for (i in 1:length(out)) {
    write.table(out[[i]], row.names = F, col.names = TRUE,
                quote = F, sep = "\t",
                file = paste0("DEgenes_limma_", names(out)[i], ".txt"))
  }

  return(limma_results)
}
