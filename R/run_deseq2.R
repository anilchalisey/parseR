#' Run DESeq2 analysis
#'
#' Function to run DESeq2 analysis on a DDS object and to extract the desired
#' contrasts.
#'
#' @param dds DESeq2 dds object created by
#'            \code{\link[DESeq2]{DESeqDataSetFromMatrix}}.
#' @param p.value p-value for determining signficance for differential
#'   expression.
#' @param adjust.method method for p-value correction for multiple testing
#'
#' @return
#' A DDS object in which has undergone dispersion estimates and model fitting.
#'
#' @importFrom DESeq2 DESeq plotDispEsts results counts
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dds <- run_deseq2(dds = ddsMat)
#' }

run_deseq2 <- function(dds = dds,
                       p.value = 0.05,
                       adjust.method = "BH") {
  # run the DESEq analysis
  dds <- DESeq2::DESeq(dds, betaPrior = FALSE)
  ppi <- 600
  png(filename = "DispPlot.png",
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)

  DESeq2::plotDispEsts(dds)
  graphics.off()

  contrast.levels <- factor(unique(colData(dds)$Contrast), ordered = TRUE)
  combinations <- combn(contrast.levels, 2, simplify = FALSE)
  combinations <- lapply(combinations, as.character)
  combinations <- lapply(combinations, function(x) c("Contrast", x))
  res <- lapply(combinations, function(x) {
    DESeq2::results(dds, contrast = c(x[1], x[3], x[2]), alpha = p.value,
                    pAdjustMethod = adjust.method,
                    independentFiltering = FALSE)
  })
  res <- lapply(res, function(x) x[order(x$padj), ])

  names(res) <- paste0(lapply(combinations, function(x) {
    paste(x[3], x[2], sep = "-")
  }))


  summar <- do.call(cbind, lapply(res, function(x) {
    rbind(sum(x$log2FoldChange < 0 & x$padj < 0.05, na.rm = T),
          sum(x$padj >= 0.05, na.rm = T),
          sum(x$log2FoldChange > 0 & x$padj < 0.05, na.rm = T))
  }))
  rownames(summar) <- c(-1, 0, 1)
  colnames(summar) <- names(res)

  deseq2_results <- list()
  for (i in 1:length(res)) {
    deseq2_results[[i]] <- merge(as.data.frame(res[[i]]),
                                 (as.data.frame(DESeq2::counts(dds))),
                                 by.x = 0, by.y = 0, all.x = T)
    names(deseq2_results[[i]]) <- c("genes", "baseMean", "logFC", "lfcse",
                                    "stat", "pvalue", "FDR", colnames(dds))
    deseq2_results[[i]] <-
      deseq2_results[[i]][order(deseq2_results[[i]]$FDR), ]
  }
  names(deseq2_results) <- names(res)
  create_maplot(res, analysis = "deseq2")
  return(list(summar, deseq2_results))
}
