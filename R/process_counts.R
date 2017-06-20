#' Pre-processing of counts
#'
#' @description  Function to proces raw counts before differential analysis.
#' Removes counts below the counts-per-million threshold set by the user and
#' generates diagnostic plots showing the count distribution before and after.
#'
#' @param se SummarizedExperiment object containing count matrix.
#' @param cpm.cutoff The counts-per-million threshold for filtering counts.
#'
#' @return
#' SummarizedExperiment object containing the raw and filtered counts.
#'
#' @export
#' @import ggplot2
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata
#'
#' @examples
#' \dontrun{
#' filtered_se <- process_counts(se = se, cpm.cutoff = 1)
#' }

process_counts <- function(se = NULL,
                           cpm.cutoff = 1) {

  counts <- SummarizedExperiment::assay(se)
  samples <- SummarizedExperiment::colData(se)
  cpm <- edgeR::cpm(counts)
  lcpm <- edgeR::cpm(counts, log = T)
  minSampNumber <- min(table(samples$Contrast))
  keep <- rowSums(cpm > cpm.cutoff) >= minSampNumber
  counts_fil <- counts
  counts_fil[!keep, ] <- NA

  lcpm <- suppressMessages(reshape2::melt(lcpm))
  colnames(lcpm) <- c("GeneID", "Sample", "Count")
  lcpm$Sample <- paste0(lcpm$Sample, "     ")

  lcp <-
    ggplot(lcpm, aes(x = Count, group = Sample,
                     fill = Sample, colour = Sample)) +
    geom_density(alpha = 0.3, size = 1) + .theme_Publication() +
    xlab(expression(log[2]*~cpm)) + ggtitle("A. Raw data") +
    theme(legend.key.size = unit(0.7, 'lines'),
          legend.direction = "horizontal",
          legend.title = element_blank()) +
    scale_fill_manual(values = col.palette) +
    scale_colour_manual(values = col.palette) +
    guides(fill = guide_legend(ncol = min(ncol(cpm), 5)))

  fil <- counts_fil[rowSums(is.na(counts_fil)) != ncol(counts_fil),]
  lcpm_fil <- edgeR::cpm(fil, log = T)
  lcpm_fil <- suppressMessages(reshape2::melt(lcpm_fil))
  colnames(lcpm_fil) <- c("GeneID", "Sample", "Count")
  lcpm_fil$Sample <- paste0(lcpm_fil$Sample, "     ")

  lcp_fil <-
    ggplot(lcpm_fil, aes(x = Count, group = Sample,
                     fill = Sample, colour = Sample)) +
    geom_density(alpha = 0.2, size = 1) + .theme_Publication() +
    xlab(expression(log[2]*~cpm)) + ggtitle("B. Filtered data") +
    theme(legend.key.size = unit(0.7, 'lines'),
          legend.direction = "horizontal",
          legend.title = element_blank()) +
    scale_fill_manual(values = col.palette) +
    scale_colour_manual(values = col.palette) +
    guides(fill = guide_legend(ncol = min(ncol(cpm), 5)))

  ppi <- 600
  png(filename = "ProcessCounts.png",
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)

  .grid_arrange_shared_legend(lcp, lcp_fil,
                              nrow = 1, ncol = 2,
                              position = "bottom")

  graphics.off()

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts,
                  counts_fil = counts_fil),
    colData = SummarizedExperiment::colData(se),
    metadata = S4Vectors::metadata(se))

  write.table(fil, file = "filtered_counts.txt", col.names = TRUE,
              row.names = TRUE, quote = FALSE, sep = "\t")

  cat("\n", sum(rowSums(is.na(counts_fil)) == ncol(counts_fil)),
      "features were removed by filtering.\n",
      sum(rowSums(is.na(counts_fil)) != ncol(counts_fil)),
      "features remain.\n\n")
  Sys.sleep(1)
  cat("\nHead of filtered counts matrix:\n")
  print(head(fil))

  Sys.sleep(1)
  cat("\nTail of filtered counts matrix:\n")
  print(tail(fil))
  return(se)
}
