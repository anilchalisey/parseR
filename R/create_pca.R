#' Create principal component analysis plot from DDS object
#'
#' Function to create a PCA plot from a DDS object
#'
#' @param dds DESeq2 dds object created by
#'            \code{\link[DESeq2]{DESeqDataSetFromMatrix}}.
#'
#' @return
#' PCA plot in PNG and PDF formats.
#'
#' @importFrom DESeq2 rlog plotPCA
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' load(system.file("data", "ddsMat.rda", package = "dealr"))
#' create_pca(ddsMat)
#' }

create_pca <- function(dds) {

  rld <- DESeq2::rlog(dds, blind = FALSE)
  data <- DESeq2::plotPCA(rld,
                          intgroup = "Contrast",
                          returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))

  abs = range(data$PC1)
  abs = abs(abs[2] - abs[1])/25
  ord = range(data$PC2)
  ord = abs(ord[2] - ord[1])/25

  ppi <- 600
  png(filename = "PCAPlot.png",
      width = 8.4*ppi,
      height = 6.5*ppi,
      res = ppi)

  pcap <- data.frame(x = data$PC1, y = data$PC2,
                     contrast = colData(dds)$Contrast)
  rownames(pcap) <- colnames(dds)
  pcapp <-
    ggplot(pcap, aes(x, y, group = contrast)) +
    geom_point(aes(fill = contrast, shape = contrast),
               colour = "black", pch = 21, size = 4) +
    .theme_Publication() + scale_fill_manual(values = col.palette) +
    xlab(paste0("PC1: ", percentVar[1], "% of variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% of variance")) +
    geom_vline(xintercept = 0, lty = 3, alpha = 0.5) +
    geom_hline(yintercept = 0, lty = 3, alpha = 0.5) +
    geom_text(aes(label = rownames(pcap)),
              x = pcap$x - ifelse(pcap$x > 0, abs, -abs),
              y = pcap$y - ifelse(pcap$y > 0, ord, -ord)) +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(ncol = min(ncol(dds), 5)))
  print(pcapp)
  graphics.off()
}

