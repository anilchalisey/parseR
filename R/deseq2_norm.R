#' Count normalisation using DESeq2
#'
#' Function to show the effect of the count normalisation that will be
#' performed as part of the DESeq2 analysis.
#'
#' @param dds DESeq2 dds object created by
#'            \code{\link[DESeq2]{DESeqDataSetFromMatrix}}.
#'
#' @return
#' Boxplots showing the effect of normalisation in PNG and PDF format.
#'
#' @importFrom DESeq2 estimateSizeFactors counts
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' load(system.file("data", "ddsMat.rda", package = "dealr"))
#' deseq2_norm(dds)
#' }

deseq2_norm <- function(dds = dds) {
  ppi <- 600
  png(filename = "Normalisation.png",
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)

  normdds <- DESeq2::estimateSizeFactors(dds)
  lcpm <- log2(DESeq2::counts(normdds, normalized = FALSE) + 0.5)
  lcpm <- suppressMessages(reshape2::melt(lcpm))
  lcpm$Var2 <- paste0(lcpm$Var2, "     ")
  lcp <-
    ggplot(lcpm, aes(x = Var2, y = value, fill = Var2)) +
    geom_boxplot() + .theme_Publication() +
    ggtitle("A. Before normalisation") +
    ylab(expression(log[2]*~cpm)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key.size = unit(0.7, 'lines'),
          legend.direction = "horizontal",
          legend.title = element_blank()) +
    scale_fill_manual(values = col.palette) +
    guides(fill = guide_legend(ncol = min(ncol(normdds), 5)))

  lcpm_norm <- log2(DESeq2::counts(normdds, normalized = TRUE) + 0.5)
  lcpm_norm <- suppressMessages(reshape2::melt(lcpm_norm))
  lcpm_norm$Var2 <- paste0(lcpm_norm$Var2, "     ")
  lcp_norm <-
    ggplot(lcpm_norm, aes(x = Var2, y = value, fill = Var2)) +
    geom_boxplot() + .theme_Publication() +
    ggtitle("B. After normalisation") +
    ylab(expression(log[2]*~cpm)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.key.size = unit(0.7, 'lines'),
          legend.direction = "horizontal",
          legend.title = element_blank()) +
    scale_fill_manual(values = col.palette) +
    guides(fill = guide_legend(ncol = min(ncol(normdds), 5)))

  .grid_arrange_shared_legend(lcp, lcp_norm,
                              nrow = 1, ncol = 2,
                              position = "bottom")
  graphics.off()
  print(data.frame(Sample = colnames(normdds),
                   NormFacs = normdds@colData$sizeFactor))
}
