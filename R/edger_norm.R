#' Count normalisation using edgeR
#'
#' Function to perform count normalisation prior to differential analysis,
#' with generation of before and after plots.
#'
#' @param dge DGE object.
#' @param norm.method Normalisation method to use - may be "TMM", "RLE",
#'   "upperquartile" or "none".
#'
#' @return
#' Normalised DGE object.
#'
#' @import ggplot2
#' @import edgeR
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dge <- edger_norm(dge, norm.method = "TMM")
#' }

edger_norm <- function(dge = NULL, norm.method = "TMM") {

  lcpm <- edgeR::cpm.default(dge, log = TRUE)
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
    guides(fill = guide_legend(ncol = min(ncol(dge), 5)))

  dge_norm <- edgeR::calcNormFactors(dge, method = norm.method)
  lcpm_norm <- edgeR::cpm(dge_norm, log = TRUE)
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
    guides(fill = guide_legend(ncol = min(ncol(dge_norm), 5)))

  ppi <- 600
  png(filename = "Normalisation.png",
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)

  .grid_arrange_shared_legend(lcp, lcp_norm,
                              nrow = 1, ncol = 2,
                              position = "bottom")
  graphics.off()

  print(data.frame(Sample = colnames(dge_norm),
                   NormFacs = dge_norm$samples$norm.factors))
  return(dge_norm)
}




