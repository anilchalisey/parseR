#' Create multi-dimensional scaling plot from DGE object
#'
#' Function to create a MDS plot from a DGE object.
#'
#' @param dge DGE object produced by \code{\link[edgeR]{DGEList}}.
#' @param main Plot title.
#' @param dim.plot Dimensions to plot - default are dimensions 1 and 2.
#' @param n Number of pairwise genes to consider when plotting MDS.  Default is
#'   500.
#'
#' @return
#' MDS plot.
#'
#' @importFrom limma plotMDS
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' create_mds(dge)
#' }

create_mds <- function(dge = dge,
                       main = "",
                       dim.plot = c(1,2),
                       n = min(500, nrow(dge$counts))) {

  plotmds.invisible <- function(...) {
    ff <- tempfile()
    png(filename = ff)
    res <- limma::plotMDS(...)
    dev.off()
    unlink(ff)
    res
  }

  mds <- plotmds.invisible(dge, top = n, dim.plot = dim.plot)
  abs <- range(mds$x)
  abs <- abs(abs[2] - abs[1])/25
  ord <- range(mds$y)
  ord <- abs(ord[2] - ord[1])/25

  ppi <- 600
  png(filename = "MDSPlot.png",
      width = 8.4*ppi,
      height = 6.5*ppi,
      res = ppi)

  mdsp <- data.frame(x = mds$x, y = mds$y,
                     contrast = dge$samples$group)
  mdspp <-
    ggplot(mdsp, aes(x, y, group = contrast)) +
    geom_point(aes(fill = contrast, shape = contrast),
               colour = "black", pch = 21, size = 4) +
    .theme_Publication() + scale_fill_manual(values = col.palette) +
    xlab(paste0("Leading logFC dimension ", dim.plot[1])) +
    ylab(paste0("Leading logFC dimension ", dim.plot[2])) +
    geom_vline(xintercept = 0, lty = 3, alpha = 0.5) +
    geom_hline(yintercept = 0, lty = 3, alpha = 0.5) +
    geom_text(aes(label = rownames(mdsp)),
              x = mdsp$x - ifelse(mdsp$x > 0, abs, -abs),
              y = mdsp$y - ifelse(mdsp$y > 0, ord, -ord)) +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(ncol = min(ncol(dge), 5)))
  print(mdspp)
  graphics.off()
}

