#' Heatmap of most differential genes.
#'
#' Heatmap showing centred and normalised expression across each sample
#' for the most differentially-expressed genes in each contrast.
#'
#' @param results Results slot from object created by \code{\link{edger_pipe}},
#'                \code{\link{limma_pipe}}, or \code{\link{deseq2_pipe}}.
#' @param n Topmost number of genes to plot.
#' @param main Title of plot.
#' @param show_row_names Boolean indicating whether to show gene names in the
#'                       heatmap.
#'
#' @return
#' Heatmap in PNG format.
#'
#' @import ComplexHeatmap
#'
#' @export
#' @examples
#' \dontrun{
#' load(system.file("data", "er.rda", package = "dealr"))
#' create_hmp(er$results)
#' }

create_hmp <- function(results,
                       n = 100,
                       main = "Heatmap of top DEGs",
                       show_row_names = FALSE) {

  exp <- results
  if ("adj.P.Val" %in% colnames(exp[[1]])) {
    exp <- lapply(exp, function(x) {
      rownames(x) <- x$genes
      x[order(x$`adj.P.Val`), ]
      x[1:n, -which(names(x) %in%
                      c("lfcse", "stat", "baseMean", "logFC",
                        "AveExpr", "t", "P.Value", "adj.P.Val",
                        "B", "logCPM", "F", "PValue", "genes"))]
    })
  } else {
    exp <- lapply(exp, function(x) {
      rownames(x) <- x$genes
      x[order(x$FDR), ]
      x[1:n, -which(names(x) %in%
                      c("lfcse", "stat", "baseMean", "logFC",
                        "AveExpr", "t", "P.Value", "FDR",
                        "B", "logCPM", "F", "PValue", "genes"))]
    })
  }

  exp <- lapply(exp, function(x) {
    x <- as.matrix(sapply(x, as.numeric))
    x <- t(scale(t(x)))
  })

  for (i in 1:length(exp)) {
    ht <- ComplexHeatmap::Heatmap(exp[[i]], width = unit(120, "mm"),
                  show_row_names = show_row_names,
                  col = c("yellow2", "firebrick2"),
                  name = "row z score",
                  column_title = gsub("Contrast", "", names(exp)[i]),
                  column_title_gp = gpar(font = rep(2, length(exp)),
                                         fontsize = 10),
                  column_names_gp = gpar(fontsize = 6),
                  row_names_gp = gpar(font = rep(1, n), fontsize = 7),
                  show_heatmap_legend = TRUE)

    ppi <- 600
    png(filename = paste0("TopHeatMap_", names(exp)[i], ".png"),
        width = 8.4*ppi, height = 6.5*ppi, res = ppi)

    ComplexHeatmap::draw(ht)
    graphics.off()
  }
}
