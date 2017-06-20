#' Create MA plot
#'
#' Function to create Mean-Average (Bland-Altmann) plots showing the log-fold
#' changes in gene expression versus the average expression of that gene.
#'
#' @param res The output of either \code{\link{run_edger}},
#'   \code{\link{run_limma}}, or \code{\link{run_deseq2}}.
#' @param analysis The analysis method - one of "edger", "voom", or "deseq2".
#'
#' @return
#' MA plot in PNG and PDF formats.
#'
#' @importFrom limma plotMD decideTests
#' @importFrom edgeR decideTestsDGE
#' @import DESeq2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' create_maplot(results, analysis = "edger")
#' }

create_maplot <- function(res, analysis = "limma") {

  ppi <- 600
  png(filename = "MAPlot.png",
      width = 8.4*ppi,
      height = 6.5*ppi,
      res = ppi)
  if (analysis == "limma") {
    if (ncol(res) < 3) {
      par(mfrow = c(ceiling(ncol(res)/3), ceiling(ncol(res))),
          oma = c(0, 0, 2, 0))
    } else {
      par(mfrow = c(ceiling(ncol(res)/3), 3),
          oma = c(0, 0, 2, 0))
    }
    for (i in 1:ncol(res)) {
      limma::plotMD(res, column = i, status = limma::decideTests(res)[, i],
                    xlab = "Mean log-CPM across all samples",
                    ylab = "Expression log ratio", legend = FALSE,
                    main = gsub("Contrast", "", colnames(res))[[i]],
                    cex = 0.5, hl.col = "red3")
      abline(h = c(-1, 0, 1), col = "blue", lty = 2, lwd = 0.6)
        mtext("Mean-Difference Plot", outer = T)
      }
    }
    if (analysis == "edger") {
      if (length(res) < 3) {
        par(mfrow = c(ceiling(length(res)/3), ceiling(length(res))),
            oma = c(0, 0, 2, 0))
      } else {
        par(mfrow = c(ceiling(length(res)/3), 3),
            oma = c(0, 0, 2, 0))
      }
      for (i in 1:length(res)) {
        limma::plotMD(res[[i]], status = edgeR::decideTestsDGE(res[[i]]),
                      xlab = "Mean log-CPM across all samples",
                      ylab = "Expression log ratio", legend = FALSE,
                      main = gsub("Contrast", "", names(res)[[i]]),
                      cex = 0.5, hl.col = "red3")
        abline(h = c(-1, 0, 1), col = "blue", lty = 2, lwd = 0.6)
        mtext("Mean-Difference Plot", outer = T)
      }
    }
    if (analysis == "deseq2") {
      if (length(res) < 3) {
        par(mfrow = c(ceiling(length(res)/3), ceiling(length(res))),
            oma = c(0, 0, 2, 0))
      } else {
        par(mfrow = c(ceiling(length(res)/3), 3),
            oma = c(0, 0, 2, 0))
      }
      for (i in 1:length(res)) {
        DESeq2::plotMA(res[[i]], cex = 0.5, main = names(res)[i])
        abline(h = c(-1, 0, 1), col = "blue", lty = 2, lwd = 0.6)
        mtext("Mean-Difference Plot", outer = T, cex = 1.5, font = 2)
      }
    }
    graphics.off()
}


