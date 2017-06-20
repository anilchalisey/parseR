#' Run edgeR-QLF
#'
#' Function to run edgeR-QLF on a DGE object.
#'
#' @param dge DGE object produced by \code{\link[edgeR]{DGEList}}.
#' @param design Design matrix produced by \code{create_design}.
#' @param contrast Contrast matrix produced by \code{create_contrast}.
#' @param adjust.method Method to be used for adjustment of nominal p-values.
#'   May be one of "BH", "bonferroni", "holm", "hochberg", "hommel", "BY".
#' @param p.value Value between 0 and 1.  Adjusted p-value for the differential
#'   expression analysis.
#'
#' @return
#' Results of edgeR analysis
#'
#' @import edgeR
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qlf <- run_edger(dge, design, contrast.matrix)
#' }

run_edger <- function(dge = dge,
                      design = design,
                      contrast = contrast.matrix,
                      adjust.method = "BH",
                      p.value = 0.05) {

  y <- edgeR::estimateDisp(dge, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  ppi <- 600
  png(filename = "DispersionPlots.png",
      width = 8.4*ppi,
      height = 6.5*ppi,
      res = ppi)

  par(mfrow = c(1,2))
  edgeR::plotBCV(y)
  edgeR::plotQLDisp(fit)

  graphics.off()

  qlf <- list()
  for (i in seq_along(colnames(contrast))) {
    qlf[[i]] <- edgeR::glmQLFTest(fit, contrast = contrast[, i])
  }
  names(qlf) <- colnames(contrast)

  summar <- do.call(cbind, lapply(qlf, function(x) {
    summary(edgeR::decideTestsDGE(object = x,
                                  adjust.method = adjust.method,
                                  p.value = p.value, lfc = 0))
  }))
  colnames(summar) <- gsub("Contrast", "", names(qlf))

  results <- list()
  for (i in 1:length(qlf)) {
    results[[i]] <- edgeR::topTags(qlf[[i]], n = Inf)$table
  }
  names(results) <- c(gsub("Contrast", "", names(qlf)))
  for (i in 1:(length(results))) {
    results[[i]] <- merge(results[[i]], (as.data.frame(dge$counts)),
                          by.x = 1, by.y = 0, all.x = T)
    results[[i]] <- results[[i]][order(results[[i]]$FDR), ]
  }
  create_maplot(qlf, analysis = "edger")
  return(list(summar, results))
}
