#' Run limma+voom
#'
#' Function to run voom normalisation (with quality weights)
#' followed by limma on a DGE object.
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
#' Results of limma analysis
#'
#' @importFrom limma lmFit contrasts.fit eBayes decideTests
#' @importFrom limma voomWithQualityWeights topTable
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' qlf <- run_limma(dge, design, contrast.matrix)
#' }

run_limma <- function(dge = dge,
                      design = design,
                      contrast = contrast,
                      adjust.method = "BH",
                      p.value = 0.05) {

  ppi <- 600
  png(filename = "limma.png",
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)
  elist <- limma::voomWithQualityWeights(dge, design, plot = T)
  graphics.off()

  vfit <- limma::lmFit(elist, design)
  vfit <- limma::contrasts.fit(vfit, contrasts = contrast)
  efit <- limma::eBayes(vfit)

  summar <- summary(limma::decideTests(efit))
  names(summar) <- gsub("Contrast", "", names(summar))

  limma_results <- list()
  for (i in 1:ncol(efit)) {
    limma_results[[i]] <- limma::topTable(efit, coef = i, n = Inf)
  }
  names(limma_results) <- c(gsub("Contrast", "", colnames(efit)))
  for (i in 1:(length(limma_results))) {
    limma_results[[i]] <- merge(limma_results[[i]], (as.data.frame(dge$counts)),
                          by = 0, all.x = T)
    limma_results[[i]]$Row.names <- NULL
    limma_results[[i]] <-
      limma_results[[i]][order(limma_results[[i]]$adj.P.Val), ]
  }
  create_maplot(efit, analysis = "limma")
  return(list(summar, limma_results))
}
