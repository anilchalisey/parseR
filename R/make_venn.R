#' Make venn diagrams showing overlapping gene lists
#'
#' @param results Object produced by \code{run_dealr()} function.
#' @param p.value p-value for determining signficance for differential
#'   expression.
#'
#' @return
#' Venn diagrams showing overlap of gene lists between edgeR, limma, and DESeq2
#'
#' @export

make_venn <- function(results = results, p.value = 0.05) {
  DE <- results$diff_genes
  for (i in seq_along(DE)) {
    if (names(DE)[i] != "limma") {
      DE[[i]] <- lapply(DE[[i]], function(x) x <- x[x$FDR <= p.value, ])
    }
    if (names(DE)[i] == "limma") {
      DE[[i]] <- lapply(DE[[i]], function(x) x <- x[x$`adj.P.Val` <= p.value, ])
    }
  }
  cDE <- do.call(Map, c(list, DE))
  cDE <- lapply(cDE, function(x) {
    lapply(x, function(y) {
      y$genes
    })
  })

  venns <- list()
  if (length(cDE) == 2) {
    for (i in seq_along(cDE)) {
      ppi <- 600
      png(filename = paste0(names(cDE)[i], "_venn.png"),
          width = 8.4*ppi, height = 6.5*ppi, res = ppi)
      venns[[i]] <- .venn2(cDE[[i]][[1]], cDE[[i]][[2]],
                           names(cDE[[i]]), names(cDE)[i])
      graphics.off()
    }
  }
  if (length(cDE) == 3) {
    for (i in seq_along(cDE)) {
      ppi <- 600
      png(filename = paste0(names(cDE)[i], "_venn.png"),
          width = 8.4*ppi, height = 6.5*ppi, res = ppi)
      venns[[i]] <- .venn3(cDE[[i]][[1]], cDE[[i]][[2]], cDE[[i]][[3]],
                           names(cDE[[i]]), names(cDE)[i])
      graphics.off()
    }
  }
  return(venns)
}
