#' Gene ontology and pathway analysis of DE genes
#'
#' @description \code{run_goseq} is a function to perform gene ontology analysis
#' using the goseq tool.
#'
#' @param universe List of all genes tested for differential expression.
#' @param degs List of differentially expressed genes.
#' @param test.cats GO categories to test [default = c("GO:BP", "GO:MF")]
#' @param numInCat Maximum number of genes in a GO term.
#' @param qv Threshold q-value for significance.
#' @param out Output directory.
#'
#' @return
#' \code{run_goseq} returns a dataframe and tab-delimited file of significant
#' GO terms.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' universe <- unlist(as.list(rownames(counts)))
#' de.genes <- er$results$`OX-BUFF`$genes
#' gor <- run_goseq(universe, degs = de.genes)
#'
#' gomdbr <- run_goseqmsigdb(universe, degs = de.genes, symbol = "HGNC")
#' }

run_goseq <- function(universe = universe,
                      degs = degs,
                      test.cats = c("GO:BP", "GO:MF"),
                      numInCat = 1000,
                      qv = 0.05,
                      out = ".") {

  .create_dir(out)

  y <- as.integer(universe %in% degs)
  y <- setNames(y, universe)

  ppi <- 600
  png(filename = file.path(out, "PWFPlot.png"),
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)
  suppressMessages(suppressWarnings(
    pwf <- goseq::nullp(y, "hg19", "geneSymbol", plot.fit = T)
    ))
  graphics.off()

  suppressMessages(suppressWarnings(
    gw <- goseq::goseq(pwf, "hg19", "geneSymbol", use_genes_without_cat = TRUE,
                       test.cats = test.cats)
    ))

  gw$qValue <- p.adjust(gw$over_represented_pvalue,
                        method = "BH", n = nrow(gw))
  gw <- gw[gw$numInCat < numInCat & gw$qValue < qv, ]
  rownames(gw) <- NULL
  colnames(gw) <- c("GO:ID", "pValue", "underpValue", "DEGenes",
                    "numInCat", "GO:term", "Ontology", "adj.pValue")
  gw <- gw[, c(1, 6, 7, 4, 5, 2, 8)]
  write.table(gw, file = file.path(out, "goseqanalysis.txt"), row.names = F,
              col.names = T, sep = "\t", quote = F)
  return(gw)
}

#' @description \code{run_goseqmdb} is a function to perform pathway analysis
#'   using the goseq tool.
#'
#' @param symbol Character specifying whether the gene identifiers are hgnc gene
#'   symbols ("HGNC") or entrezid ("ENTREZ").
#'
#' @return
#' \code{run_goseqmsigdb} returns a dataframe and tab-delimited file of
#' significant terms from the Broad Institute's Molecular Signatures Database.
#'
#' @rdname run_goseq
#'
#' @export

run_goseqmsigdb <- function(universe = universe,
                            degs = degs,
                            numInCat = 1000,
                            qv = 0.05,
                            out = "",
                            symbol = "HGNC") {

  .create_dir(out)
  if (symbol == "HGNC") {
    universe <- geneidcon[geneidcon$Approved.Symbol %in% universe, ]$entrez_id
    degs <- geneidcon[geneidcon$Approved.Symbol %in% degs, ]$entrez_id
  }

  y <- as.integer(universe %in% degs)
  y <- setNames(y, universe)

  ppi <- 600
  png(filename = file.path(out, "PWFPlot.png"),
      width = 8.4*ppi, height = 6.5*ppi, res = ppi)
  suppressMessages(suppressWarnings(
    pwf <- goseq::nullp(y, "hg19", "knownGene", plot.fit = T)
  ))
  graphics.off()

  suppressMessages(suppressWarnings(
    gw <- goseq::goseq(pwf, "hg19", "knownGene", gene2cat = MSigDB,
                       use_genes_without_cat = TRUE)))

  gw$qValue <- p.adjust(gw$over_represented_pvalue,
                        method = "BH", n = nrow(gw))
  gw <- gw[gw$numInCat < numInCat & gw$qValue < qv, ]
  rownames(gw) <- NULL
  colnames(gw) <- c("Term", "pValue", "underpValue", "DEGenes",
                    "numInCat", "adj.pValue")
  gw <- gw[, c(1, 4, 5, 2, 6)]
  write.table(gw, file = file.path(out, "msigdb_goseqanalysis.txt"),
              row.names = F, col.names = T, sep = "\t", quote = F)
  return(gw)
}
