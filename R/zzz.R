.onLoad <- function(libname = find.package("dealr"), pkgname = "dealr"){

  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      c("design.matrix", "contrast.levels", "samples.data",
        "contrast_levels", "counts.matrix", "B", "Bcount",
        "Bpos", "Bstd", "Consensus", "Guess", "Log-odds",
        "Log-pval", "Multiplicity", "Name", "Occurence", "P",
        "Placeholder", "Statistics", "StrandBias", "Tcount", "Tpos",
        "Tstd", "files", "path_to_samples", "contrast.matrix", "value",
        "variable", "ggplot2", "grid", "gridExtra", "arrangeGrob", "plot_kmer",
        "par", "mtext", "Count", "FAIL", "PASS", "Sample", "Var2", "WARN",
        "abline", "annotation", "combn", "contrast", "dev.off", "formula",
        "graphics.off", "head", "model.matrix", "p.adjust", "png", "read.table",
        "sambamba", "setNames", "status", "tail", "write.table", "x", "y",
        "vennCounts", ".", "out", "divide_by")
    )
  invisible()
}
