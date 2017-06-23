#' Run the differential analysis pipeline
#'
#' Function to run the differential analysis portion of the pipeline.  There are
#' two possible starting points:
#' \itemize{
#'   \item{BAM files}{If starting with aligned BAM files, then set count = T.
#'         The pipeline will then begin with read counting before progressing to
#'         differential analysis.}
#'   \item{Count files}{If starting with count files generated externally, then
#'         set count = F and specify the text file containing the matrix of
#'         counts.  The pipeline will proceed directly to the differential
#'         analysis. The row names of the matrix should indicate the gene and
#'         the column name the library.  It is important that column names match
#'         the baseline in the sample metadata file.}
#' }
#'
#' @param threads Number of threads to utilise where parallel
#'   processing is possible.
#' @param experimentTitle Name of the experiment.
#' @param modules The differential expression analysis method
#'   to be used. May be any combination of 'L' (limma), 'E' (edgeR), and 'D'
#'   (DESeq2).
#' @param p.value Value between 0 and 1 specifying the adjusted
#'   p.value threshold for significance.
#' @param sample.path Full path to the sample metadata file.
#' @param contrast.column  Column in sample metadata file containing
#'   information about the contrasts of interest.
#' @param block.column Column in sample metadata file containing
#'   information about blocking factors or batch effect.
#' @param contrast.levels The order in which the contrasts should be evaluated.
#' @param count Boolean (TRUE or FALSE) indicating whether to perform counting.
#' @param count.file The path to the file containing the count matrix
#'   (if count = F).  The first column should be the gene/feature names and the
#'   subsequent columns the counts of reads in each sample.  The sample column
#'   names must match the names in the sample metadata table.
#' @param featurecounts The path to featureCounts (if not in executable path).
#' @param bamfiles Vector of full path to bam files to be counted (if count = T).
#' @param annotationFile Path to annotation file.  If NULL, the default hg19
#'   genome is used.
#' @param annotationFormat Specify the format of the annotation file. Acceptable
#'   formats include ‘GTF’ and ‘SAF’.  The in-built annotation is 'SAF'.
#' @param requireBothEndsMapped Logical.  If TRUE, only fragments that have both
#'   ends successfully aligned will be considered for summarization. This option
#'   should be used together with pairedEnd = TRUE.
#' @param excludeChimeric Logical.  If TRUE, the chimeric fragments (those
#'   fragments that have their two ends aligned to different chromosomes) will
#'   NOT be counted. This option should be used together with pairedEnd = TRUE.
#' @param pairedEnd Logical. If TRUE,  fragments (or templates) will be counted
#'   instead of reads. This option is only applicable for paired-end reads.
#' @param countMultiMapping If TRUE, multi-mapping reads/fragments will be
#'   counted.
#' @param multiFeatureReads Reads/fragments overlapping with more than one
#'   meta-feature/feature will be counted more than once. Note that when
#'   performing meta-feature level summarization, a read (or fragment) will
#'   still be counted once if it overlaps with multiple features within the same
#'   meta-feature (as long as it does not overlap with other metafeatures).
#' @param ignoreDup Logical.  If TRUE, reads that were marked as duplicates will
#'   be ignored.  In paired end data, the entire read pair will be ignored if at
#'   least one end is found to be a duplicate read.
#' @param outname Character string.  Name of the output file. The output file
#' contains the number of reads assigned to each meta-feature or feature.
#' @param cpm.cutoff The counts-per-million threshold for filtering counts.
#' @param design Custom design matrix that may be used in place of default
#'    matrix created by pipeline.  Only used for edgeR and limma analyses.
#' @param contrast Custom contrast matrix that may be used in place of default
#'    matrix created by pipeline.  Only used for edgeR and limma analyses.
#' @param norm.method Normalisation method used by edgeR/limma analyses.  May
#'    be "TMM", "RLE", "upperquartile" or "none".
#' @param adjust.method Method to be used for adjustment of nominal p-values.
#'   May be one of "BH", "bonferroni", "holm", "hochberg", "hommel", "BY".
#'
#' @return
#' Results of differential analysis.
#'
#' @export

run_diff <- function(
  # General options - important
  threads = 1,
  experimentTitle = "",
  modules = "LED",
  p.value = 0.05,
  # sample metadata options - important
  sample.path = NULL,
  contrast.column = NULL,
  block.column = NULL,
  contrast.levels = NULL,
  # Counting options - important
  count = TRUE,
  count.file = NULL,
  featurecounts = "featureCounts",
  bamfiles = NULL,
  # Counting options - default options suitable for most
  annotationFile = NULL,
  annotationFormat = NULL,
  requireBothEndsMapped = TRUE,
  excludeChimeric = TRUE,
  pairedEnd = TRUE,
  countMultiMapping = FALSE,
  multiFeatureReads = TRUE,
  ignoreDup = FALSE,
  outname = "counts.txt",
  # Count filtering options - default option suitable for most
  cpm.cutoff = 1,
  # Experimental design options - default option suitable for most
  design = NULL,
  contrast = NULL,
  norm.method = "TMM",
  adjust.method = "BH") {

  if (!isTRUE(count) & is.null(count.file)) {
    stop("count is FALSE, but a count file has not been specified.")
  }
  if (isTRUE(count) & is.null(bamfiles)) {
    stop("count is TRUE, but bam files have not been specified.")
  }
  if (isTRUE(count) & .check_cmd(featurecounts) == "Not Found") {
    stop("the featureCounts path is invalid or it is not installed correctly.")
  }

  ############
  # Load up the sample metadata
  ############
  sample.data <- load_samples(sample.path = sample.path,
                              contrast.column = contrast.column,
                              block.column = block.column,
                              contrast.levels = contrast.levels)

  ############
  # Create tree folder structure
  ############
  .create_folders(modules = modules)

  ############
  # Raw count data
  ############
  if (isTRUE(count)) {
    counts <- count_reads(featurecounts = featurecounts,
                          annotationFile = annotationFile,
                          annotationFormat = annotationFormat,
                          requireBothEndsMapped = requireBothEndsMapped,
                          excludeChimeric = excludeChimeric,
                          pairedEnd = pairedEnd,
                          countMultiMapping = countMultiMapping,
                          multiFeatureReads = multiFeatureReads,
                          threads = threads, ignoreDup = ignoreDup,
                          outname = outname,
                          alignments = bamfiles)
  } else {
    counts <- data.table::fread(input = count.file, data.table = FALSE)
    rownames(counts) <- as.list(unlist(counts[1]))
    counts <- as.matrix(counts[, -1])
  }

  ############
  # Create SummarizedExperiment
  ############
  se <- create_se(samples = sample.data, counts = counts,
                  experimentTitle = experimentTitle)

  ############
  # Filter the data and show the effect
  ############
  se <- process_counts(se = se, cpm.cutoff = 1)

  ############
  # Some tidying up
  ############
  unlink("hg19_genes.saf")
  tomove <- c("ProcessCounts.png", "filtered_counts.txt")
  if (isTRUE(count)) tomove <- c(tomove, "counts.txt", "counts.txt.full",
                                 "counts.txt.full.summary")
  lapply(tomove, function(x) {
    file.rename(x, file.path("DE_results/Common_results/", x))
  })

  ############
  # Setting up the differential analysis
  ############
  results <- list()
  results$samples <- colnames(se)
  results$summary <- list()
  results$diff_genes <- list()

  modules_selected <- unlist(strsplit(modules, ""))

  ############
  # limma analysis
  ############

  if ("L" %in% modules_selected) {
    cat(paste('\n Gene expression analysis is performed with limma + voom.\n\n'))
    limma_results <- limma_pipe(se = se, design = design, contrast = contrast,
                                block.column = block.column,
                                norm.method = norm.method,
                                adjust.method = adjust.method,
                                p.value = p.value)
    results$summary$limma <- limma_results[[1]]
    results$diff_genes$limma <- limma_results[[2]]
    create_hmp(results = limma_results[[2]])

    ############
    # Tidying up
    ############

    patterns <- c("DEgenes*", "*.png")
    tomove <- lapply(patterns, function(x) list.files(pattern = x))
    lapply(tomove, function(x) {
      file.rename(from = x, to = file.path("DE_results/limma_results", x))
    })
  }

  ############
  # edgeR analysis
  ############

  if ("E" %in% modules_selected) {
    cat(paste('\n Gene expression analysis is performed with edgeR.\n\n'))
    edger_results <- edger_pipe(se = se, design = design, contrast = contrast,
                                block.column = block.column,
                                norm.method = norm.method,
                                adjust.method = adjust.method,
                                p.value = p.value)
    results$summary$edger <- edger_results[[1]]
    results$diff_genes$edger <- edger_results[[2]]
    create_hmp(results = edger_results[[2]])

    ############
    # Tidying up
    ############

    patterns <- c("DEgenes*", "*.png")
    tomove <- lapply(patterns, function(x) list.files(pattern = x))
    lapply(tomove, function(x) {
      file.rename(from = x, to = file.path("DE_results/edger_results", x))
    })
  }

  ############
  # DESeq2 analysis
  ############

  if ("D" %in% modules_selected) {
    cat(paste('\n Gene expression analysis is performed with DEseq2.\n\n'))
    deseq2_results <- deseq2_pipe(se = se,
                                  block.column = block.column,
                                  adjust.method = adjust.method,
                                  p.value = p.value)
    results$summary$deseq2 <- deseq2_results[[1]]
    results$diff_genes$deseq2 <- deseq2_results[[2]]
    create_hmp(results = deseq2_results[[2]])

    ############
    # Tidying up
    ############

    patterns <- c("DEgenes*", "*.png")
    tomove <- lapply(patterns, function(x) list.files(pattern = x))
    lapply(tomove, function(x) {
      file.rename(from = x, to = file.path("DE_results/deseq2_results", x))
    })
  }

  ############
  # Combined analyses
  ############

  # genes that are identified in all modules for each comparison
  if (length(modules_selected) > 1) {
    cat(paste('\n Identifying genes determined as DEGs by all analyses.\n\n'))
    common_genes <- .extract_common(results = results, p.value = 0.05)
    for (i in seq_along(common_genes)) {
      write.table(common_genes[[i]],
                  file = paste0(names(common_genes)[i], "_common.txt"),
                  row.names = FALSE, col.names = TRUE, quote = FALSE,
                  sep = "\t")
    }
    vennplots <- make_venn(results = results, p.value = p.value)

    ############
    # Tidying up
    ############

    tomove <- list.files(pattern = "common.txt|venn.png")
    lapply(tomove, function(x) {
      file.rename(from = x, to = file.path("DE_results/Common_results/", x))
    })
  }

  cat("Analysis complete!\n\n")
  return(results)
}

