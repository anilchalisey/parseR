#' Count reads in bam files using featureCounts
#'
#' @description  Function to count reads mapping to a reference genome.  Users
#' can either supply their own reference or use the built-in one for hg19.
#'
#' @param featurecounts Character string.  Path to featureCounts
#' @param annotationFile Character string. Path to annotation file.  If NULL, the
#'   default hg19 genome is used.
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
#' @param pairInsertRange Vector of two integers specifying the minimum and
#'   maximum fragment/template lengths.  The default values are c(50, 600).
#'   Must be used together with pairedEnd = TRUE.
#' @param countMultiMapping If TRUE, multi-mapping reads/fragments will be
#'   counted.
#' @param multiFeatureReads Reads/fragments overlapping with more than one
#'   meta-feature/feature will be counted more than once. Note that when
#'   performing meta-feature level summarization, a read (or fragment) will
#'   still be counted once if it overlaps with multiple features within the same
#'   meta-feature (as long as it does not overlap with other metafeatures).
#' @param minQual The minimum mapping quality score a read must satisfy in
#'   order to be counted. For paired-end reads, at least one end should satisfy
#'   this criteria. 0 by default.
#' @param stranded Indicate if strand-specific read counting should be
#'   performed.  Acceptable values: 0 (unstranded), 1 (stranded) and 2
#'   (reversely stranded). 0 by default. For paired-end reads, strand of the
#'   first read is taken as the strand of the whole fragment.
#' @param groupBy Specify the attribute type used to group features (eg. exons)
#'   into meta-features (eg. genes) when GTF annotation is provided. Default is
#'   'gene id’.
#' @param featureType  Specify the feature type. Only rows which have the
#'   matched feature type in the provided GTF annotation file will be included
#'   for read counting. Default is ‘exon’.
#' @param threads Number of the threads. The value should be between 1 and
#'   32. 1 by default.
#' @param ignoreDup Logical.  If TRUE, reads that were marked as duplicates will
#'   be ignored.  In paired end data, the entire read pair will be ignored if at
#'   least one end is found to be a duplicate read.
#' @param outname Character string.  Name of the output file. The output file
#' contains the number of reads assigned to each meta-feature or feature.
#' @param alignments Character vector.  Paths to BAM files.
#'
#' @return
#' Tab-delimited file containing the number of reads assigned to each
#' meta-feature or feature and a matrix with the same.
#'
#' @importFrom data.table fwrite fread
#'
#' @export
#'
#' @examples
#' \dontrun{
#' alignments <- list.files(path = "bam_files", pattern = "*.bam$")
#' threads <- parallel::detectCores() - 1
#' counts <- count_reads(featurecounts = "featureCounts", threads = threads,
#'                       alignments = alignments)
#' }



count_reads <- function(featurecounts = "featureCounts",
                        annotationFile = NULL,
                        annotationFormat = "SAF",
                        requireBothEndsMapped = TRUE,
                        excludeChimeric = TRUE,
                        pairedEnd = TRUE,
                        pairInsertRange = NULL,
                        countMultiMapping = FALSE,
                        multiFeatureReads = FALSE,
                        minQual = 0,
                        stranded = 0,
                        groupBy = NULL,
                        featureType = NULL,
                        threads = 1,
                        ignoreDup = FALSE,
                        outname = NULL,
                        alignments = NULL) {

  if (is.null(annotationFile)) {
    annotationFormat <- "SAF"
    data.table::fwrite(x = annotation, file = "hg19_genes.saf",
                      col.names = T, sep = "\t", row.names = F, quote = F)
    annotationFile <- "hg19_genes.saf"
  }

  if (annotationFormat == "SAF") {
    groupBy <- NULL
    featureType <- NULL
  }

  if (is.null(outname)) outname <- "raw_counts.txt"

  feature.counts <- paste(featurecounts,
                          "-a", annotationFile,
                          "-F", annotationFormat,
                          if (isTRUE(requireBothEndsMapped)) {"-B"},
                          if (isTRUE(excludeChimeric)) {"-C"},
                          if (isTRUE(pairedEnd)) {"-p"},
                          if (!is.null(pairInsertRange)) {
                            paste("-P","-d", pairInsertRange[1],
                                  "-D", pairInsertRange[2])},
                          if (isTRUE(countMultiMapping)) {"-M"},
                          if (isTRUE(multiFeatureReads)) {"-O"},
                          "-Q", minQual,
                          "-s", stranded,
                          if (!is.null(groupBy)) {paste("-g", groupBy)},
                          if (!is.null(featureType)) {paste("-t", featureType)},
                          "-T", threads,
                          if (isTRUE(ignoreDup)) {"--ignoreDup"},
                          "-o", paste0(outname, ".full"),
                          paste(alignments, collapse = " "))
    res <- .cmd(feature.counts)
    results <- as.data.frame(data.table::fread(input = paste0(outname, ".full"),
                                               skip = 1))
    results <- results[, -c(2:6)]
    colnames(results) <- .remove_ext(colnames(results))
    rownames(results) <- results$Geneid
    results <- as.matrix(results[, -1])

    cat("\nTotal number of features:\n", nrow(results))
    Sys.sleep(1)
    cat("\nHead of counts matrix:\n")
    print(head(results))

    Sys.sleep(1)
    cat("\nTail of counts matrix:\n")
    print(tail(results))

    write.table(results, file = outname, col.names = NA, row.names = TRUE,
                sep = "\t", quote = FALSE)

    return(results)
}
