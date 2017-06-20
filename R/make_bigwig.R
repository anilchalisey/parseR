#' Create bigwig files from BAM
#'
#' @description Script to generate bigwig files with the option for
#' normalisation from a BAM file.  Currently will only work with paired end
#' BAM files.
#'
#' @param bamfile Path to a BAM file.
#' @param scale The number of reads to scale to.  If NULL no scaling performed.
#'
#' @return
#' Coverage file in bigwig format.
#'
#' @importFrom Rsamtools idxstatsBam BamFile
#' @importFrom GenomicAlignments readGAlignmentPairs coverage
#' @importFrom rtracklayer export.bw
#'
#' @export
#'
#' @examples
#' \dontrun{
#' make_bigwig(bamfile = "HB1_sample.bam", scale = 20000000)
#' }
#'
make_bigwig <- function(bamfile = NULL,
                        scale = 20000000) {

  idxstat <- Rsamtools::idxstatsBam(bamfile)
  read_count <- idxstat[, 3]
  read_count <- sum(read_count)
  ratio <- round(scale/read_count, digits = 2)

  # Name for bigwig file
  bname <- .remove_ext(bamfile)
  idx <- paste0(bamfile, ".bai")
  bigwig <- paste0(bname, ".bw")

  # Convert BAM to RLE
  bam <- Rsamtools::BamFile(file = bamfile,
                            index = idx,
                            yieldSize = 100000,
                            asMates = TRUE)

  open(bam)
  cvg <- NULL
  repeat {
    chunk <- GenomicAlignments::readGAlignmentPairs(bam)
    if (length(chunk) == 0L)
      break
    chunk_cvg <- GenomicAlignments::coverage(chunk)
    if (is.null(cvg)) {
      cvg <- chunk_cvg
    } else {
      cvg <- cvg + chunk_cvg
    }
  }
  close(bam)
  if (!is.null(scale)) cvg <- cvg*ratio
  rtracklayer::export.bw(object = cvg, con = bigwig)
}
