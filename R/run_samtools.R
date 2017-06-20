#' Wrapper scripts for SAMtools functions
#'
#' @description \code{run_samsort} sorts alignment (SAM/BAM/CRAM) files either by position
#' or read name using \code{samtools sort}.
#'
#' @param file A vector of characters specifying the path to the bam files.
#' @param samtools The path to samtools (if not in executable path).
#' @param outformat String specifying output format: ('SAM'/'BAM'/'CRAM').  For
#'   \code{run_samview} this should be SAM or 'BAM' only.
#' @param threads A positive integer specifying the number of sorting and
#'   compression threads
#' @param memory String specifying maximum memory per thread; suffix K/M/G
#'   recognized.
#' @param sortbyname Boolean. If TRUE, reads are sorted by name. If FALSE,
#' reads are sorted by chromosome/position.
#' @param suffix Suffix to add to the basename of the file e.g, _psort (for
#' position-sorted).
#' @param keep Boolean.  If TRUE, the input file is kept.  If FALSE, the input
#' file is deleted after a successful sort.
#'
#' @return \code{run_samsort} returns a sorted alignment file.
#'
#' @examples
#' \dontrun{
#' run_samsort(file = "HB1_sample.sam", samtools = "samtools",
#'             outformat = "BAM", threads = (parallel::detectCores() - 1),
#'             memory = "1G", sortbyname = FALSE, suffix = "", keep = TRUE)
#'
#' }
#'
#' @export

run_samsort <- function(file = NULL,
                        samtools = "samtools",
                        outformat = "BAM",
                        threads = 1,
                        memory = "768M",
                        sortbyname = FALSE,
                        suffix = "",
                        keep = TRUE) {

  if (!(outformat %in% c("SAM","BAM","CRAM"))) {
    stop("output format must be \"SAM\" \"BAM\" or \"CRAM\"")
  }

  if (tolower(outformat) == tools::file_ext(file) && suffix == "") {
    stop("output and input files must have different names - would you like to
         add a suffix?")
  }

  bname <- .remove_ext(file)
  outfile <- paste0(bname, suffix, ".", tolower(outformat))
  temp <- paste0(outfile, "_tmp")

  if (isTRUE(sortbyname)) {
    sam.sort <- sprintf("%s sort -n -@ %s -m %s -T %s -O %s -o %s %s",
                        samtools, threads, memory, temp, outformat,
                        outfile, file)
  } else {
    sam.sort <- sprintf("%s sort -@ %s -m %s -T %s -O %s -o %s %s",
                        samtools, threads, memory, temp, outformat,
                        outfile, file)
  }

  for (i in seq_along(sam.sort)) {
    .cmd(sam.sort[i])
  }

  if (!isTRUE(keep) && file.exists(outfile)) {unlink(file)}
  return(outfile)
}

#' @description \code{run_samindex} indexes sorted BAM files using
#' \code{samtools index}.
#' @param bamfile Vector of characters specifying the path to sorted BAM files.
#' @rdname run_samsort
#' @return \code{run_samindex} returns a BAI-format index file for sorted
#' BAM files.
#' @examples
#' \dontrun{
#' run_samindex(samtools = "samtools", bamfile = "HB1_sample.bam",
#' threads = (parallel::detectCores() - 1))
#' }
#'
#' @export

run_samindex <- function(samtools = "samtools",
                         bamfile = NULL,
                         threads = 1) {
  sam.index <- sprintf("%s index -b -@ %s %s",
                       samtools, threads, bamfile)

  for (i in seq_along(sam.index)) {
    .cmd(sam.index[[i]])
  }
}

#' @description \code{run_samview} outputs all alignments matching the flag and
#'   region filters specified in either SAM or BAM format using
#'   \code{samtools view}.
#' @param regions Either path to a single BED file containing regions of
#'   interest or regions defined in the format "RNAME:START[-END]".
#' @param chrom.sizes Path to chromosome size file. Required if converting to
#'   SAM format and including header.  If not available, then one will be
#'   generated automatically from the existing header.  This is a tab-delimited
#'   text file where each line contains the reference name in the first column
#'   and length in the second column.
#' @param include.flag Integer value.  Include output of alignments with any
#'   bits set in this value present in the FLAG field.  See
#'   \url{https:://broadinstitute.github.io}.
#' @param exclude.flag Integer value.  Exclude output of alignments with any
#'   bits set in this value present in the FLAG field.  See
#'   \url{https:://broadinstitute.github.io}.
#' @param minQual Skip alignments with MAPQ smaller than this value.
#' @param outname Name for output file.
#' @param include.header Boolean.  If TRUE header is included.  This option is
#'   only necessary if outformat = "SAM" as header always included in "BAM".
#' @param count Boolean.  If true, only count alignments instead of printing.
#' @param subsample Float value - if set then subsampling is performed.  The
#'   integer part is used to seed the random number generator and the part after
#'   the decimal sets the fraction of reads to subsample.
#' @param keep.paired Boolean.  If TRUE, keep paired reads only.
#' @param keep.proper.pair Boolean.  If TRUE, keep concordantly paired reads
#'   only.
#' @param remove.unmapped Boolean.  If TRUE, remove unmapped reads.
#' @param remove.not.primary Boolean.  If TRUE, remove reads mapped as secondary
#'   alignments.
#' @param remove.duplicates Boolean.  If TRUE, remove reads marked as optical or
#'   PCR duplicates.
#' @param remove.supplementary.alignment Boolean.  If TRUE, remove reads marked
#'   as supplementary alignment.
#' @param remove.mitochondrial Character string.  If set, this will remove reads
#'   mapping to the mitochondrial genome.  The string should match the reference
#'   name for the mitochindrial genome in the alignment file.  Examples include
#'   "ChrM", "M" and "MT".
#' @rdname run_samsort
#' @return \code{run_samview} returns aligned reads according to the regions and
#' filters specified in BAM or SAM format.
#' @examples
#' \dontrun{
#' run_samview(samtools = "samtools", file = "HB1_sample.bam", regions = NULL,
#' chrom.sizes = NULL, include.flag = NULL, exclude.flag = NULL,
#' minQual = NULL, outformat = "BAM", outname = NULL, include.header = FALSE,
#' count = FALSE, threads = (parallel::detectCores() - 1), subsample = NULL,
#' keep.paired = TRUE, keep.proper.pair = TRUE, remove.unmapped = TRUE,
#' remove.not.primary = TRUE, remove.duplicates = TRUE,
#' remove.supplementary.alignment = TRUE, remove.mitochondrial = "ChrM")
#' }
#'
#' @export


run_samview <- function(samtools = "samtools",
                        file = NULL,
                        regions = NULL,
                        chrom.sizes = NULL,
                        include.flag = NULL,
                        exclude.flag = NULL,
                        minQual = NULL,
                        outformat = NULL,
                        outname = NULL,
                        include.header = FALSE,
                        count = FALSE,
                        threads = 1,
                        subsample = NULL,
                        keep.paired = TRUE,
                        keep.proper.pair = TRUE,
                        remove.unmapped = FALSE,
                        remove.not.primary = FALSE,
                        remove.duplicates = FALSE,
                        remove.supplementary.alignment = FALSE,
                        remove.mitochondrial = NULL) {

  if (!(outformat %in% c("SAM","BAM","CRAM"))) {
    stop("output format must be \"SAM\" \"BAM\" or \"CRAM\"")
  }

  if (isTRUE(include.header) && is.null(chrom.sizes) && outformat == "BAM") {
    sam.idx <- sprintf("%s idxstats %s > chrom.sizes", samtools, file)
    chrom.sizes <- "chrom.sizes"
    on.exit(unlink("chrom.sizes"))
  }

  bname <- .remove_ext(file)
  if (is.null(outname)) outname <- paste0(bname, "_filtered.",
                                          tolower(outformat))

  if (!is.null(regions)) {
    if (file.exists(regions) & length(regions) == 1) {
      oneRegion = FALSE
      rstring <- paste("-L", regions)
    } else if (grepl(":", regions) & grepl("-", regions)) {
        oneRegion = TRUE
        rstring <- regions
      } else {
        stop("The region argument must  be a single bed file")
      }
    } else {
        oneRegion <- FALSE
      }

  sam.view <- sprintf("%s view", samtools)
  if (outformat == "BAM") sam.view <- sprintf("%s -b", sam.view)
  if (isTRUE(include.header)) sam.view <- sprintf("%s -h -t %s",
                                                  sam.view, chrom.sizes)
  if (isTRUE(count)) sam.view <- sprintf("%s -c", sam.view)
  if (!is.null(minQual)) sam.view <- sprintf("%s -q %s", sam.view, minQual)
  if (!is.null(include.flag)) sam.view <- sprintf("%s -f %s", sam.view,
                                                  include.flag)
  if (!is.null(exclude.flag)) sam.view <- sprintf("%s -F %s", sam.view,
                                                 exclude.flag)
  if (isTRUE(keep.paired)) sam.view <- sprintf("%s -f 2", sam.view)
  if (isTRUE(keep.proper.pair)) sam.view <- sprintf("%s -f 3", sam.view)
  if (isTRUE(remove.unmapped)) sam.view <- sprintf("%s -F 4", sam.view)
  if (isTRUE(remove.not.primary)) sam.view <- sprintf("%s -F 256", sam.view)
  if (isTRUE(remove.duplicates)) sam.view <- sprintf("%s -F 1024", sam.view)
  if (isTRUE(remove.supplementary.alignment)) sam.view <- sprintf("%s -F 2056",
                                                                  sam.view)
  if (!is.null(subsample)) sam.view <- sprintf("%s -s %s", sam.view, subsample)
  if (!is.null(regions) & !isTRUE(oneRegion)) sam.view <- sprintf("%s %s",
                                                                  sam.view,
                                                                  rstring)
  sam.view <- sprintf("%s %s", sam.view, file)
  if (!is.null(regions) & isTRUE(oneRegion)) sam.view <- sprintf("%s %s",
                                                                  sam.view,
                                                                  rstring)
  sam.view <- sprintf("%s > %s", sam.view, outname)
  if (!is.null(remove.mitochondrial)) {
    sam.view <- sprintf("%s idxstats %s | cut -f 1 | grep -v %s | xargs %s",
                        samtools, file, shQuote(remove.mitochondrial), sam.view)
  }


  if (isTRUE(count)) {
    res <- vector("list", length(sam.view))
    for (i in seq_along(sam.view)) {
      res[[i]] <- .cmd(cmd = sam.view[[i]], intern = TRUE)
      res[[i]] <- as.numeric(unlist(res[[i]]))
    }
    return(res)
  } else {
    for (i in seq_along(sam.view)) {
      .cmd(cmd = sam.view[[i]])
    }
    return(outname)
  }
}

#' @description \code{run_samflagstat} uses \code{samtools flagstat} to
#' calculate and print statistics from a BAM file.  It provides counts for 13
#' categories of reads: total, secondary, supplementary, duplicates, mapped,
#' paired in sequencing, read1, read2, properly paired, with itself and its mate
#' mapped, singletons, mate mapped to a different chromosome, with mate mapped
#' to a different chromosome with MAPQ>5.
#' @rdname run_samsort
#' @return \code{run_samflagstat} returns a dataframe of alignment statistics.
#' @examples
#' \dontrun{
#' run_samflagstat(samtools = "samtools", threads = parallel::detectCores(),
#'                 bamfile = "HB1_sample.bam")
#' }
#'
#' @export

run_samflagstat <- function(samtools = "samtools",
                            threads = 1,
                            bamfile = NULL) {
  sam.flag <- sprintf("%s flagstat -@ %s %s", samtools, threads, bamfile)
  res <- vector("list", length(sam.flag))
  for (i in seq_along(sam.flag)) {
    res[[i]] <- .cmd(sam.flag[[i]], intern = TRUE)
    res[[i]] <- lapply(res[[i]], .extract_numbers)
    res[[i]] <- lapply(res[[i]], as.integer)
    Flag <- c("Total alignments + unaligned", "Marked secondary",
              "Marked supplementary", "Marked duplicate", "Mapped",
              "Paired in sequencing", "Mate 1",
              "Mate 2", "Aligned as proper pair",
              "Both mates mapped", "Aligned as singleton",
              "Mate mapped to a different chr",
              "Mate mapped to a different chr and MAPQ>5")
    res[[i]] <- data.frame(Flag = Flag,
                           QC.passed = sapply(res[[i]], function(x) x[1]),
                           QC.failed = sapply(res[[i]], function(x) x[2]))
    total <- res[[i]][1, 2] + res[[i]][1, 3]
    res[[i]] <- c(total, res[[i]][,2])
    names(res[[i]]) <- c(Flag[1], "QC-passed", Flag[2:13])
  }
  if (length(res) == 1) {
    res <- as.data.frame(res[[1]])
  } else {
    res <- as.data.frame(do.call(cbind, res))
  }
  names(res) <- bamfile
  return(res)
}

