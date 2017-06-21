#' Alignment of RNA-seq reads
#'
#' @param threads A positive integer specifying the number of sorting and
#'   compression threads
#' @param plot Logical.  If TRUE, a PNG summarising the alignment statistics is
#'   produced
#' @param out.dir Character string.  Directory to save results in.  If the
#'   directory does not exist, it will be created
#' @param bigwig Logical indicating whether bigwig files should be created.
#'   This is memory and time-consuming.
#' @param hisat2 The path to hisat2 (if not in executable path).
#' @param samtools The path to samtools (if not in executable path).
#' @param sambamba The path to sambamba (if not in executable path).
#' @param idx The basename of the index for the reference genome. The basename
#'   is the name of any of the index files up to but not including the final
#'   .1.ht2, etc.
#' @param reads.dir Character string.  Directory containing the raw reads.  If
#'   this is specified, then reads1 and reads2 should be set as NULL.
#' @param reads1 Character vector of mate1 reads.  If specified, then reads.dir
#'   must be NULL.
#' @param reads2 Character vector of mate2 reads.  If specified, then reads.dir
#'   must be NULL.  Must be the same length as mate1.  If single-end sequencing,
#'   then should be left as NULL.
#' @param fastq Logical indicating if reads are FASTQ files.
#' @param fasta Logical indicating if reads are FASTA files.
#' @param softClipPenalty Sets the maximum (MX) and minimum (MN) penalties for
#' soft-clipping per base, both integers.  Must be given in the format "MX,MN".
#' @param noSoftClip Logical indicating whether to disallow soft-clipping.
#' @param knownSplice Path to text file containing known splice sites.
#' @param strand Specify strand-specific information.  Default is unstranded.
#' @param tmo Logical indicating whether to report only those reads aligning to
#' known transcripts.
#' @param maxAlign Integer indicating the maximum number of distinct primary
#' alignments to search for each read.
#' @param secondary Logical indicating whether to report secondary alignments.
#' @param nomixed By default, when hisat2 cannot find a concordant or discordant
#' alignment for a pair, it then tries to find alignments for the individual
#' mates. If TRUE, this option disables that behavior.
#' @param nodiscordant By default, hisat2 looks for discordant alignments if it
#' cannot find any concordant alignments.  If true, this option disables that
#' behavior.
#' @param rgid Character string, to which the read group ID is set.
#' @param quiet If TRUE, print nothing except alignments and serious errors.
#' @param non_deterministic When set to TRUE, HISAT2 re-initializes its
#' pseudo-random generator for each read using the current time.#' @param rgid
#' @param memory String specifying maximum memory per thread; suffix K/M/G
#'   recognized.
#' @param remove.mitochondrial Character string.  If set, this will remove reads
#'   mapping to the mitochondrial genome.  The string should match the reference
#'   name for the mitochindrial genome in the alignment file.  Examples include
#'   "ChrM", "M" and "MT".
#' @param hash_table Size of hash table for finding read pairs (default is
#' 262144 reads); will be rounded down to the nearest power of two.  For best
#' performance should be > (average coverage) * (insert size).
#' @param overflow_size Size of the overflow list where reads, thrown out of
#' the hash table, get a second chance to meet their pairs (default is
#' 200000 reads); increasing the size reduces the number of temporary files
#' created.
#' @param io_buffer Controls sizes of the two buffers (in MB) used for reading
#' and writing BAM during the second pass (default is 128).
#' @param scale The number of reads to scale to.  If NULL no scaling performed.
#'
#' @return Raw and filtered BAM files +/- bigwig files
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' \dontrun{
#' align_rna(threads = 2, plot = TRUE, out.dir = ".", bigwig = FALSE,
#'           hisat2 = "hisat2", samtools = "samtools", sambamba = "sambamba",
#'           idx = "../prana/data-raw/index/UCSC.hg19",
#'           reads.dir = "../prana/data-raw/seqFiles", fastq = TRUE)
#' }

align_rna <- function(
  ## Important - general options
  threads = 1,
  plot = TRUE,
  out.dir = ".",
  bigwig = FALSE,
  hisat2 = "hisat2",
  samtools = "samtools",
  sambamba = "sambamba",
  ## Important - align functions
  idx = NULL,
  reads.dir = NULL,
  reads1 = NULL,
  reads2 = NULL,
  ## Can usually be left as deafult - align functions
  fastq = TRUE,
  fasta = FALSE,
  softClipPenalty = NULL,
  noSoftClip = FALSE,
  knownSplice = NULL,
  strand = NULL,
  tmo = FALSE,
  secondary = FALSE,
  maxAlign = NULL,
  nomixed = FALSE,
  nodiscordant = FALSE,
  rgid = NULL,
  quiet = FALSE,
  non_deterministic = TRUE,
  # Can usually be left as default - samtools functions
  memory = "1G",
  remove.mitochondrial = "ChrM",
  # Can usually be left as default - sambamba function
  hash_table = 262144,
  overflow_size = 200000,
  io_buffer = 128,
  # Can usually be left as default - make_bigwig function
  scale = 20000000) {

  if (!is.null(reads.dir) & !(is.null(reads1))) {
    stop("Cannot specify both a reads directory and reads1 - please choose one
         or the other.")
  }

  if (.check_cmd(hisat2) == "Not Found") {
    stop("the hisat2 path is invalid or it is not installed correctly.")
  }

  if (.check_cmd(samtools) == "Not Found") {
    stop("the samtools path is invalid or it is not installed correctly.")
  }

  if (.check_cmd(sambamba) == "Not Found") {
    stop("the sambamba path is invalid or it is not installed correctly.")
  }

  # Get the file names
  if (!is.null(reads.dir)) {
    pattern <- "\\.(fastq|fq|fasta|fa).*"
    mate1   <- list.files(path = paste0(reads.dir, "/"),
                          pattern = paste0("*_1", pattern),
                          full.names = T)
    mate2   <- list.files(path = paste0(reads.dir, "/"),
                          pattern = paste0("*_2", pattern),
                          full.names = T)
    if (length(mate2) == 0) mate2 <- NULL
  } else {
    mate1 <- reads1
    mate2 <- reads2
  }

  if (!is.null(mate2)) {
    if (length(mate1) != length(mate2)) {
      stop("The number of reads1 files and reads2 files do not match.")
    }
  }

  # Align
  sam.files <-
  run_hisat2(hisat2 = hisat2, idx = idx, mate1 = mate1, mate2 = mate2,
             fastq = fastq, fasta = fasta, softClipPenalty = softClipPenalty,
             noSoftClip = noSoftClip, noSplice = FALSE,
             knownSplice = knownSplice, strand = strand, tmo = tmo,
             maxAlign = maxAlign, secondary = secondary, minInsert = NULL,
             maxInsert = NULL, nomixed = nomixed,
             nodiscordant = nodiscordant, threads = threads, rgid = rgid,
             quiet = quiet, non_deterministic = non_deterministic)

  # Convert to BAM
  bam.files <-
  run_samsort(samtools = samtools, file = sam.files, outformat = "BAM",
              threads = threads, memory = memory, sortbyname = FALSE,
              suffix = "", keep = FALSE)

  # Mark duplicates and index
  run_sambambadup(sambamba = "sambamba", bamfile = bam.files,
                  outfile = NULL, remove = FALSE,
                  threads = threads, hash_table = hash_table,
                  overflow_size = overflow_size, io_buffer = io_buffer)

  unlink(bam.files)
  torename <- paste0(.remove_ext(bam.files), "_markdup.bam")
  file.rename(from = torename, to = bam.files)
  file.rename(from = paste0(torename, ".bai"), to = paste0(bam.files, ".bai"))

  # Get some statistics for the BAM files
  diag.stats <-
    run_samflagstat(samtools = samtools, threads = threads, bamfile = bam.files)

  data.m <- t(diag.stats)
  data.m <- apply(data.m, 2, function(x) x/data.m[, 1])
  data.m <- data.frame(cbind(rownames(data.m), data.m))
  names(data.m) <- c("sample", rownames(diag.stats))
  data.m <- data.m[, c(1, 2, 3, 7, 11)]
  data.m <- reshape2::melt(data.m, id.vars = "sample")
  data.m$value <- as.numeric(data.m$value)

  p <-
  ggplot(data.m, aes(sample, value)) +
    geom_bar(aes(fill = variable), colour = "black",
             position = "dodge", stat = "identity") +
    .theme_Publication() +
    labs(title = "Alignment statistics",
         subtitle = "Mapping statistics for all samples evaluated",
         y = "Proportion", x = "") +
    guides(fill = guide_legend(title = NULL)) +
    theme(plot.title = element_text(hjust = 0),
          legend.position = "right", legend.direction = "vertical",
          legend.key = element_rect(size = 5),
                   legend.key.size = unit(1.5, 'lines'))

  if (isTRUE(plot)) {ggsave(plot = p, filename = "alignmentstats.png")}

  # Filter out mitochondrial and improperly aligned reads
  filtered.bam.files <-
  run_samview(samtools = samtools, file = bam.files, regions = NULL,
              chrom.sizes = NULL, include.flag = NULL, exclude.flag = NULL,
              minQual = NULL, outformat = "BAM", outname = NULL,
              include.header = FALSE, count = FALSE, threads = threads,
              subsample = FALSE, keep.paired = TRUE, keep.proper.pair = TRUE,
              remove.unmapped = TRUE, remove.not.primary = TRUE,
              remove.duplicates = FALSE, remove.supplementary.alignment = TRUE,
              remove.mitochondrial = remove.mitochondrial)

  # create index
  run_samindex(samtools = samtools, bamfile = filtered.bam.files,
               threads = threads)

  # Create bigwig files
  if (isTRUE(bigwig)) {
    lapply(filtered.bam.files, function(x) make_bigwig(x, scale = 20000000))
  }

  bamdir <- file.path(out.dir, "bam")
  dir.create(bamdir, recursive = TRUE, showWarnings = FALSE)
  lapply(bam.files, function(x) file.rename(x, file.path(bamdir, x)))
  bam.bai <- paste0(bam.files, ".bai")
  lapply(bam.bai, function(x) file.rename(x, file.path(bamdir, x)))
  lapply(list.files(pattern = "*.log"),
         function(x) file.rename(x, file.path(bamdir, x)))
  filtereddir <- file.path(out.dir, "filteredbam")
  dir.create(filtereddir, recursive = TRUE, showWarnings = FALSE)
  filteredbam <- list.files(pattern = "filtered.bam")
  lapply(filteredbam, function(x) file.rename(x, file.path(filtereddir, x)))

  if (isTRUE(bigwig)) {
    bwdir <- file.path(out.dir, "bw")
    dir.create(bwdir, recursive = TRUE, showWarnings = FALSE)
    lapply(list.files(pattern = "*.bw$"), function(x) {
      file.rename(x, file.path(bwdir, x))
    })
  }
}

