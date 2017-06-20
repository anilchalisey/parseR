#' Mark duplicates in BAM file
#'
#' @description Wrapper script to mark duplicates and optionally remove them in
#'   a BAM file using \code{Sambamba}.
#'
#' @param sambamba Path to Sambamba.
#' @param bamfile Vector of characters specifying path to BAM files.
#' @param outfile Name of output file.  If left as NULL, the suffix _markdup or
#'   _dedup will be appended to the input name to indicate marking only or
#'   removal of duplicates.
#' @param remove Boolean. If TRUE, duplicate reads are removed.
#' @param threads Number of threads to use.
#' @param hash_table Size of hash table for finding read pairs (default is
#' 262144 reads); will be rounded down to the nearest power of two.  For best
#' performance should be > (average coverage) * (insert size).
#' @param overflow_size Size of the overflow list where reads, thrown out of
#' the hash table, get a second chance to meet their pairs (default is
#' 200000 reads); increasing the size reduces the number of temporary files
#' created.
#' @param io_buffer Controls sizes of the two buffers (in MB) used for reading
#' and writing BAM during the second pass (default is 128).
#'
#' @return A BAM file in which duplicate reads have been marked or removed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_sambambadup(sambamba = "sambamba", bamfile = "HB1_sample.bam",
#'                 outfile = "HB1_sample_markdup.bam", remove = FALSE,
#'                 threads = (parallel::detectCores() - 1),
#'                 hash_table = 1000000, overflow_size = 1000000)
#' }

run_sambambadup <- function(sambamba = "sambamba",
                            bamfile = NULL,
                            outfile = NULL,
                            remove = FALSE,
                            threads = 1,
                            hash_table = 262144,
                            overflow_size = 200000,
                            io_buffer = 128) {

  bname <- .remove_ext(bamfile)
  if (is.null(outfile)) {
    if (isTRUE(remove)) {
      outfile <- paste0(bname, "_dedup.bam")
    } else {
      outfile <- paste0(bname, "_markdup.bam")
    }
  }

  hash_table <- as.integer(hash_table)
  overflow_size <- as.integer(overflow_size)
  io_buffer <- as.integer(io_buffer)

  run.sambambadup <- vector("list", length(bamfile))

  for (i in seq_along(bamfile)) {
    run.sambambadup[[i]] <- sprintf("%s markdup --hash-table-size=%s \\
                                    --overflow-list-size=%s \\
                                    --io-buffer-size=%s \\
                                    --nthreads=%s", sambamba, hash_table,
                                    overflow_size, io_buffer, threads)

    if (isTRUE(remove)) run.sambambadup <- sprintf("%s --remove-duplicates",
                                                   run.sambambadup[[i]])
    run.sambambadup[[i]] <- sprintf("%s %s %s", run.sambambadup[[i]],
                                    bamfile[i], outfile[i])
  }

  for (i in seq_along(run.sambambadup)) {
    .cmd(run.sambambadup[[i]])
  }
}

