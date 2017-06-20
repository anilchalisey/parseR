#' Aggregate and summarise FastQC reports for multiple samples
#'
#' @description Aggregate and summarise multiple FastQC reports.
#'
#' @param object A fqcRes.multi object produced by `read_multifastqc()`.
#'
#' @return
#' \itemize{
#'     \item \strong{aggregate_fqc()} returns a list containing two data-frames.
#'     The first data-frame has the following columns:
#'          \itemize{
#'          \item sample: sample names
#'          \item tot.seq: total sequences (i.e.: the number of reads)
#'          \item low.qual: the number of sequences flagged as poor quality
#'          \item seq.len: sequence length
#'          \item pct.gc: \% of GC content
#'          \item pct.dup: \% of duplicate reads
#'          }
#'      In the second datframe, each row is a sample and the columns are the
#'      FASTQC modules.  Each cell indicates the outcome (pass, warn or fail)
#'      for that module by each sample.
#'    \item \strong{module_fqc()}: Generates a list of summarising the
#'    module results.
#'          \itemize{
#'          \item The first list item is a dataframe, in which each row is
#'    a sample and the columns are the FASTQC modules.  Each cell indicates the
#'    outcome (pass, warn or fail) for that module by each sample.
#'          \item The second item summarises the module results per module,
#'          indicating for each module, how many samples passed, warned or
#'          failed.  In addition, individual summaries for warned and failed are
#'          provided.
#'          \item The third item summrises the module results per sample,
#'          indicating for each sample, how many modules passed, warned or
#'          failed.  In addition, individual summaries for warned and failed are
#'          provided.
#'          }
#'    \item \strong{stats_fqc()}: returns a data frame containing
#'    basic metrics of the fastqc reports.  Columns are: sample, tot.seq,
#'    low.qual.seq, seq.len, pct.gc, and pct.dup.
#'    }
#' @importFrom stats aggregate
#' @importFrom reshape2 melt
#' @importFrom tidyr spread
#' @import dplyr
#'
#' @export

aggregate_fqc <- function(object) {

  if (!inherits(object, "fqcRes_multi"))
    stop("An object of class fqcRes_multi required.")

  metrics <- as.data.frame(do.call(
    rbind, lapply(object, function(x) x$basic_statistics[4:7,2])))
  colnames(metrics) <- c("tot.seq", "low.qual", "seq.len", "pct.gc")
  metrics$`pct.dup` <- unname(sapply(object, function(x) x$pct_duplication))
  metrics$sample <- gsub("_fastqc.zip", "", rownames(metrics))
  metrics <- metrics[, c(6, 1:5)]
  rownames(metrics) <- NULL

  modules <- as.data.frame(do.call(
    rbind, lapply(object, function(x) x$Summary)))
  modules <- reshape2::dcast(modules,
                             formula = sample ~ module, value.var = "status")
  modules$sample <- gsub("_fastqc.gz", "", modules$sample)
  res <- list("metrics" = metrics, "module_summary" = modules)
  return(res)
}

#' @rdname aggregate_fqc
#'
#' @export

module_fqc <- function(object) {
  res <- aggregate_fqc(object)
  ms <- list()
  ms$all <- res$module_summary
  ms$all$sample <- gsub(".fastq.gz", "", ms$all$sample)
  ms$by_module <- list()
  ms$by_module$summary <- ms$all %>%
    reshape2::melt(id.vars = "sample") %>%
    dplyr::group_by(variable, value) %>%
    dplyr::summarise(count = n()) %>%
    tidyr::spread(key = "value", value = "count", fill = 0) %>%
    dplyr::select(Module = variable, PASS, WARN, FAIL) %>%
    as.data.frame()
  ms$by_module$WARN <- .get_problems(ms$all, "module", "WARN")
  ms$by_module$FAIL <- .get_problems(ms$all, "module", "FAIL")
  ms$by_sample <- list()
  ms$by_sample$summary <- ms$all %>%
    reshape2::melt(id.vars = "sample") %>%
    dplyr::group_by(sample, value) %>%
    dplyr::summarise(count = n()) %>%
    tidyr::spread(key = "value", value = "count", fill = 0) %>%
    dplyr::select(sample, PASS, WARN, FAIL) %>%
    as.data.frame()
  ms$by_sample$WARN <- .get_problems(ms$all, "sample", "WARN")
  ms$by_sample$FAIL <- .get_problems(ms$all, "sample", "FAIL")
  return(ms)
}

#' @rdname aggregate_fqc
#'
#' @export

stats_fqc <- function(object) {
  res <- aggregate_fqc(object)
  res$metrics
}


.get_problems <- function(x, by, status) {

  if (by == "module") .formula <- sample~module
  if (by == "sample") .formula <- module~sample
  probs <- x %>%
    reshape2::melt(id.vars = "sample") %>%
    dplyr::filter(value == status) %>%
    dplyr::select(sample, module = variable, status = value) %>%
    stats::aggregate(.formula, data = ., paste, collapse = ", ")
  nb <- x %>%
    reshape2::melt(id.vars = "sample") %>%
    dplyr::filter(value == status) %>%
    dplyr::select(sample, module = variable, status = value) %>%
    dplyr::group_by_(by) %>%
    dplyr::summarise(number = n())
  res <- merge(nb, probs, by = by)
  return(res)
}
