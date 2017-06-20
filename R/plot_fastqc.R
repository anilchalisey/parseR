#' Plot FastQC Results
#' @description Plot FastQC data
#' @param fqcRes An object of class fqcRes or a FASTQC zipped result file.
#' @param modules Character vector containing the names of fastqc modules for
#'   which you want to import the data. Default is all. Allowed values include
#'   one or the combination of:
#'   \itemize{
#'   \item "summary",
#'   \item "basic_statistics",
#'   \item "per_base_sequence_quality",
#'   \item "per_sequence_quality_scores",
#'   \item "per_base_sequence_content",
#'   \item "per_sequence_gc_content",
#'   \item "per_base_n_content",
#'   \item "sequence_duplication_levels",
#'   \item "overrepresented_sequences",
#'   \item "adapter_content",
#'   \item "kmer_content"
#'   }
#' Partial match of module names are allowed and case is ignored.  Also, white
#' space may be used in place of underscore. For example, you can use
#' modules = "GC content", instead of "per_sequence_gc_content").
#'
#' @return Returns a list of ggplots (+/- dataframes) containing plots of the
#'   specified modules.
#'
#' @import ggplot2
#' @import ggthemes
#' @import grid
#' @importFrom grDevices windowsFont windowsFonts
#'
#' @export

plot_fqc <- function(fqcRes, modules = "all") {

  if (inherits(fqcRes, "character"))
    fqcRes <- read_fastqc(fqcRes)
  if (!inherits(fqcRes, "fqcRes"))
    stop("data should be an object of class fqcRes")

  modules <- .valid_modules(modules)
  res <- list()

  if ("Summary" %in% modules)
    res$Summary <- .plot_msum(fqcRes)
  if ("Basic_Statistics" %in% modules)
    res$Basic_Statistics <- .plot_bstats(fqcRes)
  if ("Per_base_sequence_quality" %in% modules)
    res$Per_base_sequence_quality <- .plot_bq(fqcRes)
  if ("Per_sequence_quality_scores" %in% modules)
    res$Per_sequence_quality_scores <- .plot_sq(fqcRes)
  if ("Per_base_sequence_content" %in% modules)
    res$Per_base_sequence_content <- .plot_sc(fqcRes)
  if ("Per_sequence_GC_content" %in% modules)
    res$Per_sequence_GC_content <- .plot_gc(fqcRes)
  if ("Per_base_N_content" %in% modules)
    res$Per_base_N_content <- .plot_nc(fqcRes)
  if ("Sequence_Duplication_Levels" %in% modules)
    res$Sequence_Duplication_Levels <- .plot_dl(fqcRes)
  if ("Overrepresented_sequences" %in% modules)
    res$Overrepresented_sequences <- .plot_os(fqcRes)
  if ("Adapter_Content" %in% modules)
    res$Adapter_Content <- .plot_ad(fqcRes)
  if ("Kmer_content" %in% modules)
    res$Kmer_content <- .plot_kmer(fqcRes)

  p <- res
  return(p)
}

## Hidden plot functions

# Basic statistics - returns a dataframe of the basic statistics module
.plot_bstats <- function(fqcRes){
  d <- fqcRes$basic_statistics
  as.data.frame(d)
}

# Summary - returns a dataframe of the module summary
.plot_msum <- function(fqcRes){
  d <- fqcRes$Summary
  as.data.frame(d)
}

# Per sequence GC content
.plot_gc <- function(fqcRes){
  status <- .get_status(fqcRes, "GC content")
  d <- fqcRes$per_sequence_gc_content
  if (nrow(d) == 0) return(NULL)
  ggplot(d, aes(d$`GC Content`, d$Count)) +
    geom_line(size = 1) +
    labs(title = "Per sequence GC content",
         x = "Mean GC Content (%)",
         y = "Count",
         caption = paste0("Status: ", status)) +
    .theme_Publication() +
    theme(plot.title = element_text(hjust = 0),
          plot.caption = element_text(color = switch(status, PASS = "blue",
                                                     WARN = "orange",
                                                     FAIL = "red")))
}

# Per base N content
.plot_nc <- function(fqcRes) {
  status <- .get_status(fqcRes, "N content")
  d <- fqcRes$per_base_n_content

  if (class(d$Base) == "character") {
    .final <- as.numeric(d$Base[length(d$Base)])
    d$Base <- c(1:9, seq(10, .final - 1, 2), .final)
  }

  ggplot(d, aes(d$Base, d$`N-Count`, group = 1)) +
    geom_line(size = 1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = "Per base N content",
         x = "Position in read (bp)",
         y = "Frequency (%)",
         subtitle = "N content across all bases",
         caption = paste0("Status: ", status)) +
    .theme_Publication() +
    theme(plot.title = element_text(hjust = 0),
                   plot.caption = element_text(color = switch(
                     status,
                     PASS = "blue",
                     WARN = "orange",
                     FAIL = "red")))
}

# Per base sequence quality
.plot_bq <- function(fqcRes) {
  status <- .get_status(fqcRes, "Per base sequence quality")
  d <- fqcRes$per_base_sequence_quality

  if (class(d$Base) == "character") {
    .final <- as.numeric(d$Base[length(d$Base)])
    d$Base <- c(1:9, seq(10, .final - 1, 2), .final)
  }
  ggplot(data = d, aes(x = d$Base,
                       ymin = d$`10th Percentile`,
                       lower = d$`Lower Quartile`,
                       middle = d$Median,
                       upper = d$`Upper Quartile`,
                       ymax = d$`90th Percentile`)) +
    geom_boxplot(stat = "identity", fill = "gold",
                 colour = "black", lwd = 0.2, fatten = 2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,
                      ymax = 20, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20,
                      ymax = 27, alpha = 0.1, fill = "orange") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 27,
                      ymax = Inf, alpha = 0.1, fill = "green") +
    labs(title = "Per base sequence quality",
                  x = "Position in read (bp)",
                  y = "Median quality scores",
                  subtitle = "Green: high quality zone        Amber: medium quality zone        Red: low quality zone",
                  caption = paste0("Status: ", status)) +
    .theme_Publication() +
    theme(plot.title = element_text(hjust = 0),
                   plot.caption = element_text(color = switch(
                     status,
                     PASS = "blue",
                     WARN = "orange",
                     FAIL = "red")))
}

# Per sequence quality scores
.plot_sq <- function(fqcRes) {
  status <- .get_status(fqcRes, "Per sequence quality scores")
  d <- fqcRes$per_sequence_quality_scores

  ggplot(d, aes(d$Quality, d$Count)) +
    geom_line(size = 1) +
    labs(title = "Per sequence quality scores",
         subtitle = "Quality score distribution over all sequences",
         x = "Mean Sequence Quality (Phred Score)",
         y = "Count",
         caption = paste0("Status: ", status)) +
    .theme_Publication() +
    theme(plot.title = element_text(hjust = 0),
                   plot.caption = element_text(color = switch(
                     status,
                     PASS = "blue",
                     WARN = "orange",
                     FAIL = "red")))
}

# Per base sequence content
.plot_sc <- function(fqcRes) {
  status <- .get_status(fqcRes, "Per base sequence content")
  d <- fqcRes$per_base_sequence_content

  if (class(d$Base) == "character") {
    .final <- as.numeric(d$Base[length(d$Base)])
    d$Base <- c(1:9, seq(10, .final - 1, 2), .final)
  }
  d <- reshape2::melt(d, id.vars = "Base")

  ggplot(d, aes(d$Base, d$value, group = d$variable,
                colour = d$variable)) +
    geom_line(size = 0.8) +
    labs(title = "Per base sequence content",
                  x = "Position in read (bp)",
                  y = "Nucleotide frequency (%)",
                  subtitle = "Sequence content across all bases",
                  caption = paste0("Status: ", status)) +
    .theme_Publication() +
    theme(legend.key = element_rect(size = 5),
                   legend.key.size = unit(1.5, 'lines')) +
    guides(color = guide_legend(title = NULL)) +
    theme(plot.title = element_text(hjust = 0),
                   plot.caption = element_text(color = switch(
                     status,
                     PASS = "blue",
                     WARN = "orange",
                     FAIL = "red")))
}

# Sequence duplication levels
.plot_dl <- function(fqcRes) {
  status <- .get_status(fqcRes, "Sequence Duplication Levels")
  d <- fqcRes$sequence_duplication_levels
  d$`Duplication Level` <- factor(d$`Duplication Level`,
                                  levels = d$`Duplication Level`)
  d <- reshape2::melt(d, id.vars = "Duplication Level")

  ggplot(d, aes(d$`Duplication Level`, d$value)) +
    geom_bar(stat = "identity", position = "dodge", colour = "black",
             aes_string(group = d$variable, fill = d$variable)) +
    labs(title = "Sequence duplication level",
                  x = "Sequence Duplication Level",
                  y = "Percentage (%)",
                  subtitle = paste0("Percentage of distinct reads: ",
                                    fqcRes$pct_duplication, "%"),
                  caption = paste0("Status: ", status)) +
    .theme_Publication() +
    scale_fill_manual(values = c("#D55E00", "#0072B2"),
                               labels = c("Deduplicated     ", "Total     ")) +
    theme(legend.key = element_rect(size = 5),
                   legend.key.size = unit(1.5, 'lines')) +
    guides(fill = guide_legend(title = NULL)) +
    theme(plot.title = element_text(hjust = 0),
                   plot.caption = element_text(color = switch(
                     status,
                     PASS = "blue",
                     WARN = "orange",
                     FAIL = "red")))
}

# Over-represented sequences
.plot_os <- function(fqcRes) {
  status <- .get_status(fqcRes, "Overrepresented sequences")
  d <- fqcRes$overrepresented_sequences

  if (nrow(d) == 0 ) {
    ggplot(d) +
      labs(title = "Overrepresented sequences") +
      annotate("text", x = 0.5, y = 0.5,
                        label = "No overrepresented sequences were identified",
                        size = 6, color = "black", fontface = 2) +
      theme_void() +
      theme(plot.caption = element_text(color = switch(status, PASS = "blue",
                                                       WARN = "orange",
                                                       FAIL = "red")))
  } else {
    as.data.frame(d)
  }

}

# Adapter content
.plot_ad <- function(fqcRes) {
  status <- .get_status(fqcRes, "Adapt")
  d <- fqcRes$adapter_content

  if (class(d$Position) == "character") {
    .final <- as.numeric(d$Position[length(d$Position)])
    d$Position <- c(1:9, seq(10, .final - 1, 2), .final)
  }
  d <- reshape2::melt(d, id.vars = "Position")

  ggplot(d, aes_string(x = "Position",
                       y = "value",
                       group = "variable",
                       colour = "variable")) +
    geom_line(size = 1) +
    labs(title = "Adapter content",
                  x = "Position in read (bp)",
                  y = "Adapter frequency (%)",
                  subtitle = "Sequence content across all bases",
                  caption = paste0("Status: ", status)) +
    .theme_Publication() +
    theme(legend.key = element_rect(size = 5),
                   legend.key.size = unit(1.5, 'lines')) +
    guides(color = guide_legend(title = NULL, ncol = 3)) +
    theme(plot.title = element_text(hjust = 0),
                   plot.caption = element_text(color = switch(
                     status,
                     PASS = "blue",
                     WARN = "orange",
                     FAIL = "red")))
}

# Kmers
.plot_kmer <- function(fqcRes) {
  d <-  fqcRes$kmer_content

  if (nrow(d) == 0) {
    ggplot(d) +
      labs(title = "Kmer content") +
      annotate("text", x = 0.5, y = 0.5,
                        label = "No overrepresented kmers",
                        size = 6, colour = "black", fontface = 2) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0),
                     plot.caption = element_text(color = switch(
                       status,
                       PASS = "blue",
                       WARN = "orange",
                       FAIL = "red")))
  } else {
    as.data.frame(d)
  }
}

