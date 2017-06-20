#' Run HOMER for promoter motif enrichment
#'
#' Function to run the HOMER tool for motif finding.  This function searches
#' motifs enriched in the promoters of the gene list provided.
#'
#' @param genelist Vector of the genes in which to look for enrichment.
#' @param genome The reference genome [default = "human"].
#' @param output_dir The directory in which to output the results.
#' @param len Length of motifs to be found  [default = c(8, 10, 12)].
#' @param bg File of gene IDs to use as background.  By default HOMER uses all other
#'           promoters as the background set.
#' @param start Number of base pairs upstream of the promoter start site to begin
#'              looking [default = -300].
#' @param end Number of base pairs downstream of the promoter start site to stop
#'            looking [default = +50].
#' @param nomask Use unmasked version of the genome [default = FALSE].
#' @param S Number of motifs of each length to find [default = 25].
#' @param mis Number of mismatches allowed [default = 2].
#' @param noconvert Do not convert input files into unigene ids [default = FALSE].
#' @param norevopp Do not search the revese strand for motifs [default = FALSE].
#' @param nomotif Do not search for de novo motifs [default = FALSE].
#' @param basic Do not check de novo motifs dor similarity to known [default = FALSE].
#' @param bits Scale sequence logos by information content [default = FALSE].
#' @param nocheck Do not check for similarity between de novo and known
#'                motifs [default = FALSE].
#' @param noknown Do not check for known motif enrichment [default = FALSE].
#' @param b Use binomial distribution to calculate p-values [default = FALSE;
#'          uses hypergeometric instead].
#' @param nogo Do not search for GO enrichment [default = TRUE].
#' @param humanGO Convert IDs to human for GO analysis [default = TRUE].
#' @param noweight Do not perform GC correction [default = FALSE].
#' @param noredun Do not remove predetermined redundant
#'                promoters/sequences [default = FALSE].
#' @param g Input file is a group file, i.e. 1st column = id,
#'          2nd = 0 or 1 [1=target,0=back] [default = FALSE].
#' @param cpg Use CpG\% instead of GC\% for sequence normalization [default = FALSE].
#' @param rand randomise labels for target and backgound sequences [default = FALSE].
#' @param fdr Numerical indicating number of randomisations to perform for
#'            empirical FDR calculation for de novo discovery.
#'            [default = NULL; no fdr calculation]
#' @param p Number of processors to use [default = 1]
#'
#' @return
#' Output from HOMER, and also an R list object of results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' motif_results <- homer_motifs(de_genes, genome = "human", output_dir = "motifs")
#' }

run_homer <- function(genelist,
                      genome = "human",
                      output_dir,
                      # basic options
                      len = NULL,
                      bg = NULL,
                      start = NULL,
                      end = NULL,
                      nomask = FALSE,
                      S = NULL,
                      mis = NULL,
                      noconvert = FALSE,
                      norevopp = FALSE,
                      nomotif = FALSE,
                      # known motif options
                      basic = FALSE,
                      bits = FALSE,
                      nocheck = FALSE,
                      noknown = FALSE,
                      # advanced options
                      b = FALSE,
                      nogo = TRUE,
                      humanGO = TRUE,
                      noweight = FALSE,
                      noredun = FALSE,
                      g = FALSE,
                      cpg = FALSE,
                      rand = FALSE,
                      fdr = NULL,
                      # other
                      p = NULL) {

  dir.create(output_dir, showWarnings = FALSE)
  write.table(genelist, "genelist.txt", row.names = F, col.names = F,
              quote = F, sep = "\t")

  genelist <- "genelist.txt"

  homer <- sprintf("findMotifs.pl %s %s %s -nofacts",
                   genelist, genome, output_dir)
  if (!is.null(len)) {homer   <- paste(homer, "-len", len)}
  if (!is.null(bg)) {homer    <- paste(homer, "-bg", bg)}
  if (!is.null(start)) {homer <- paste(homer, "-start", start)}
  if (!is.null(end)) {homer   <- paste(homer, "-end", end)}
  if (nomask) {homer          <- paste(homer, "-nomask")}
  if (!is.null(S)) {homer     <- paste(homer, "-S", S)}
  if (!is.null(mis)) {homer   <- paste(homer, "-mis", mis)}
  if (noconvert) {homer       <- paste(homer, "-noconvert")}
  if (norevopp) {homer        <- paste(homer, "-norevopp")}
  if (nomotif) {homer         <- paste(homer, "-nomotif")}
  if (noconvert) {homer       <- paste(homer, "-noconvert")}
  if (basic) {homer           <- paste(homer, "-basic")}
  if (bits) {homer            <- paste(homer, "-bits")}
  if (nocheck) {homer         <- paste(homer, "-nocheck")}
  if (noknown) {homer         <- paste(homer, "-noknown")}
  if (b) {homer               <- paste(homer, "-b")}
  if (nogo) {homer            <- paste(homer, "-nogo")}
  if (humanGO) {homer         <- paste(homer, "-humanGO")}
  if (noweight) {homer        <- paste(homer, "-noweight")}
  if (noredun) {homer         <- paste(homer, "-noredun")}
  if (g) {homer               <- paste(homer, "-g")}
  if (cpg) {homer             <- paste(homer, "-cpg")}
  if (rand) {homer            <- paste(homer, "-rand")}
  if (!is.null(fdr)) {homer   <- paste(homer, "-fdr", fdr)}
  if (!is.null(p)) {homer     <- paste(homer, "-p", p)}

  cat(homer)
  .cmd(homer)
  unlink("genelist.txt")
  res <- list(command = homer)

  if (noknown == FALSE) {
    res$knownmotifs <- homer_known(motif_dir = output_dir)
  }

  if (nomotif == FALSE) {
    tmp              <- homer_denovo(motif_dir = output_dir)
    res$homer_motifs <- tmp$homer_motifs
    res$homer_PWMs   <- tmp$homer_PWMs
  }

  return(res)
}
