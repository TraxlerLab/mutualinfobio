#' Example FASTA Alignment from Weigt et al. (2008)
#'
#' This alignment file comes from the supplementary materials of the open access paper by Skerker et al. (2008),
#' "Rewiring the Specificity of Two-Component Signal Transduction Systems" (Cell, 2008).
#' It contains multiple sequence alignments used in mutual information analysis.
#'
#' The original file is licensed under a Creative Commons Attribution (CC BY) license,
#' and is included here with attribution to the original source.
#'
#' @source Skerker et al. (2008) Cell. \doi{10.1016/j.cell.2008.04.040}
#' @format A FASTA-formatted file located at \code{inst/extdata/cell3925mmc4.fasta}
#' @examples
#' fasta_path <- system.file("extdata", "cell3925mmc4.fasta", package = "mutualinfobio")
#' aln <- Alignment(alignment_path = fasta_path)
#' @name cell3925mmc4.fasta
#' @rdname cell3925mmc4.fasta
#' @keywords datasets
