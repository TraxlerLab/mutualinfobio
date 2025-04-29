#' Class Union: numericORNULL (Internal)
#'
#' Internal class union for a slot that can contain either a numeric value or NULL.
#'
#' @docType class
#' @keywords internal
setClassUnion('numericORNULL', members=c('numeric','NULL'))

#' Alignment Class
#'
#' An S4 class to represent a multiple sequence alignment and store associated metadata and precomputed metrics.
#'
#' @slot alignmentID A character string giving a unique name or label for the alignment.
#' @slot alignment_path A character string with the path to the input alignment file (FASTA or plain text).
#' @slot alignment_data A data frame with columns \code{lineID} and \code{sequence}, parsed from the alignment file.
#' @slot alignmentLines A character vector of aligned sequence strings (raw).
#' @slot numColumns A numeric value representing the number of alignment columns.
#' @slot counts A data frame with residue pair counts used in scoring.
#' @slot frequencies A numeric vector of residue frequencies.
#' @slot totalValid A numeric vector representing the number of valid (non-gap) residues per column.
#' @slot columnStrings A character vector where each element is a string of residues from one alignment column.
#' @slot mi_data A data frame containing mutual information scores (i, j, score).
#' @docType class
#' @export
setClass("Alignment",
                     slots=c(
                       alignmentID = "character",
                       alignment_path ="character",
                       alignment_data = "data.frame",
                       alignmentLines="vector",
                       numColumns="numeric",
                       counts = "data.frame",
                       frequencies = "vector",
                       totalValid = "vector",
                       columnStrings = "vector",
                       mi_data = "data.frame"  # load mi data
                     )
                     )
