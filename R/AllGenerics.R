#' Check That All Sequences Are Equal Length
#'
#' Generic function to check if all aligned sequences in an \code{Alignment} object have the same length.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments (unused).
#' @return An integer: sequence length if equal; otherwise, \code{-1}.
#' @importFrom stringr str_length
#' @export
setGeneric("assertAllAlignmentsEqualLength", function(object, ...) {
  standardGeneric("assertAllAlignmentsEqualLength")
})
#' Extract Column from Alignment as String
#'
#' Generic function to retrieve a specific column from an alignment object as a character string.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments, such as the column index.
#' @return A character string representing one alignment column.
#' @export
setGeneric("getColumnAsString", function(object, ...) {
  standardGeneric("getColumnAsString")
})

#' Get Covariance Score
#'
#' Generic function to compute a covariance score between residues in an alignment.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments (e.g., residue positions).
#' @return A numeric covariance score.
#' @export
setGeneric("getScore", function(object, ...) {
  standardGeneric("getScore")
})

#' Get Covariance Score (Version 3)
#'
#' Alternate version of \code{getScore}, that may include perfectly conserved residues.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments.
#' @return A numeric covariance score.
#' @export
setGeneric("getScore3", function(object, ...) {
  standardGeneric("getScore3")
})


#' Get Mutual Information Score
#'
#' Generic function to compute a mutual information-based covariance score between residue columns.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments (e.g., residue positions).
#' @return A numeric mutual information score.
#' @export
setGeneric("getMICovarianceScore", function(object, ...) {
  standardGeneric("getMICovarianceScore")
})

#' Get MI Covariance Score (Version 3)
#'
#' Alternate version of \code{getMICovarianceScore}, get the Covariance Scores based on getScore3 which may include perfectly conserved residues.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments.
#' @return A numeric mutual information score.
#' @export
setGeneric("getMICovarianceScore3", function(object, ...) {
  standardGeneric("getMICovarianceScore3")
})

#' Get MI Triangle Matrix
#'
#' Generic function to compute a lower-triangular matrix of mutual information scores for all pairwise combinations of alignment positions.
#'
#' @param object An \code{Alignment} object.
#' @param ... Additional arguments.
#' @return A numeric matrix of mutual information scores.
#' @export
setGeneric("getMICovarianceScoreTriangleMatrix", function(object, ...) {
  standardGeneric("getMICovarianceScoreTriangleMatrix")
})
