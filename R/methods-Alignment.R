#' Construct an Alignment Object
#'
#' This function initializes an S4 \code{Alignment} object from a given file path. The file can be a plain text file or a FASTA file. It reads the alignment data using the appropriate parser based on the file extension.
#'
#' @param alignmentID A character string giving a unique identifier for the alignment. Defaults to \code{"AlignmentID"}.
#' @param alignment_path A character string indicating the path to the alignment file. Accepted formats: \code{.txt} or \code{.fasta}.
#'
#' @return An S4 object of class \code{Alignment} containing the parsed alignment data.
#'
#' @examples
#' \dontrun{
#' # Create an alignment object from a FASTA file
#' aln <- Alignment(alignment_path = "example.fasta")
#'
#' # Create an alignment object from a tab-delimited text file
#' aln <- Alignment(alignment_path = "alignment.txt")
#'}
#' @export
Alignment <-function(alignmentID = "AlignmentID", alignment_path = "cvnAF2_1025_translated_50percent_alignment_nogaps.txt"){
  .Object <- new('Alignment',alignmentID=alignmentID,alignment_path = alignment_path)
  .Object@alignment_path = alignment_path

  if(tolower(sub(pattern = "^(.*\\.|[^.]+)(?=[^.]*)", replacement = "", alignment_path, perl = TRUE)) == "txt"){
    .Object@alignment_data = getAlignmentLinesFromFile(alignment_path)
  }

  if(tolower(sub(pattern = "^(.*\\.|[^.]+)(?=[^.]*)", replacement = "", alignment_path, perl = TRUE)) == "fasta"){
    .Object@alignment_data <- getAlignmentLinesFromFASTAFile(alignment_path)
  }

  if(tolower(sub(pattern = "^(.*\\.|[^.]+)(?=[^.]*)", replacement = "", alignment_path, perl = TRUE)) == "fas"){
    .Object@alignment_data <- getAlignmentLinesFromFASTAFile(alignment_path)
  }

  return(.Object)
}

#' @describeIn assertAllAlignmentsEqualLength
#' Checks whether all sequences in the alignment have equal length. Returns the length if equal, or \code{-1} if they differ.
#' @return An integer: the shared alignment length if consistent, or \code{-1} if inconsistent.
#' @export
setMethod("assertAllAlignmentsEqualLength", signature(object = "Alignment"),
          function(object) {

            alength <- str_length(object@alignmentLines[1])

            for(i in 1:length(object@alignmentLines)) {
              if(alength != str_length(object@alignmentLines[i])){
                return(-1)
              }
              alength <- str_length(object@alignmentLines[i])
            }
            return(alength)
          })


#' @describeIn getColumnAsString
#' Extracts a single column from the alignment and returns it as a concatenated string of residues from all sequences.
#'
#' @param column An integer specifying the column index to extract.
#' @return A character string representing the residues at the specified alignment column.
#' @export
setMethod("getColumnAsString", signature(object = "Alignment"),
          function(object,column) {
            return (paste((substr(object@alignment_data$sequence,column,column)),collapse="" ))
          }
          )

#' @describeIn getScore
#' Computes the covariance score between two columns in an \code{Alignment} object. Frequencies and residue pairings are calculated to produce a score based on variability and co-occurrence.
#'
#' @param i An integer. Index of the first alignment column.
#' @param j An integer. Index of the second alignment column.
#'
#' @return A numeric covariance score between the specified columns.
#' @export
setMethod("getScore", signature(object = "Alignment"),
          function(object,i,j) {
            iString <- getColumnAsString(object,i)
            jString <- getColumnAsString(object,j)

            iFrequency <- getResidueFrequencies(iString)
            jFrequency <- getResidueFrequencies(jString)

            pairs <- getPairs(iString,jString,iFrequency,jFrequency)

            if(length(pairs$num) == 1){
              return(0)
            }else{
              score <- getCovarianceScore(pairs,nchar(iString),iFrequency,jFrequency)
            }
            return(score)
          })


#' @describeIn getMICovarianceScore
#' Calculates mutual information-based covariance scores for all pairs of alignment columns in an \code{Alignment} object. Scores are computed using pairwise comparisons in parallel and returned as a symmetric data frame.
#' @param cores An integer. Number of cores to use. Default = 1.
#' @return A data frame with three columns: \code{i} (first column index), \code{j} (second column index), and \code{score} (the mutual information score for the pair).
#' @importFrom stringr str_length
#' @export
setMethod("getMICovarianceScore", signature(object = "Alignment"),
          function(object, cores = 1) {
            max_col_length <- stringr::str_length(object@alignment_data[1,"sequence"]);
            scoreij <-0;
            line <- ""
            max_rows_in_scores <-0;
            max_rows_in_scores<-(max_col_length*(max_col_length-1))/2

            scores <- expand.grid(i = 1:max_col_length, j = 1:max_col_length)

            scores_lower <- scores[scores$i > scores$j, ]

            # Define a function to calculate scores
            calculate_score <- function(i, j) {
              getScore(object, i, j)
            }

            # Set the number of cores to use
            num_cores <- min(cores, parallel::detectCores())

            # Use mcmapply for parallel computation for the lower triangular part
            scores_lower$score <- mcmapply(calculate_score, scores_lower$i, scores_lower$j, mc.cores = num_cores)

            # Mirror the lower triangular matrix to obtain a symmetric matrix
            scores_upper <- scores_lower
            scores_upper[c("i", "j")] <- scores_lower[c("j", "i")]
            scores_upper <- scores_upper[, c("j", "i", "score")]

            # Combine the upper and lower triangular parts to get the symmetric matrix
            scores <- rbind(scores_lower, scores_upper)
            return(scores)

          })





#' @describeIn getMICovarianceScoreTriangleMatrix
#' Computes mutual information-based covariance scores for all unique pairs of alignment columns (i < j), and returns them as a data frame in triangular format.
#'
#' @return A data frame with columns: \code{i} (first column index), \code{j} (second column index), and \code{score} (mutual information score between columns). Only unique column pairs (upper triangle) are included.
#' @importFrom stringr str_length
#' @export
setMethod("getMICovarianceScoreTriangleMatrix", signature(object = "Alignment"),

          function(object) {
            max_col_length <- stringr::str_length(object@alignment_data[1,"sequence"]);
            scoreij <-0;
            line <- ""
            max_rows_in_scores <-0;
            for(i in 1:max_col_length){
              for(j in (i+1):max_col_length){
                if(j > max_col_length){
                  break;
                }
                max_rows_in_scores <- max_rows_in_scores + 1;
              }
            }
            scores <- matrix(, nrow = max_rows_in_scores, ncol = 3);
            visiting_row <- 1;

            for(i in 1:(max_col_length)){
              for(j in (i+1):max_col_length){
                if(j > max_col_length){
                  break;
                }
                scores[visiting_row,1] <- i;
                scores[visiting_row,2] <- j;
                scores[visiting_row,3] <- getScore(object,i,j);
                visiting_row <- visiting_row +1;
              }
            }
            colnames(scores)<- c("i","j","score")
            scores_table <- data.frame(scores)

            return(scores_table)
          })

#' Get Sequence Names
#'
#' Returns the sequence names (headers) from the alignment.
#'
#' @param object An \code{Alignment} object.
#' @return A character vector of sequence names.
#' @export
get_sequence_names <- function(object) {
  stopifnot(is(object, "Alignment"))
  return(object@alignment_data$lineID)
}
#' Get a Single Sequence by Name
#'
#' Returns the amino acid sequence corresponding to a given name.
#'
#' @param object An \code{Alignment} object.
#' @param name A character string specifying the sequence name.
#' @return A character string of the sequence, or \code{NULL} if not found.
#' @export
get_sequence <- function(object, name) {
  stopifnot(is(object, "Alignment"))
  row <- object@alignment_data[object@alignment_data$lineID == name, ]
  if (nrow(row) == 0) return(NULL)
  return(row$sequence)
}

#' Get Alignment Length
#'
#' Returns the length of the aligned sequences (including gaps).
#'
#' @param object An \code{Alignment} object.
#' @return An integer representing the alignment length, or 0 if empty.
#' @export
get_alignment_length <- function(object) {
  stopifnot(is(object, "Alignment"))
  if (nrow(object@alignment_data) > 0) {
    return(nchar(object@alignment_data$sequence[1]))
  }
  return(0)
}

#' Read MI CSV File
#'
#' Reads a mutual information score file with columns i, j, and score.
#'
#' @param path Path to the CSV file
#' @return A data frame with columns i, j, score
#' @export
readMICSV <- function(path) {
  read.csv(path)
}

#' Attach MI CSV Data to Alignment Object
#'
#' Loads raw MI data into the Alignment object (as a data frame).
#'
#' @param object An \code{Alignment} object.
#' @param path Path to the MI CSV file.
#' @return Updated \code{Alignment} object with \code{mi_data} slot populated.
#' @export
addMICSVData <- function(object, path) {
  stopifnot(is(object, "Alignment"))
  object@mi_data <- read.csv(path)
  return(object)
}

#' Get Top MI Score for a Given Position
#'
#' Returns the row from the MI data frame with the highest score for position \code{i}.
#'
#' @param object An \code{Alignment} object.
#' @param i Integer index of the query position.
#' @return A list with elements \code{j} and \code{score}, or \code{NULL} if not found.
#' @export
getTopMIFor <- function(object, i) {
  stopifnot(is(object, "Alignment"))

  df <- object@mi_data
  hits <- df[df$i == i, ]
  if (nrow(hits) == 0) return(NULL)

  top <- hits[which.max(hits$score), ]
  return(list(j = top$j, score = top$score))
}

#' Grep Sequence Names from Alignment
#'
#' Returns sequence names that match the given regex pattern.
#'
#' @param object An \code{Alignment} object.
#' @param pattern A regular expression pattern (as a character string).
#' @return A character vector of matching sequence names.
#' @export
grepNames <- function(object, pattern) {
  stopifnot(is(object, "Alignment"))
  names <- get_sequence_names(object)
  return(names[grepl(pattern, names)])
}

#' Get Scored Residues for a Sequence
#'
#' Given a sequence name, return a list of tuples:
#' (residue_index, alignment_position, residue, score)
#' Skips gaps and uses 0.0 if no MI score is available.
#'
#' @param object An \code{Alignment} object
#' @param name A sequence name
#' @return A data frame with columns: residue_index, alignment_position, residue, score
#' @export
getScoredResidues <- function(object, name) {
  stopifnot(is(object, "Alignment"))

  seq <- get_sequence(object, name)
  if (is.null(seq)) stop(paste("Sequence", name, "not found in alignment."))

  residue_index <- 0
  results <- list()

  for (i in seq_along(strsplit(seq, "")[[1]])) {
    res <- substr(seq, i, i)
    if (res == "-") next

    residue_index <- residue_index + 1
    score <- object@mi_top_scores[[as.character(i - 1)]]  # Python-style 0-based
    score_val <- if (!is.null(score)) score$score else 0.0

    results[[length(results) + 1]] <- list(
      residue_index = residue_index,
      alignment_position = i,
      residue = res,
      score = score_val
    )
  }

  # Convert to data frame
  do.call(rbind.data.frame, results)
}

#' Get Scored Residues for a Sequence
#'
#' For a given sequence name, return a data frame with:
#' - residue_index (non-gap count)
#' - alignment_position (with gaps)
#' - residue (amino acid)
#' - score (from getTopMIFor; defaults to 0.0 if not found)
#'
#' @param object An \code{Alignment} object.
#' @param name The sequence name to extract.
#' @return A data frame with columns: residue_index, alignment_position, residue, score
#' @export
getScoredResidues <- function(object, name) {
  stopifnot(is(object, "Alignment"))

  seq <- get_sequence(object, name)
  if (is.null(seq)) stop(paste("Sequence", name, "not found in alignment."))

  aa <- strsplit(seq, "")[[1]]
  residue_index <- 0
  results <- list()

  for (i in seq_along(aa)) {
    res <- aa[i]
    if (res == "-") next

    residue_index <- residue_index + 1
    top <- getTopMIFor(object, i - 1)  # 0-based like Python
    score_val <- if (!is.null(top)) top$score else 0.0

    results[[length(results) + 1]] <- list(
      residue_index = residue_index,
      alignment_position = i,
      residue = res,
      score = score_val
    )
  }

  return(do.call(rbind.data.frame, results))
}


#' Read Alignment Column Mask
#'
#' Reads a mask file containing a line like '#ColumnsMap\t586, 587, 588, ...'
#' and returns the list of column indices.
#'
#' @param mask_path Path to the mask file
#' @return An integer vector of column indices
#' @export
readMask <- function(mask_path) {
  lines <- readLines(mask_path)

  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, "#ColumnsMap")) {
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) != 2) stop("Malformed #ColumnsMap line.")
      return(as.integer(trimws(unlist(strsplit(parts[2], ",")))))
    }
  }

  stop("No valid #ColumnsMap line found in mask file.")
}










