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
    #print("is fasta")
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
            #if(length(object@alignmentLines) == 0){s
            #  return(0)
            #}

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

            #return (toString((substr(object@alignment_data$sequence,column,column)),sep="" ))
            return (paste((substr(object@alignment_data$sequence,column,column)),collapse="" ))
            #return (substr(object@alignment_data$sequence,column,column))
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
            #print("Executing get score ")
            iString <- getColumnAsString(object,i)
            jString <- getColumnAsString(object,j)
            # print("istring and jstring")
            # print(iString)
            # print(jString)
            iFrequency <- getResidueFrequencies(iString)
            jFrequency <- getResidueFrequencies(jString)
            # print("frequencies")
            # print(iFrequency)
            # print(jFrequency)
            pairs <- getPairs(iString,jString,iFrequency,jFrequency)
            # print("Result from getPairs")
            # print(pairs)
            if(length(pairs$num) == 1){
              return(0)
              #print("There was only one pair and is perfectly conserved")
            }else{
              # print(paste("how big is this string?",nchar(iString)))
              # print(paste("pairs:",pairs))
              score <- getCovarianceScore(pairs,nchar(iString),iFrequency,jFrequency)
            }
            #score <-0
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

            #print(paste0("max_rows_in_scores:",max_rows_in_scores));
            # scores <- matrix(0, nrow = max_rows_in_scores, ncol = 3);
            # visiting_row <- 1;
            #
            # for(i in 1:(max_col_length)){
            #   for(j in (i+1):max_col_length){
            #     if(j > max_col_length){
            #       break;
            #     }
            #     scores[visiting_row,1] <- i;
            #     scores[visiting_row,2] <- j;
            #     #print("Getting regular score")
            #     scores[visiting_row,3] <- getScore(object,i,j);
            #     visiting_row <- visiting_row +1;
            #   }
            # }
            # #print(paste("how many visited rows?:",visiting_row))
            # #print(paste("how many supposed to visit?:",max_rows_in_scores))
            # colnames(scores)<- c("i","j","score")
            # scores_table <- data.frame(scores)
            # full_scores <-project_geometrical_matrix(scores_table)
            # return(full_scores)

            scores <- expand.grid(i = 1:max_col_length, j = 1:max_col_length)
            #print(scores)
            # Filter for the lower triangular part
            # scores <- scores[scores$i > scores$j, ]
            # # Define a function to calculate scores
            # #print("print about to define a function")
            # calculate_score <- function(i, j) {
            #   #print("used the function parallel")
            #   getScore(object, i, j)
            # }
            # #print("before running mcmapply")
            # num_cores <- detectCores()
            # scores$score <- mcmapply(calculate_score, scores$i, scores$j, mc.cores = num_cores)
            # #print(scores)
            # #scores_table <- data.frame(scores)
            # #full_scores <-project_geometrical_matrix(scores_table)
            # return(scores)
            scores_lower <- scores[scores$i > scores$j, ]

            # Define a function to calculate scores
            calculate_score <- function(i, j) {
              getScore(object, i, j)
            }

            # Set the number of cores to use
            num_cores <- min(cores, parallel::detectCores())
            #print("we are in the getMICovarianceScore normal")
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
