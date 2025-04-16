#' Standard amino acid residue list (internal).
#' @keywords internal
residues <- c('A','C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
#' Residue list length (internal).
#' @keywords internal
residues_length <- 20
#' Draw Mutual Information Heatmap
#'
#' Plots a mutual information (MI) heatmap from pairwise MI scores using `ggplot2`. The heatmap includes annotations for alignment domains and optional reference lines between regions.
#'
#' @param data A data frame with columns \code{i}, \code{j}, and \code{score}, representing mutual information scores between alignment positions.
#' @param alignment_labels A character vector of domain or region labels corresponding to alignment segments.
#' @param alignment_label_locations A numeric vector with x/y coordinates for placing \code{alignment_labels} on the plot.
#' @param alignment_domain_residues A numeric vector indicating the boundaries of alignment domains or regions.
#' @param heatmap_lines Logical. If \code{TRUE}, adds separating lines between regions on the heatmap.
#'
#' @return A \code{ggplot} object representing the MI heatmap.
#' @export
draw_MI_heatmap <- function(data, alignment_labels,alignment_label_locations,alignment_domain_residues,heatmap_lines){
    mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

    if(heatmap_lines == TRUE){
      mycol <- c("white", "navy")
    }
    MIplot <- ggplot2::ggplot(data, ggplot2::aes(i, j, fill=score))

    MIplot <- MIplot +ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = mycol)+
    ggplot2::theme_void()+
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))

    heatmap_max_height = max(alignment_domain_residues)
    MIplot <- MIplot + ggplot2::annotate("segment", x = 0, xend = 0, y = 0, yend = heatmap_max_height,colour = "black")
    MIplot <- MIplot + ggplot2::annotate("segment", y = 0, yend = 0, x = 0, xend = heatmap_max_height,colour = "black")
    for(i in 1:length(alignment_labels)) {
      if(i == 1){
        MIplot <- MIplot + ggplot2::annotate("segment", x = 0, xend = alignment_domain_residues[i], y = -10, yend = -10,colour = "black")
      }else{
        MIplot <- MIplot + ggplot2::annotate("segment", x = alignment_domain_residues[i-1]+1, xend = alignment_domain_residues[i], y = -10, yend = -10,colour = "black")
      }
      MIplot <- MIplot + ggplot2::annotate("text", x = alignment_label_locations[i], y = -20, label = alignment_labels[i])
      MIplot <- MIplot + ggplot2::annotate("text", x = alignment_domain_residues[i], y = -25, label = alignment_domain_residues[i],size=3)

      if(heatmap_lines == TRUE){
        MIplot <- MIplot + ggplot2::annotate("segment", x = alignment_domain_residues[i], xend = alignment_domain_residues[i], y = 0, yend = heatmap_max_height,colour = "black")
      }
    }

    for(i in 1:length(alignment_labels)) {
      if(i == 1){
        MIplot <- MIplot + ggplot2::annotate("segment", x = -10, xend = -10, y = 0, yend = alignment_domain_residues[i],colour = "black")
      }else{
        MIplot <- MIplot + ggplot2::annotate("segment", x = -10, xend = -10, y = alignment_domain_residues[i-1]+1, yend = alignment_domain_residues[i],colour = "black")
      }
      MIplot <- MIplot + ggplot2::annotate("text", x = -20, y = alignment_label_locations[i], label = alignment_labels[i], angle=90)
      MIplot <- MIplot + ggplot2::annotate("text", x = -25,  y = alignment_domain_residues[i],  label = alignment_domain_residues[i],size=3)

      if(heatmap_lines == TRUE){
        MIplot <- MIplot + ggplot2::annotate("segment", y = alignment_domain_residues[i], yend = alignment_domain_residues[i], x = 0, xend = heatmap_max_height,colour = "black")
      }
    }

  return (MIplot)
}
#' Project Lower-Triangle Matrix Internally
#'
#' Used internally to mirror mutual information scores across the diagonal for visualization purposes.
#'
#' @param raw_data A data frame with columns \code{i}, \code{j}, and \code{score}, usually output from \code{getMICovarianceScore}.
#' @return A symmetric data frame with columns \code{i}, \code{j}, and \code{score}.
#' @keywords internal
project_geometrical_matrix <-function(raw_data){
  adjusted_raw_data <- raw_data %>%
    dplyr::mutate(new_j = j-1)%>%
    dplyr::select(i,new_j, score) %>%
    dplyr::rename(j=new_j)

  bottom_left <- adjusted_raw_data %>%
    dplyr::select (i,j,score)

  top_left <- adjusted_raw_data %>%
    dplyr::select (j,i,score) %>%
    dplyr::rename (i=j,j=i)

  total <- rbind(bottom_left,top_left)
  return (total)
}

#' Check if Residue is Valid (Internal)
#'
#' Internal helper function that checks whether a character matches one of the standard 20 amino acid residue codes.
#'
#' @param letter A single-character string representing a residue (e.g., \code{"A"}, \code{"W"}).
#'
#' @return \code{TRUE} if the residue is valid, otherwise \code{FALSE}.
#' @keywords internal
isValidResidue <- function(letter) {
  #residues <- getResidues()
  #residues_length
  for(residue in residues){
    if(letter == residue){
      return (TRUE)
    }
  }
  #for(i in 1:residues_length){
  #  if (letter == residues[i]){
  #    return(TRUE)
  #  }
  #}

  return(FALSE)
}

#' Get Standard Amino Acid Residue List (Internal)
#'
#' Returns the list of amino acid residue codes used internally.
#'
#' @return A character vector of valid residue codes.
#' @keywords internal
getResidues <-function(){
  return (residues)
}

#' Get Residue Letter from Index (Internal)
#'
#' Internal helper function that retrieves a residue letter from the standard amino acid list based on a 1-based index.
#'
#' @param i An integer index from 1 to 20.
#' @return A single-character string representing the residue (e.g., \code{"A"}, \code{"W"}).
#' @keywords internal
getResidueFromList<-function (i){
  #residues <- getResidues()
  return(residues[i])
}

#' Read Aligned FASTA File
#'
#' Reads a multi-sequence FASTA file and returns a data frame containing sequence names and aligned sequences.
#'
#' @param file A character string specifying the path to a FASTA file.
#'
#' @return A data frame with two columns: \code{name} (sequence header without \code{>}) and \code{sequence} (the aligned sequence in uppercase letters).
#'
#' @examples
#' \dontrun{
#' ReadFasta("alignment.fasta")
#' }
#' @export
ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))

  for(i in 1:length(ind)) {
    seqs[i]<-toupper(paste(fasta[s$from[i]:s$to[i]], collapse=""))
  }
  # Create a data frame
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

#' Read Alignment Lines from Plain Text File (Internal)
#'
#' Reads a non-FASTA alignment file (e.g., tab-delimited) into a data frame with two columns: sequence ID and aligned sequence.
#'
#' @param file_str A character string specifying the path to the input file.
#' @return A data frame with columns \code{lineID} and \code{sequence}.
#' @keywords internal
getAlignmentLinesFromFile <- function (file_str){
  MI_raw_alignment <- read.delim(file_str, header = FALSE, col.names = c("lineID","sequence"), sep="")
  return(MI_raw_alignment)
}

#' Read Alignment Lines from FASTA File (Internal)
#'
#' Internal helper that reads a multi-sequence FASTA file and returns a standardized data frame compatible with alignment tools.
#'
#' @param file_str A character string specifying the path to the FASTA file.
#' @return A data frame with columns \code{lineID} and \code{sequence}.
#' @keywords internal
getAlignmentLinesFromFASTAFile <- function (file_str){
  #print("executingGetAlignmentLinesFromFasta")
  #MI_raw_alignment <- read.delim(file_str, header = FALSE, col.names = c("lineID","sequence"), sep="")
  dataframeFASTA <- ReadFasta(file_str)
  colnames(dataframeFASTA)[1] <- "lineID"
  colnames(dataframeFASTA)[2] <- "sequence"

  return(dataframeFASTA)
}

#' Check Residue and Return Index (Internal)
#'
#' Internal helper function that returns the index of a valid residue in the standard amino acid list. Returns 0 if the residue is not found or is a gap.
#'
#' @param letter A single-character string representing an amino acid residue.
#' @return An integer index from 1–20 if the residue is valid; otherwise, 0.
#' @keywords internal
isInResidueList <- function(letter){
  #residues <- getResidues()
  #print(paste("checking letter",letter))
  if (letter == "-"){
    return(0)
  }
  for(i in 1:residues_length){
    if (letter == residues[i]){
      return(i)
    }
  }
  return(0)
}

#' Get Residue Frequencies from a Sequence
#'
#' Computes the frequency of each standard amino acid residue in a given input sequence.
#'
#' @param sequence A character string representing a protein sequence (e.g., \code{"MKLIVT..."}) using one-letter residue codes.
#'
#' @return A numeric matrix (20x1) where each row corresponds to the count of a standard residue in the input sequence, in the same order as returned by \code{getResidues()}.
#'
#' @examples
#' getResidueFrequencies("ACDEFGHIKLMNPQRSTVWY")
#' @export
getResidueFrequencies <-function(sequence){
  #residues <- getResidues()
  #residue_frequencies = numeric(residues_length)
  residue_frequencies = matrix(0, nrow = residues_length, ncol = 1)
  rf_seq_len = nchar(sequence)
  for (i in 1:rf_seq_len){
    found_number  <- isInResidueList(substr(sequence,i,i))
    if(found_number>0){
      residue_frequencies[found_number] <- residue_frequencies[found_number]+1
    }
  }
  return(residue_frequencies)
}

#' Concatenate Two Residues into a Pair String (Internal)
#'
#' Internal helper function that joins two one-letter residue codes into a pair string (e.g., \code{"A"} and \code{"L"} become \code{"AL"}).
#'
#' @param ichar A single-character string representing the first residue.
#' @param jchar A single-character string representing the second residue.
#'
#' @return A two-character string representing the residue pair.
#' @keywords internal
makePair <-function(ichar,jchar){
  paste0(ichar,jchar)
}
#' Generate and Count Residue Pairs (Internal)
#'
#' Constructs all valid residue-residue pairs from two aligned sequence strings and counts their co-occurrence across positions. Used in mutual information and covariance scoring.
#'
#' @param istring A character string representing the first aligned sequence column.
#' @param jstring A character string representing the second aligned sequence column.
#' @param ifreq A numeric vector of residue frequencies in \code{istring}.
#' @param jfreq A numeric vector of residue frequencies in \code{jstring}.
#'
#' @return A data frame with columns \code{pair} (two-letter residue string) and \code{num} (occurrence count).
#' @keywords internal
getPairs <- function(istring,jstring,ifreq,jfreq){
 count <- 0

 len_ifreq = length(ifreq)
 len_jfreq = length(jfreq)

 count = sum(ifreq > 0) * sum (jfreq>0)

 letter_pairs <-  matrix(, nrow = count, ncol = 1)
 letter_counts <-  matrix(0, nrow = count, ncol = 1)

 if(count==0){
  return(letter_pairs)
 }

count <- 0
 for(i in 1:len_ifreq){
   if(ifreq[i]>0){
     for(j in 1:len_jfreq){
       if(jfreq[j]>0){
         count <- count +1;
         pair <- makePair(getResidueFromList(i),getResidueFromList(j))
         letter_pairs[count] <- pair
       }
     }
   }
 }


for(i in 1:nchar(istring)){
  pair<-paste(substr(istring,i,i),substr(jstring,i,i),sep = "")
    for(j in 1:length(letter_pairs)){
      if(pair == letter_pairs[j]){
        letter_counts[j] <- letter_counts[j]+1
        break;
      }
    }
}

return(data.frame(pair=letter_pairs,num=letter_counts, stringsAsFactors=FALSE))
}

#' Generate and Count Residue Pairs (Internal)
#'
#' Constructs all valid residue-residue pairs from two aligned sequence strings and counts their co-occurrence across positions. Used in mutual information and covariance scoring.
#'
#' @param istring A character string representing the first aligned sequence column.
#' @param jstring A character string representing the second aligned sequence column.
#' @param ifreq A numeric vector of residue frequencies in \code{istring}.
#' @param jfreq A numeric vector of residue frequencies in \code{jstring}.
#'
#' @return A data frame with columns \code{pair} (two-letter residue string) and \code{num} (occurrence count).
#' @keywords internal
getPairs3 <- function(istring, jstring, ifreq, jfreq) {
  # Find non-zero indices in ifreq and jfreq
  nonzero_ifreq <- which(ifreq > 0)
  nonzero_jfreq <- which(jfreq > 0)

  # Create all possible pairs
  pairs <- expand.grid(ifreq = residues[nonzero_ifreq], jfreq = residues[nonzero_jfreq])
  pair_map <- setNames(1:length(pairs$ifreq), paste0(pairs$ifreq, pairs$jfreq))

  # Initialize counts to zero
  letter_counts <- rep(0, length(pairs$ifreq))

  # Iterate over pairs in the strings
  for (i in 1:nchar(istring)) {
    pair <- paste0(substr(istring, i, i), substr(jstring, i, i))
    if (pair %in% names(pair_map)) {
      index <- pair_map[[pair]]
      letter_counts[pair_map[[pair]]] <- letter_counts[pair_map[[pair]]] + 1
    }
  }
  pairs$counts <- letter_counts
  return(pairs)
}




#' Compute Covariance Score from Residue Pair Counts (Internal)
#'
#' Calculates a mutual information–based covariance score given observed residue pair counts and their marginal frequencies. Used by internal scoring methods.
#'
#' @param pairs A data frame with columns \code{pair} (residue pair strings) and \code{num} (observed counts).
#' @param numValidSequences The number of aligned sequences used to normalize the score.
#' @param iFrequencies A numeric vector of residue frequencies in the first alignment column.
#' @param jFrequencies A numeric vector of residue frequencies in the second alignment column.
#'
#' @return A numeric value representing the covariance score.
#' @keywords internal
getCovarianceScore <-function(pairs,numValidSequences,iFrequencies,jFrequencies){
  covarianceScore <- 0
  #print("we are in the old getCovarianceScore Function")
  for(ipair in 1:length(pairs$num)){
    pair_str <-pairs[ipair,"pair"]
    pair_num <-pairs[ipair,"num"]
    if(pair_num>0){
    fxi <- (iFrequencies[isInResidueList(substring(pairs[ipair,"pair"], 1, 1))]) / numValidSequences
    fxj <- (jFrequencies[isInResidueList(substring(pairs[ipair,"pair"], 2, 2))]) / numValidSequences
    jointProb <- pair_num/numValidSequences
    covarianceScore <- covarianceScore + (jointProb * ( log( jointProb/ (fxi*fxj) )))
    }
  }
  return(covarianceScore)
}



#' Annotate MI scores by gene regions
#'
#' @param iscores Data frame of MI scores.
#' @param igenes A vector of gene labels (e.g., c("gene1", "gene2", "interface")).
#' @param imidway Integer index for gene boundary.
#' @return A grouped tibble by gene.
#' @importFrom dplyr mutate group_by
getHistWithLabels <- function(iscores,igenes,imidway){

  iscores_hist <- iscores %>%
    dplyr::mutate(gene = ifelse(i < imidway & j < imidway+1,
                         igenes[1],
                         ifelse(i>imidway+1 & j>imidway+1,
                                igenes[2],
                                ifelse(i<imidway & j>imidway+1,
                                       igenes[3],"NA")))) %>%
    dplyr::group_by(gene)
  return(iscores_hist)
}
