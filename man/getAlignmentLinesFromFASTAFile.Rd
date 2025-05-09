% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeatMapTools.R
\name{getAlignmentLinesFromFASTAFile}
\alias{getAlignmentLinesFromFASTAFile}
\title{Read Alignment Lines from FASTA File (Internal)}
\usage{
getAlignmentLinesFromFASTAFile(file_str)
}
\arguments{
\item{file_str}{A character string specifying the path to the FASTA file.}
}
\value{
A data frame with columns \code{lineID} and \code{sequence}.
}
\description{
Internal helper that reads a multi-sequence FASTA file and returns a standardized data frame compatible with alignment tools.
}
\keyword{internal}
