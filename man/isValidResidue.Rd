% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeatMapTools.R
\name{isValidResidue}
\alias{isValidResidue}
\title{Check if Residue is Valid (Internal)}
\usage{
isValidResidue(letter)
}
\arguments{
\item{letter}{A single-character string representing a residue (e.g., \code{"A"}, \code{"W"}).}
}
\value{
\code{TRUE} if the residue is valid, otherwise \code{FALSE}.
}
\description{
Internal helper function that checks whether a character matches one of the standard 20 amino acid residue codes.
}
\keyword{internal}
