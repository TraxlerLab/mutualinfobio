% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeatMapTools.R
\name{getResidueFromList}
\alias{getResidueFromList}
\title{Get Residue Letter from Index (Internal)}
\usage{
getResidueFromList(i)
}
\arguments{
\item{i}{An integer index from 1 to 20.}
}
\value{
A single-character string representing the residue (e.g., \code{"A"}, \code{"W"}).
}
\description{
Internal helper function that retrieves a residue letter from the standard amino acid list based on a 1-based index.
}
\keyword{internal}
