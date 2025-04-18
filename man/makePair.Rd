% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeatMapTools.R
\name{makePair}
\alias{makePair}
\title{Concatenate Two Residues into a Pair String (Internal)}
\usage{
makePair(ichar, jchar)
}
\arguments{
\item{ichar}{A single-character string representing the first residue.}

\item{jchar}{A single-character string representing the second residue.}
}
\value{
A two-character string representing the residue pair.
}
\description{
Internal helper function that joins two one-letter residue codes into a pair string (e.g., \code{"A"} and \code{"L"} become \code{"AL"}).
}
\keyword{internal}
