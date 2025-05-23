% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-Alignment.R
\name{get_sequence}
\alias{get_sequence}
\title{Get a Single Sequence by Name}
\usage{
get_sequence(object, name)
}
\arguments{
\item{object}{An \code{Alignment} object.}

\item{name}{A character string specifying the sequence name.}
}
\value{
A character string of the sequence, or \code{NULL} if not found.
}
\description{
Returns the amino acid sequence corresponding to a given name.
}
