% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HeatMapTools.R
\name{draw_MI_heatmap}
\alias{draw_MI_heatmap}
\title{Draw Mutual Information Heatmap}
\usage{
draw_MI_heatmap(
  data,
  alignment_labels,
  alignment_label_locations,
  alignment_domain_residues,
  heatmap_lines
)
}
\arguments{
\item{data}{A data frame with columns \code{i}, \code{j}, and \code{score}, representing mutual information scores between alignment positions.}

\item{alignment_labels}{A character vector of domain or region labels corresponding to alignment segments.}

\item{alignment_label_locations}{A numeric vector with x/y coordinates for placing \code{alignment_labels} on the plot.}

\item{alignment_domain_residues}{A numeric vector indicating the boundaries of alignment domains or regions.}

\item{heatmap_lines}{Logical. If \code{TRUE}, adds separating lines between regions on the heatmap.}
}
\value{
A \code{ggplot} object representing the MI heatmap.
}
\description{
Plots a mutual information (MI) heatmap from pairwise MI scores using `ggplot2`. The heatmap includes annotations for alignment domains and optional reference lines between regions.
}
